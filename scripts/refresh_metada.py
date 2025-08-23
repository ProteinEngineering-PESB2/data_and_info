# tools/enrich_catalog_metadata.py
from __future__ import annotations

"""
Enrich and normalize dataset metadata.json files under a datasets root.

This tool:
- Ensures a consistent schema for catalogs (fields, license, source info).
- (Optionally) infers artifacts and computes sha256/size.
- Adds useful stats (class_stats), sizes, and a basic table schema.
- Validates table invariants (row count, unique primary key, AA-only sequences).
- Writes updated metadata.json (with optional backup) or runs in dry-run mode.

Conventions
-----------
- Expects per-dataset folders like: datasets/<slug>_v1/ containing:
    - <table>.csv|.tsv
    - <sequences>.fasta (optional)
    - metadata.json (optional; created if missing)
- Primary key column is "accession"; sequence column is "sequence".
- Label column defaults to "ec" (multi-label with ';' separator).

Notes
-----
- No network calls; all information is derived from local files and CLI flags.
- Logging uses project-wide bioclust logging if available.
"""

import hashlib
import json
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import pandas as pd
import typer

# --- Optional unified logging (falls back to stdlib if unavailable) ---
try:
    from bioclust.logging import configure_logging, get_logger  # type: ignore
    _HAS_BIOC_LOG = True
except Exception:  # pragma: no cover
    _HAS_BIOC_LOG = False
    def configure_logging(level=logging.INFO, **_):  # type: ignore
        logging.basicConfig(level=level, format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s")
    def get_logger(name: str) -> logging.Logger:  # type: ignore
        return logging.getLogger(name)

log = get_logger("catalogs.enrich")


app = typer.Typer(
    name="enrich-catalog-metadata",
    help="Normalize and enrich metadata.json files for ML-ready datasets.",
    no_args_is_help=True,
)


# --------------------------- helpers ---------------------------

AA_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


@dataclass
class EnrichOptions:
    datasets_root: Path
    license_id: str = "LicenseRef-UniProt-Terms"
    license_url: str = "https://www.uniprot.org/help/license"
    database_release: Optional[str] = None
    retrieved_at: Optional[str] = None  # ISO 8601 string
    infer_artifacts: bool = True
    recompute_artifacts: bool = False
    backup: bool = True
    dry_run: bool = False
    primary_key: str = "accession"
    sequence_col: str = "sequence"
    label_cols: List[str] = field(default_factory=lambda: ["ec"])
    label_sep: str = ";"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _snake(s: str) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_").lower()
    return re.sub(r"_+", "_", s)


def _detect_table_file(d: Path) -> Optional[Path]:
    for suf in (".csv", ".tsv", ".tab"):
        cands = list(d.glob(f"*{suf}"))
        if cands:
            return cands[0]
    return None


def _detect_fasta_file(d: Path) -> Optional[Path]:
    for suf in (".fasta", ".fa", ".faa", ".fst"):
        cands = list(d.glob(f"*{suf}"))
        if cands:
            return cands[0]
    return None


def _read_metadata(path: Path) -> Dict[str, Any]:
    if not path.exists():
        return {}
    return json.loads(path.read_text(encoding="utf-8"))


def _write_metadata(path: Path, data: Dict[str, Any], backup: bool) -> None:
    if backup and path.exists():
        bak = path.with_suffix(".json.bak")
        path.replace(bak)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")


def _class_stats(df: pd.DataFrame, seq_col: str, label_cols: List[str]) -> Dict[str, Any]:
    stats: Dict[str, Any] = {}
    # Basic counts
    stats["rows"] = int(len(df))
    if "length" in df.columns:
        lens = pd.to_numeric(df["length"], errors="coerce").dropna().astype(int)
    else:
        lens = df[seq_col].astype(str).str.len()
    if len(lens) > 0:
        stats["length_avg"] = float(round(lens.mean(), 6))
        stats["length_min"] = int(lens.min())
        stats["length_max"] = int(lens.max())

    # Label columns non-null counts
    for c in label_cols:
        if c in df.columns:
            nonnull = df[c].astype(str).str.strip().replace({"": None}).notna().sum()
            stats[f"{c}_nonnull"] = int(nonnull)

    # Common external refs
    for c in ("go_id", "xref_pdb", "xref_alphafolddb"):
        if c in df.columns:
            nonnull = df[c].astype(str).str.strip().replace({"": None}).notna().sum()
            stats[f"{c}_nonnull"] = int(nonnull)
    return stats


def _infer_label_sep(df: pd.DataFrame, label_col: str, default: str) -> str:
    if label_col not in df.columns:
        return default
    sample = df[label_col].dropna().astype(str).head(50)
    if any(";" in s for s in sample):
        return ";"
    if any("," in s for s in sample):
        return ","
    return default


def _build_schema(primary_key: str, sequence_col: str) -> Dict[str, Any]:
    return {
        "primaryKey": primary_key,
        "fields": [
            {"name": "accession", "type": "string", "constraints": {"required": True, "unique": True}},
            {"name": "id", "type": "string"},
            {"name": "protein_name", "type": "string"},
            {"name": "gene_primary", "type": "string"},
            {"name": "organism_name", "type": "string"},
            {"name": "length", "type": "integer", "constraints": {"minimum": 1}},
            {"name": sequence_col, "type": "string", "constraints": {"pattern": "^[ACDEFGHIKLMNPQRSTVWY]+$"}},
            {"name": "ec", "type": "string", "description": "Semicolon-separated EC numbers"},
            {"name": "go_id", "type": "string", "description": "Semicolon-separated GO IDs"},
            {"name": "protein_existence", "type": "string"},
            {"name": "xref_pdb", "type": "string"},
            {"name": "xref_alphafolddb", "type": "string"},
        ],
    }


def _infer_filters_from_query(query: Optional[str]) -> Optional[Dict[str, Any]]:
    if not query:
        return None
    filters: Dict[str, Any] = {}
    # Heuristic extraction
    m = re.search(r"taxonomy_id\s*:\s*(\d+)", query, flags=re.I)
    if m:
        filters["taxonomy_id"] = int(m.group(1))
    m = re.search(r"ec\s*:\s*([0-9]+\.\*|\d+\.\d+\.\d+\.\d+|\d+\.\d+\.\*)", query, flags=re.I)
    if m:
        filters["ec_prefix"] = m.group(1)
    m = re.search(r"fragment\s*:\s*(true|false)", query, flags=re.I)
    if m:
        filters["fragment"] = m.group(1).lower() == "true"
    m = re.search(r"length\s*:\s*\[(\d+)\s+TO\s+(\d+)\]", query, flags=re.I)
    if m:
        filters["length_range"] = [int(m.group(1)), int(m.group(2))]
    return filters or None


def _validate_df(df: pd.DataFrame, opts: EnrichOptions) -> List[str]:
    problems: List[str] = []
    # Primary key existence & uniqueness
    if opts.primary_key not in df.columns:
        problems.append(f"Missing primary key column '{opts.primary_key}'.")
    else:
        pk = df[opts.primary_key].astype(str)
        if pk.isna().any() or (pk == "").any():
            problems.append("Primary key contains empty values.")
        if pk.duplicated().any():
            problems.append("Primary key contains duplicates.")
    # Sequence AA-only check (best-effort)
    if opts.sequence_col in df.columns:
        bad = df[opts.sequence_col].astype(str).str.match(AA_PATTERN).fillna(False)
        if (~bad).any():
            problems.append("Sequence column contains non-canonical residues.")
    return problems


def _load_table(path: Path) -> pd.DataFrame:
    suf = path.suffix.lower()
    if suf in {".csv"}:
        return pd.read_csv(path)
    if suf in {".tsv", ".tab"}:
        return pd.read_csv(path, sep="\t")
    raise ValueError(f"Unsupported table format: {suf}")


def _ensure_artifacts(meta: Dict[str, Any], folder: Path, opts: EnrichOptions) -> Tuple[Dict[str, Any], int]:
    """
    Ensure 'artifacts' section exists and contains sha256/size.
    Returns (updated_meta, total_size_bytes).
    """
    artifacts: List[Dict[str, Any]] = list(meta.get("artifacts") or [])
    by_path = {a.get("path"): a for a in artifacts if "path" in a}

    table_path = _detect_table_file(folder)
    fasta_path = _detect_fasta_file(folder)
    total = 0

    if opts.infer_artifacts:
        # Ensure table artifact exists
        if table_path is not None and (table_path.name not in by_path):
            artifacts.append({
                "kind": "table",
                "format": table_path.suffix.lstrip(".").lower(),
                "path": table_path.name,
                "primary": True,
                "media_type": "text/csv" if table_path.suffix.lower() == ".csv" else "text/tab-separated-values",
            })
        # Ensure fasta artifact exists if file present
        if fasta_path is not None and (fasta_path.name not in by_path):
            artifacts.append({
                "kind": "sequences",
                "format": fasta_path.suffix.lstrip(".").lower(),
                "path": fasta_path.name,
                "media_type": "chemical/seq+fasta",
            })

    # Recompute sha256/size as needed
    for art in artifacts:
        ap = folder / art["path"]
        if not ap.exists():
            log.warning("Artifact path not found: %s", ap)
            continue
        if opts.recompute_artifacts or ("sha256" not in art):
            art["sha256"] = _sha256(ap)
        if opts.recompute_artifacts or ("size" not in art):
            art["size"] = ap.stat().st_size
        total += int(art.get("size", 0))

    meta["artifacts"] = artifacts
    meta["size_total"] = total
    return meta, total


def _enrich_one(folder: Path, opts: EnrichOptions) -> Tuple[Dict[str, Any], Dict[str, Any], Optional[Path]]:
    """
    Enrich metadata.json inside 'folder'. Returns (old_meta, new_meta, meta_path).
    If no table is found, the folder is skipped (returns old==new=={}).
    """
    meta_path = folder / "metadata.json"
    old = _read_metadata(meta_path)

    table_file = _detect_table_file(folder)
    if table_file is None:
        log.info("Skipping (no CSV/TSV found): %s", folder)
        return {}, {}, None

    # Load table
    df = _load_table(table_file)

    # Base fields
    name = old.get("name") or _snake(folder.name.split("_v")[0])
    slug = old.get("slug") or name
    version = old.get("version") or "1.0.0"
    description = old.get("description") or f"Dataset '{name}'"
    query = old.get("query")

    # Labels & separator
    label_cols = old.get("labels") or opts.label_cols
    label_sep = old.get("label_sep") or _infer_label_sep(df, label_cols[0], opts.label_sep)

    # Compute/update artifacts
    new = dict(old)  # start from old
    new["name"] = name
    new["slug"] = slug
    new["version"] = version
    new["description"] = description
    if query:
        new["query"] = query
        filters = _infer_filters_from_query(query)
        if filters:
            new["filters"] = filters

    # Counts & stats
    new["count"] = int(len(df))
    new["class_stats"] = _class_stats(df, opts.sequence_col, label_cols)

    # License fields
    new["license"] = old.get("license") or "UniProt Terms of Use"
    new["license_id"] = old.get("license_id") or opts.license_id
    new["license_url"] = old.get("license_url") or opts.license_url

    # Source info
    source = list(old.get("source") or [])
    if not source:
        source = [{
            "database": "UniProtKB",
            "endpoint": "https://rest.uniprot.org/uniprotkb/search",
        }]
    # Enrich first source entry
    if opts.database_release is not None:
        source[0]["database_release"] = opts.database_release
    if opts.retrieved_at is not None:
        source[0]["retrieved_at"] = opts.retrieved_at
    new["source"] = source

    # Artifacts
    new, _ = _ensure_artifacts(new, folder, opts)

    # Schema (basic)
    new["schema"] = old.get("schema") or _build_schema(opts.primary_key, opts.sequence_col)

    # Task/labels
    new["task"] = old.get("task") or "multi_label_classification"
    new["labels"] = label_cols
    new["label_sep"] = label_sep
    new["splits"] = old.get("splits") or "none"

    # Validate table invariants (non-fatal; reported in 'validation')
    new["validation"] = _validate_df(df, opts)

    return old, new, meta_path


def _diff_dict(old: Dict[str, Any], new: Dict[str, Any]) -> List[str]:
    """Human-readable line diff for console (minimal)."""
    lines: List[str] = []
    old_keys = set(old.keys())
    new_keys = set(new.keys())
    for k in sorted(old_keys | new_keys):
        ov = old.get(k, "<MISSING>")
        nv = new.get(k, "<MISSING>")
        if ov != nv:
            lines.append(f"* {k}: {ov!r} -> {nv!r}")
    return lines


# --------------------------- CLI commands ---------------------------

@app.command("run")
def run(
    datasets_root: Path = typer.Option(..., "--datasets-root", "-R", help="Root directory containing dataset folders."),
    database_release: Optional[str] = typer.Option(None, "--database-release", help="Source database release tag (e.g., 2025_03)."),
    retrieved_at: Optional[str] = typer.Option(None, "--retrieved-at", help="ISO 8601 timestamp for data retrieval."),
    license_id: str = typer.Option("LicenseRef-UniProt-Terms", "--license-id", help="SPDX or custom license identifier."),
    license_url: str = typer.Option("https://www.uniprot.org/help/license", "--license-url", help="License URL."),
    infer_artifacts: bool = typer.Option(True, "--infer-artifacts/--no-infer-artifacts", help="Infer artifacts when missing."),
    recompute_artifacts: bool = typer.Option(False, "--recompute-artifacts/--no-recompute-artifacts", help="Recompute sha256/size for all artifacts."),
    backup: bool = typer.Option(True, "--backup/--no-backup", help="Write metadata.json.bak before overwriting."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Do not write; show differences only."),
    # logging
    log_level: str = typer.Option("INFO", "--log-level", help="DEBUG|INFO|WARNING|ERROR|CRITICAL"),
    log_dir: Optional[Path] = typer.Option(None, "--log-dir", help="Directory for the log file."),
    log_file: str = typer.Option("bioclust_catalogs.log", "--log-file", help="Log filename."),
    no_rotation: bool = typer.Option(True, "--no-rotation/--rotation", help="Disable daily rotation."),
) -> None:
    """Enrich all dataset metadata.json files under --datasets-root."""
    # Configure logging first
    configure_logging(
        level=getattr(logging, log_level.upper(), logging.INFO),
        log_dir=str(log_dir) if log_dir else None,
        use_rotation=not no_rotation,
        filename=log_file,
    )

    opts = EnrichOptions(
        datasets_root=datasets_root,
        license_id=license_id,
        license_url=license_url,
        database_release=database_release,
        retrieved_at=retrieved_at,
        infer_artifacts=infer_artifacts,
        recompute_artifacts=recompute_artifacts,
        backup=backup,
        dry_run=dry_run,
    )

    if not datasets_root.exists():
        raise typer.BadParameter(f"Datasets root not found: {datasets_root}")

    folders = [p for p in datasets_root.iterdir() if p.is_dir()]
    if not folders:
        log.warning("No dataset folders found under: %s", datasets_root)

    changed = 0
    for folder in sorted(folders):
        try:
            old, new, meta_path = _enrich_one(folder, opts)
            if meta_path is None:
                continue
            if old == new:
                log.info("No changes: %s", meta_path)
                continue

            diff = _diff_dict(old, new)
            if dry_run:
                log.info("DRY-RUN changes for %s:\n  %s", meta_path, "\n  ".join(diff) or "(no diff)")
            else:
                _write_metadata(meta_path, new, backup=backup)
                log.info("Updated %s (%d changes).", meta_path, len(diff))
                changed += 1

        except Exception as e:
            log.exception("Failed to process folder %s: %s", folder, e)

    msg = f"Done. Updated {changed} metadata file(s)."
    log.info(msg)
    typer.echo(msg)


if __name__ == "__main__":
    app()
