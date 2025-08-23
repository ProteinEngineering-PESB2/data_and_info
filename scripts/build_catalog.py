#!/usr/bin/env python3
# build_catalog.py
from __future__ import annotations

"""
Build a unified catalog.json from a datasets root.

Features
--------
- Supports both layouts:
  1) datasets/<slug>_v1/ (single-level, version suffix in folder name)
  2) datasets/<slug>/<version>/ (nested)
- Reads enriched metadata.json if present (as produced by your tools),
  verifies/updates sha256 and size (optional), and merges into a unified catalog.
- Infers artifacts if metadata.json misses files; picks a primary sensibly.
- Emits: generated_at, base_url (optional), datasets[] with:
    name, slug, title, version, description, license, license_id, license_url,
    count, size_total, tasks, labels, artifacts[], source, extras (filters/schema/class_stats)
- Optionally include only the latest version per dataset (semver-aware if packaging is available).
- Optionally add absolute "url" for each artifact (base_url + relative path).

CLI
---
python build_catalog.py run \
  --datasets-root ./datasets \
  --output ./catalog.json \
  --base-url https://raw.githubusercontent.com/user/repo/main/datasets \
  --latest-only \
  --include-optional \
  --add-urls \
  --recompute-checksums \
  --log-level INFO
"""

import json
import hashlib
import logging
import re
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import typer

# ---------- Optional unified logging (fallback to stdlib) ----------
try:
    from bioclust.logging import configure_logging, get_logger  # type: ignore
    _HAS_BIOC_LOG = True
except Exception:  # pragma: no cover
    _HAS_BIOC_LOG = False
    def configure_logging(level=logging.INFO, **_):  # type: ignore
        logging.basicConfig(level=level, format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s")
    def get_logger(name: str) -> logging.Logger:  # type: ignore
        return logging.getLogger(name)

log = get_logger("catalogs.build_catalog")

# Optional semver support
try:
    from packaging.version import Version as _V  # type: ignore
    def _vkey(s: str):
        try:
            return _V(s)
        except Exception:
            return s
except Exception:
    def _vkey(s: str):
        # naive fallback: split into numeric and alpha chunks
        parts: List[Any] = []
        for tok in re.split(r"([0-9]+)", s):
            if tok.isdigit():
                parts.append(int(tok))
            elif tok:
                parts.append(tok)
        return tuple(parts)


app = typer.Typer(name="build-catalog", help="Build or update catalog.json from datasets.", no_args_is_help=True)


# ------------------------- Helpers -------------------------

def sha256sum(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def guess_kind_format(p: Path) -> Tuple[str, str]:
    """
    Map filename -> (kind, format). Extend as needed.
    """
    name = p.name.lower()
    ext = "".join(p.suffixes).lower()  # supports .fasta.gz
    if name == "metadata.json":
        return "meta", "json"
    if name.startswith(("metadata.", "table.")):
        if name.endswith(".csv"):
            return "table", "csv"
        if name.endswith(".parquet"):
            return "table", "parquet"
    if "sequence" in name or name.endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz", ".faa", ".faa.gz")):
        return "sequences", "fasta"  # keep logical format (even if gzipped)
    if name.endswith((".npy",)):
        return "embeddings", "npy"
    if name.endswith((".npz",)):
        return "embeddings", "npz"
    if name.endswith((".jsonl",)):
        return "table", "jsonl"
    if name.endswith((".json",)) and name != "metadata.json":
        return "labels", "json"
    if name.endswith((".csv",)):
        return "table", "csv"
    if name.endswith((".parquet",)):
        return "table", "parquet"
    if name.endswith((".pdb",)):
        return "structures", "pdb"
    if name.endswith((".cif", ".mmcif")):
        return "structures", "cif"
    return "artifact", ext.lstrip(".")


def pick_primary(artifacts: List[Dict[str, Any]]) -> int:
    """
    Pick the index of the most likely 'primary' artifact in absence of metadata.json.
    Preference: table.(parquet|csv) > sequences.fasta > anything else.
    """
    for i, a in enumerate(artifacts):
        if a.get("kind") == "table" and a.get("format") in ("parquet", "csv"):
            return i
    for i, a in enumerate(artifacts):
        if a.get("kind") == "sequences":
            return i
    return 0 if artifacts else -1


def _now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _snake(s: str) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_").lower()
    return re.sub(r"_+", "_", s)


def parse_name_version_from_dirname(dirname: str) -> Tuple[str, Optional[str]]:
    """
    Try to extract (name, version) from a folder like 'oxidoreductases_aspergillaceae_v1.0.0'.
    Returns (slug, version or None).
    """
    m = re.match(r"^(?P<slug>.+?)_v(?P<ver>[0-9][A-Za-z0-9\.\-\_]*)$", dirname)
    if m:
        return _snake(m.group("slug")), m.group("ver")
    return _snake(dirname), None


# ------------------------- Data models -------------------------

@dataclass
class Artifact:
    path: str
    sha256: str
    size: int
    optional: bool = False
    kind: Optional[str] = None
    format: Optional[str] = None
    primary: Optional[bool] = None
    url: Optional[str] = None

    def to_catalog(self) -> Dict[str, Any]:
        d = {
            "path": self.path,
            "sha256": self.sha256,
            "size": self.size,
        }
        if self.optional:
            d["optional"] = True
        if self.kind:
            d["kind"] = self.kind
        if self.format:
            d["format"] = self.format
        if self.primary is not None:
            d["primary"] = self.primary
        if self.url:
            d["url"] = self.url
        return d


@dataclass
class DatasetEntry:
    name: str
    slug: str
    version: str
    title: str
    description: str
    license: str
    license_id: Optional[str]
    license_url: Optional[str]
    tasks: List[str]
    labels: List[str]
    count: Optional[int]
    size_total: Optional[int]
    artifacts: List[Artifact]
    source: Optional[List[Dict[str, Any]]] = None
    extras: Optional[Dict[str, Any]] = None  # filters, schema, class_stats, etc.

    def to_catalog(self) -> Dict[str, Any]:
        d: Dict[str, Any] = {
            "name": self.name,
            "slug": self.slug,
            "title": self.title or self.name,
            "version": self.version,
            "license": self.license or "NA",
            "description": self.description or "",
            "tasks": self.tasks or [],
            "labels": self.labels or [],
            "artifacts": [a.to_catalog() for a in self.artifacts],
        }
        if self.license_id:
            d["license_id"] = self.license_id
        if self.license_url:
            d["license_url"] = self.license_url
        if self.count is not None:
            d["count"] = int(self.count)
        if self.size_total is not None:
            d["size_total"] = int(self.size_total)
        if self.source:
            d["source"] = self.source
        if self.extras:
            d["extras"] = self.extras
        return d


# ------------------------- Builders -------------------------

def load_internal_metadata(meta_path: Path) -> Dict[str, Any]:
    return json.loads(meta_path.read_text(encoding="utf-8"))


def enrich_artifacts_from_disk(
    version_dir: Path,
    artifacts_meta: List[Dict[str, Any]],
    include_optional: bool,
    recompute_checksums: bool,
) -> Tuple[List[Artifact], int]:
    """
    Build Artifact objects from metadata.json plus any extra files on disk.
    Returns (artifacts, total_size).
    """
    files = [p for p in version_dir.iterdir() if p.is_file()]
    files_map = {p.name: p for p in files}

    artifacts: List[Artifact] = []
    listed = set()

    # 1) Use metadata.json declared artifacts (if any)
    for a in artifacts_meta or []:
        rel = a.get("path")
        if not rel:
            continue
        p = files_map.get(rel)
        if not p:
            log.warning("Listed artifact not found on disk: %s", version_dir / rel)
            continue
        sha = sha256sum(p) if (recompute_checksums or not a.get("sha256")) else a.get("sha256")
        size = p.stat().st_size if (recompute_checksums or not a.get("size")) else int(a.get("size"))
        artifacts.append(
            Artifact(
                path=rel,
                sha256=str(sha),
                size=int(size),
                optional=bool(a.get("optional", False)),
                kind=a.get("kind"),
                format=a.get("format"),
                primary=a.get("primary"),
            )
        )
        listed.add(rel)

    # 2) Add any unlisted files (optional)
    for p in files:
        if p.name == "metadata.json" or p.name in listed:
            continue
        k, f = guess_kind_format(p)
        artifacts.append(
            Artifact(
                path=p.name,
                sha256=sha256sum(p) if (recompute_checksums or True) else "",  # always compute for unlisted
                size=p.stat().st_size,
                optional=True if include_optional else False,
                kind=k,
                format=f,
                primary=False,
            )
        )

    # 3) Ensure exactly one primary
    primaries = [i for i, a in enumerate(artifacts) if a.primary]
    if len(primaries) == 0 and artifacts:
        idx = pick_primary([asdict(a) for a in artifacts])
        if idx >= 0:
            for i, a in enumerate(artifacts):
                a.primary = (i == idx)
    elif len(primaries) > 1:
        first = primaries[0]
        for i, a in enumerate(artifacts):
            a.primary = (i == first)

    total = sum(int(a.size) for a in artifacts)
    return artifacts, total


def build_entry_from_version_dir(
    ds_name: str,
    slug: str,
    version_dir: Path,
    include_optional: bool,
    recompute_checksums: bool,
) -> Optional[DatasetEntry]:
    """
    Build a DatasetEntry from a version directory.
    """
    meta_json = version_dir / "metadata.json"

    # Defaults if metadata is absent
    name = ds_name
    title = ds_name
    version = version_dir.name
    license_ = "NA"
    license_id = None
    license_url = None
    description = ""
    tasks: List[str] = []
    labels: List[str] = []
    count: Optional[int] = None
    size_total: Optional[int] = None
    source: Optional[List[Dict[str, Any]]] = None
    extras: Optional[Dict[str, Any]] = None
    artifacts_meta: List[Dict[str, Any]] = []

    if meta_json.exists():
        meta = load_internal_metadata(meta_json)
        name = meta.get("name", ds_name)
        slug = meta.get("slug", slug)
        version = meta.get("version", version)
        title = meta.get("title", name)
        license_ = meta.get("license", license_)
        license_id = meta.get("license_id", None)
        license_url = meta.get("license_url", None)
        description = meta.get("description", description)
        # tasks / labels (backward compatibility for single 'task')
        if "tasks" in meta and isinstance(meta["tasks"], list):
            tasks = [str(t) for t in meta["tasks"]]
        elif "task" in meta and meta["task"]:
            tasks = [str(meta["task"])]
        labels = [str(x) for x in (meta.get("labels") or [])]
        count = meta.get("count")
        size_total = meta.get("size_total")
        source = meta.get("source")
        # extras packs useful blocks for downstream UIs
        extras = {
            "filters": meta.get("filters"),
            "schema": meta.get("schema"),
            "class_stats": meta.get("class_stats"),
            "splits": meta.get("splits"),
            "label_sep": meta.get("label_sep"),
            "query": meta.get("query"),
        }
        artifacts_meta = list(meta.get("artifacts") or [])

    # Build artifacts (verify/update checksums and size)
    artifacts, computed_total = enrich_artifacts_from_disk(
        version_dir=version_dir,
        artifacts_meta=artifacts_meta,
        include_optional=include_optional,
        recompute_checksums=recompute_checksums,
    )
    if size_total is None:
        size_total = computed_total

    return DatasetEntry(
        name=name,
        slug=slug,
        version=version,
        title=title,
        description=description,
        license=license_,
        license_id=license_id,
        license_url=license_url,
        tasks=tasks,
        labels=labels,
        count=count,
        size_total=size_total,
        artifacts=artifacts,
        source=source,
        extras=extras,
    )


def scan_datasets(
    datasets_root: Path,
    include_optional: bool,
    recompute_checksums: bool,
) -> List[Tuple[str, DatasetEntry, Path]]:
    """
    Scan datasets root and return a list of (dataset_slug, entry, version_dir_path).
    Supports:
      - single-level:   datasets/<slug>_v1/
      - nested:         datasets/<slug>/<version>/
    """
    out: List[Tuple[str, DatasetEntry, Path]] = []

    for ds_dir in sorted([p for p in datasets_root.iterdir() if p.is_dir()]):
        # Case A: nested structure with versions as subfolders
        subdirs = [p for p in ds_dir.iterdir() if p.is_dir()]
        has_nested = any((d / "metadata.json").exists() for d in subdirs) or any(
            list(d.glob("*.csv")) for d in subdirs
        )
        if has_nested:
            slug = _snake(ds_dir.name)
            name = slug
            for vdir in sorted(subdirs, key=lambda p: _vkey(p.name)):
                entry = build_entry_from_version_dir(
                    ds_name=name,
                    slug=slug,
                    version_dir=vdir,
                    include_optional=include_optional,
                    recompute_checksums=recompute_checksums,
                )
                if entry:
                    out.append((slug, entry, vdir))
            continue

        # Case B: single-level folder with version suffix (or not)
        slug_guess, ver_guess = parse_name_version_from_dirname(ds_dir.name)
        entry = build_entry_from_version_dir(
            ds_name=slug_guess,
            slug=slug_guess,
            version_dir=ds_dir,
            include_optional=include_optional,
            recompute_checksums=recompute_checksums,
        )
        if entry:
            # If metadata had a different version/name, keep that; else use parsed guess
            if entry.version == ds_dir.name and ver_guess:
                entry.version = ver_guess
            out.append((entry.slug or slug_guess, entry, ds_dir))

    return out


def select_latest_per_dataset(entries: List[Tuple[str, DatasetEntry, Path]]) -> List[Tuple[str, DatasetEntry, Path]]:
    by_slug: Dict[str, List[Tuple[str, DatasetEntry, Path]]] = {}
    for slug, entry, path in entries:
        by_slug.setdefault(slug, []).append((slug, entry, path))
    latest: List[Tuple[str, DatasetEntry, Path]] = []
    for slug, items in by_slug.items():
        items_sorted = sorted(items, key=lambda t: _vkey(t[1].version), reverse=True)
        latest.append(items_sorted[0])
    return latest


def add_urls_to_artifacts(
    entries: List[Tuple[str, DatasetEntry, Path]],
    datasets_root: Path,
    base_url: Optional[str],
) -> None:
    """
    Compute artifact.url = base_url + relative_path if base_url is provided.
    The relative path is derived from datasets_root to the artifact file.
    """
    if not base_url:
        return
    base_url = base_url.rstrip("/")
    for _, entry, vdir in entries:
        for art in entry.artifacts:
            ap = (vdir / art.path).resolve()
            try:
                rel = ap.relative_to(datasets_root.resolve()).as_posix()
            except Exception:
                # Fallback: join folder name and artifact name
                rel = f"{vdir.name}/{art.path}"
            art.url = f"{base_url}/{rel}"


# ------------------------- Catalog build/write -------------------------

def build_catalog_dict(
    entries: List[Tuple[str, DatasetEntry, Path]],
    base_url: Optional[str],
) -> Dict[str, Any]:
    datasets = [entry.to_catalog() for _, entry, _ in sorted(entries, key=lambda t: (t[0], _vkey(t[1].version)))]
    cat: Dict[str, Any] = {
        "generated_at": _now_iso(),
        "datasets": datasets,
    }
    if base_url:
        cat["base_url"] = base_url
    return cat


def write_catalog(catalog: Dict[str, Any], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(catalog, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    log.info("[OK] Wrote catalog: %s", output)


# ------------------------- CLI -------------------------

@app.command("run")
def run(
    datasets_root: Path = typer.Option(..., "--datasets-root", "-R", help="Root folder containing datasets."),
    output: Path = typer.Option(..., "--output", "-o", help="Output path for catalog.json."),
    base_url: Optional[str] = typer.Option(None, "--base-url", help="Base URL to prefix artifact urls."),
    merge: bool = typer.Option(False, "--merge", help="Merge with existing catalog.json (preserves other datasets)."),
    latest_only: bool = typer.Option(False, "--latest-only", help="Include only the latest version per dataset."),
    include_optional: bool = typer.Option(True, "--include-optional/--no-include-optional", help="Mark unlisted files as optional."),
    add_urls: bool = typer.Option(True, "--add-urls/--no-add-urls", help="Add absolute 'url' for artifacts using --base-url."),
    recompute_checksums: bool = typer.Option(False, "--recompute-checksums", help="Recompute sha256/size from disk."),
    # logging
    log_level: str = typer.Option("INFO", "--log-level", help="DEBUG|INFO|WARNING|ERROR|CRITICAL"),
    log_dir: Optional[Path] = typer.Option(None, "--log-dir", help="Directory for the log file."),
    log_file: str = typer.Option("bioclust_catalogs.log", "--log-file", help="Log filename."),
    no_rotation: bool = typer.Option(True, "--no-rotation/--rotation", help="Disable daily rotation."),
    dry_run: bool = typer.Option(False, "--dry-run", help="Print catalog to stdout and do not write file."),
) -> None:
    """Build or update catalog.json from a datasets folder."""
    configure_logging(
        level=getattr(logging, log_level.upper(), logging.INFO),
        log_dir=str(log_dir) if log_dir else None,
        use_rotation=not no_rotation,
        filename=log_file,
    )

    if not datasets_root.exists():
        raise typer.BadParameter(f"Datasets root not found: {datasets_root}")

    entries = scan_datasets(datasets_root, include_optional=include_optional, recompute_checksums=recompute_checksums)
    log.info("Found %d dataset version(s).", len(entries))

    if latest_only:
        entries = select_latest_per_dataset(entries)
        log.info("Selected latest versions only -> %d dataset(s).", len(entries))

    if add_urls:
        add_urls_to_artifacts(entries, datasets_root=datasets_root, base_url=base_url)

    catalog = build_catalog_dict(entries, base_url=base_url)

    if merge and output.exists():
        try:
            existing = json.loads(output.read_text(encoding="utf-8"))
        except Exception:
            existing = {"datasets": []}
        # Merge policy: replace entries for (name,version); keep others
        existing_by_key: Dict[Tuple[str, str], Dict[str, Any]] = {}
        for d in existing.get("datasets", []):
            existing_by_key[(d.get("name"), d.get("version"))] = d
        for d in catalog["datasets"]:
            existing_by_key[(d.get("name"), d.get("version"))] = d
        merged = {
            "generated_at": _now_iso(),
            "datasets": [existing_by_key[k] for k in sorted(existing_by_key.keys())],
        }
        if "base_url" in (catalog or {}):
            merged["base_url"] = catalog["base_url"]
        elif "base_url" in (existing or {}):
            merged["base_url"] = existing["base_url"]
        catalog = merged

    if dry_run:
        typer.echo(json.dumps(catalog, indent=2, ensure_ascii=False))
        return

    write_catalog(catalog, output)


if __name__ == "__main__":
    app()
