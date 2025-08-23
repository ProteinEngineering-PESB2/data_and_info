# tools/make_enzyme_function_dataset.py
from __future__ import annotations

"""
Build a UniProtKB-derived enzyme dataset (CSV + FASTA + metadata.json).

Features
--------
- CLI parameters for constructing the UniProt query:
  * Multiple taxonomy ids (OR-combined)
  * Multiple EC filters (OR-combined), e.g., "1.*", "2.7.11.*", "3.1.3.1"
  * Length range [min, max]
  * fragment:true|false
  * reviewed:true|false (optional)
  * Arbitrary positive/negative filters: --include / --exclude (repeated)
- Reliable cursor pagination:
  * First page via query params
  * Subsequent pages follow the `next` URL exactly (no extra params)
- Unified logging via bioclust.logging (fallback to stdlib logging)
- Writes CSV + FASTA and a rich metadata.json (license, source, filters, schema, stats)

Examples
--------
# Transferases in Eurotiomycetidae with length 100–800 aa
python tools/make_enzyme_function_dataset.py run \
  --name transferases_eurotiomycetidae \
  --version 1.0.0 \
  --out-dir ./datasets/transferases_eurotiomycetidae_v1 \
  --taxonomy-id 451871 \
  --ec 2.* \
  --length-min 100 --length-max 800 \
  --no-fragment --reviewed \
  --page-size 500 --max-total 1000 \
  --database-release 2025_03 \
  --retrieved-at "2025-08-23T13:45:00Z"

# Oxidoreductases in Aspergillaceae with extra includes/excludes
python tools/make_enzyme_function_dataset.py run \
  --name oxidoreductases_aspergillaceae \
  --version 1.0.0 \
  --out-dir ./datasets/oxidoreductases_aspergillaceae_v1 \
  --taxonomy-id 1131492 \
  --ec 1.* \
  --length-min 100 --length-max 800 \
  --no-fragment \
  --include 'keyword:Oxidoreductase' \
  --exclude 'fragment:true' \
  --page-size 500 --max-total 2000
"""

import csv
import json
import hashlib
import logging
import re
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests
import typer

# ---------- Optional unified logging (fallback to stdlib) ----------
try:
    from bioclust.logging import configure_logging, get_logger  # type: ignore
    _HAS_BIOC_LOG = True
except Exception:  # pragma: no cover
    _HAS_BIOC_LOG = False

    def configure_logging(level=logging.INFO, **_):  # type: ignore
        logging.basicConfig(
            level=level,
            format="%(asctime)s | %(levelname)-8s | %(name)s | %(message)s",
        )

    def get_logger(name: str) -> logging.Logger:  # type: ignore
        return logging.getLogger(name)


UNIPROT_ENDPOINT = "https://rest.uniprot.org/uniprotkb/search"

# Default table fields to fetch from UniProt (flattened into our schema)
FIELDS = [
    "accession",
    "id",
    "protein_name",
    "gene_primary",
    "organism_name",
    "length",
    "sequence",
    "ec",
    "go_id",
    "protein_existence",
    "xref_pdb",
    "xref_alphafolddb",
]

AA_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")

app = typer.Typer(
    name="make-enzyme-function-dataset",
    help="Query UniProtKB, export CSV/FASTA, and generate enriched metadata.json.",
    no_args_is_help=True,
)

log = get_logger("catalogs.make_dataset")


# ------------------------ Utilities ------------------------

def sha256sum(path: Path) -> str:
    """Compute SHA256 of a file."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _quote_if_needed(term: str) -> str:
    """
    Quote a term if it includes whitespace or quotes that may break the parser.
    """
    if re.search(r'\s|"', term):
        term = term.replace('"', '\\"')
        return f'"{term}"'
    return term


def _build_or_block(field: str, values: List[str]) -> Optional[str]:
    """
    Build an OR-combined block for a field, e.g., (ec:1.* OR ec:2.*).
    """
    vals = [v for v in (values or []) if str(v).strip()]
    if not vals:
        return None
    parts = [f"{field}:{_quote_if_needed(str(v).strip())}" for v in vals]
    return "(" + " OR ".join(parts) + ")" if len(parts) > 1 else parts[0]


def _build_length_block(min_len: Optional[int], max_len: Optional[int]) -> Optional[str]:
    """Build a length range query block."""
    if min_len is None and max_len is None:
        return None
    lo = min_len if min_len is not None else 0
    hi = max_len if max_len is not None else 100000
    return f"length:[{lo} TO {hi}]"
    

def _now_iso() -> str:
    """Current UTC time in ISO 8601 format."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


# ------------------------ UniProt fetch & flatten ------------------------

def fetch_uniprot(
    query: str,
    fields: List[str],
    page_size: int = 500,
    max_total: int = 1000,
    timeout: int = 60,
    max_retries: int = 3,
) -> List[Dict[str, Any]]:
    """
    Execute a UniProt query with reliable cursor pagination.

    Strategy:
    - First page: GET with query params.
    - Next pages: follow the `next` URL (or HTTP Link rel="next") exactly,
      with NO extra params. Some deployments ignore/override `cursor` when
      `query` is re-sent, causing clients to loop on the first page.
    """
    items: List[Dict[str, Any]] = []

    base_url = UNIPROT_ENDPOINT
    params = {
        "query": query,
        "format": "json",
        "fields": ",".join(fields),
        "size": min(int(page_size), 500),  # UniProt caps at 500
    }

    retries = 0
    next_url: Optional[str] = None
    seen_urls: set[str] = set()
    page_idx = 0

    headers = {
        "Accept": "application/json",
        "User-Agent": "bioclust-make-dataset/1.0",
    }

    while True:
        try:
            if next_url:
                url = next_url
                if url in seen_urls:
                    log.warning("Next URL already seen. Breaking to avoid loop: %s", url)
                    break
                seen_urls.add(url)
                r = requests.get(url, timeout=timeout, headers=headers)
            else:
                r = requests.get(base_url, params=params, timeout=timeout, headers=headers)

            if r.status_code == 429 and retries < max_retries:
                wait = 2 ** retries
                log.warning("Rate-limited (429). Backing off %ds...", wait)
                import time as _t
                _t.sleep(wait)
                retries += 1
                continue

            r.raise_for_status()
            data = r.json()
        except Exception as e:
            log.exception("Failed request/page %d: %s", page_idx, e)
            raise

        results = data.get("results", []) or []
        page_idx += 1
        log.debug("Page %d: %d records", page_idx, len(results))

        if not results:
            break

        for rec in results:
            try:
                items.append(flatten_uniprot(rec))
                if len(items) >= max_total:
                    log.debug("Reached max_total=%d. Stopping.", max_total)
                    return items[:max_total]
            except Exception as fe:
                log.debug("Skipping record due to flatten failure: %s", fe)
                continue

        # Prefer JSON 'next'; fallback to HTTP Link header
        next_url = data.get("next")
        if not next_url:
            link = r.links.get("next")
            next_url = link["url"] if link and "url" in link else None

        log.debug("Accumulated: %d | next: %s", len(items), next_url or "<none>")

        if not next_url:
            break

        retries = 0  # reset backoff once a page is successful

    return items


def flatten_uniprot(rec: Dict[str, Any]) -> Dict[str, Optional[str]]:
    """
    Flatten a UniProt record into a normalized row.
    """
    flat: Dict[str, Optional[str]] = {}

    # Accession & ID
    flat["accession"] = rec.get("primaryAccession") or rec.get("accession")
    flat["id"] = rec.get("uniProtkbId") or rec.get("id")

    # Protein name (recommended or alternatives)
    prot_desc = rec.get("proteinDescription")
    prot_name = None
    if isinstance(prot_desc, dict):
        rn = prot_desc.get("recommendedName") or {}
        full = rn.get("fullName")
        if isinstance(full, dict):
            prot_name = full.get("value")
        elif isinstance(full, str):
            prot_name = full
        if not prot_name:
            for alt in prot_desc.get("alternativeNames") or []:
                full = alt.get("fullName")
                if isinstance(full, dict) and full.get("value"):
                    prot_name = full["value"]
                    break
                if isinstance(full, str) and full:
                    prot_name = full
                    break
    else:
        prot_name = rec.get("proteinName") or rec.get("protein_name")
    flat["protein_name"] = prot_name

    # Gene primary
    gene = None
    genes = rec.get("genes")
    if isinstance(genes, list) and genes:
        gn = genes[0].get("geneName", {})
        gene = gn.get("value") if isinstance(gn, dict) else gn or None
    flat["gene_primary"] = gene

    # Organism
    org = rec.get("organism")
    flat["organism_name"] = (org.get("scientificName") if isinstance(org, dict) else None) or rec.get("organism_name")

    # Sequence & length
    seq_obj = rec.get("sequence")
    if isinstance(seq_obj, dict):
        flat["length"] = seq_obj.get("length")
        seq_val = seq_obj.get("value")
        flat["sequence"] = seq_val.replace(" ", "").replace("\n", "") if isinstance(seq_val, str) else None
    else:
        flat["length"] = rec.get("length")
        seq_val = rec.get("sequence")
        flat["sequence"] = seq_val.replace(" ", "").replace("\n", "") if isinstance(seq_val, str) else None

    # EC numbers (recommended + alternatives + flat forms)
    ecs: List[str] = []
    if isinstance(prot_desc, dict):
        rn = prot_desc.get("recommendedName") or {}
        for e in rn.get("ecNumbers") or []:
            v = e.get("value") if isinstance(e, dict) else e
            if v:
                ecs.append(str(v))
        for alt in prot_desc.get("alternativeNames") or []:
            for e in alt.get("ecNumbers") or []:
                v = e.get("value") if isinstance(e, dict) else e
                if v:
                    ecs.append(str(v))

    ec_flat = rec.get("ec") or rec.get("ecNumbers")
    if isinstance(ec_flat, str):
        ecs.extend([x.strip() for x in ec_flat.split(";") if x.strip()])
    elif isinstance(ec_flat, list):
        ecs.extend([str(x) for x in ec_flat])

    flat["ec"] = ";".join(sorted(set(ecs))) if ecs else None

    # Cross-refs: GO, PDB, AlphaFold
    go_terms, pdb_ids, afdb_ids = [], [], []
    xrefs = rec.get("uniProtKBCrossReferences")
    if isinstance(xrefs, list):
        for x in xrefs:
            db = x.get("database")
            if db == "GO":
                go_terms.append(x.get("id"))
            elif db == "PDB":
                pdb_ids.append(x.get("id"))
            elif db == "AlphaFoldDB":
                afdb_ids.append(x.get("id"))

    if not go_terms and isinstance(rec.get("go_id"), str):
        go_terms = [t for t in rec["go_id"].split(";") if t]
    if not pdb_ids and isinstance(rec.get("xref_pdb"), str):
        pdb_ids = [t for t in rec["xref_pdb"].split(";") if t]
    if not afdb_ids and isinstance(rec.get("xref_alphafolddb"), str):
        afdb_ids = [t for t in rec["xref_alphafolddb"].split(";") if t]

    flat["go_id"] = ";".join(go_terms) if go_terms else None
    flat["xref_pdb"] = ";".join(pdb_ids) if pdb_ids else None
    flat["xref_alphafolddb"] = ";".join(afdb_ids) if afdb_ids else None

    # Protein existence (PE)
    pe = rec.get("proteinExistence")
    if isinstance(pe, dict):
        flat["protein_existence"] = pe.get("evidenceCode") or pe.get("value")
    else:
        flat["protein_existence"] = pe or rec.get("protein_existence")

    return flat


def write_csv(items: List[Dict[str, Any]], path: Path) -> None:
    """Write normalized items to CSV with fixed header order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS)
        w.writeheader()
        for it in items:
            w.writerow({k: it.get(k) for k in FIELDS})


def write_fasta(items: List[Dict[str, Any]], path: Path) -> None:
    """Write sequences to FASTA; header: >{accession} {organism}."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for it in items:
            acc = it.get("accession")
            seq = (it.get("sequence") or "").strip()
            org = (it.get("organism_name") or "").strip()
            if not acc or not seq:
                continue
            f.write(f">{acc} {org}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")


# ------------------------ Stats & schema ------------------------

def _class_stats(items: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Compute quick stats used in metadata."""
    n = len(items)
    lens = [int(it.get("length") or 0) for it in items if it.get("length")]
    out: Dict[str, Any] = {"rows": n}
    if lens:
        out["length_avg"] = float(round(sum(lens) / max(len(lens), 1), 6))
        out["length_min"] = int(min(lens))
        out["length_max"] = int(max(lens))

    for c in ("ec", "go_id", "xref_pdb", "xref_alphafolddb"):
        nn = sum(1 for it in items if isinstance(it.get(c), str) and it.get(c).strip())
        out[f"{c}_nonnull"] = int(nn)
    return out


def _build_schema() -> Dict[str, Any]:
    """Return a basic Table Schema-like dictionary."""
    return {
        "primaryKey": "accession",
        "fields": [
            {"name": "accession", "type": "string", "constraints": {"required": True, "unique": True}},
            {"name": "id", "type": "string"},
            {"name": "protein_name", "type": "string"},
            {"name": "gene_primary", "type": "string"},
            {"name": "organism_name", "type": "string"},
            {"name": "length", "type": "integer", "constraints": {"minimum": 1}},
            {"name": "sequence", "type": "string", "constraints": {"pattern": "^[ACDEFGHIKLMNPQRSTVWY]+$"}},
            {"name": "ec", "type": "string", "description": "Semicolon-separated EC numbers"},
            {"name": "go_id", "type": "string", "description": "Semicolon-separated GO IDs"},
            {"name": "protein_existence", "type": "string"},
            {"name": "xref_pdb", "type": "string"},
            {"name": "xref_alphafolddb", "type": "string"},
        ],
    }


# ------------------------ Query builder ------------------------

@dataclass
class QueryParams:
    taxonomy_ids: List[int] = field(default_factory=list)
    ec_filters: List[str] = field(default_factory=list)
    length_min: Optional[int] = None
    length_max: Optional[int] = None
    fragment: Optional[bool] = None
    reviewed: Optional[bool] = None
    include_terms: List[str] = field(default_factory=list)
    exclude_terms: List[str] = field(default_factory=list)


def build_query(q: QueryParams) -> str:
    """
    Build a UniProtKB query string from structured parameters.

    Rules:
      - taxonomy_ids and ec_filters are OR-combined within the group.
      - All groups are AND-combined.
      - Excludes are added as NOT (term).
    """
    segs: List[str] = []

    if q.taxonomy_ids:
        segs.append(_build_or_block("taxonomy_id", [str(x) for x in q.taxonomy_ids]))

    if q.ec_filters:
        segs.append(_build_or_block("ec", q.ec_filters))

    rng = _build_length_block(q.length_min, q.length_max)
    if rng:
        segs.append(rng)

    if q.fragment is not None:
        segs.append(f"fragment:{str(q.fragment).lower()}")

    if q.reviewed is not None:
        segs.append(f"reviewed:{str(q.reviewed).lower()}")

    # Free-form positive includes (AND)
    for inc in q.include_terms:
        inc = inc.strip()
        if not inc:
            continue
        if ":" in inc:
            field, val = inc.split(":", 1)
            segs.append(f"{field.strip()}:{_quote_if_needed(val.strip())}")
        else:
            segs.append(_quote_if_needed(inc))

    # Negative filters (NOT ...)
    for exc in q.exclude_terms:
        exc = exc.strip()
        if not exc:
            continue
        if ":" in exc:
            field, val = exc.split(":", 1)
            segs.append(f"NOT ({field.strip()}:{_quote_if_needed(val.strip())})")
        else:
            segs.append(f"NOT ({_quote_if_needed(exc)})")

    segs = [s for s in segs if s]
    return " AND ".join(segs) if segs else "*"


# ------------------------ CLI commands ------------------------

@app.command("count")
def count(
    # Query parameters (same as run)
    taxonomy_id: List[int] = typer.Option(None, "--taxonomy-id", help="Taxonomy ID(s). Repeat to OR-combine."),
    ec: List[str] = typer.Option(None, "--ec", help="EC filter(s), e.g., 1.*  or  2.7.11.*. Repeat to OR-combine."),
    length_min: Optional[int] = typer.Option(None, "--length-min", help="Minimum sequence length."),
    length_max: Optional[int] = typer.Option(None, "--length-max", help="Maximum sequence length."),
    fragment: Optional[bool] = typer.Option(None, "--fragment/--no-fragment", help="Include fragments or not."),
    reviewed: Optional[bool] = typer.Option(None, "--reviewed/--no-reviewed", help="Restrict to reviewed entries."),
    include: List[str] = typer.Option(None, "--include", help="Additional positive filter(s), AND-combined."),
    exclude: List[str] = typer.Option(None, "--exclude", help="Negative filter(s) (NOT ...)."),
    # Fetch control
    page_size: int = typer.Option(500, "--page-size", min=1, max=500, help="UniProt page size."),
    max_total: int = typer.Option(10000, "--max-total", min=1, help="Hard cap on total results to traverse."),
    timeout: int = typer.Option(60, "--timeout", help="HTTP timeout seconds."),
    # Logging
    log_level: str = typer.Option("INFO", "--log-level", help="DEBUG|INFO|WARNING|ERROR|CRITICAL"),
    log_dir: Optional[Path] = typer.Option(None, "--log-dir", help="Directory for the log file."),
    log_file: str = typer.Option("bioclust_catalogs.log", "--log-file", help="Log filename."),
    no_rotation: bool = typer.Option(True, "--no-rotation/--rotation", help="Disable daily rotation for logs."),
) -> None:
    """
    Count results for a query without writing artifacts (fetches only 'accession').
    """
    configure_logging(
        level=getattr(logging, log_level.upper(), logging.INFO),
        log_dir=str(log_dir) if log_dir else None,
        use_rotation=not no_rotation,
        filename=log_file,
    )
    qp = QueryParams(
        taxonomy_ids=list(taxonomy_id or []),
        ec_filters=[s.strip() for s in (ec or []) if s and s.strip()],
        length_min=length_min,
        length_max=length_max,
        fragment=fragment,
        reviewed=reviewed,
        include_terms=[s for s in (include or []) if s.strip()],
        exclude_terms=[s for s in (exclude or []) if s.strip()],
    )
    query = build_query(qp)
    log.info("Counting with query: %s", query)

    # Fetch with minimal payload: only accession
    items = fetch_uniprot(
        query=query,
        fields=["accession"],
        page_size=page_size,
        max_total=max_total,
        timeout=timeout,
    )
    typer.echo(f"COUNT: {len(items)}")


@app.command("run")
def run(
    # Core dataset identity
    name: str = typer.Option(..., "--name", "-n", help="Dataset canonical name (slug-friendly)."),
    version: str = typer.Option("1.0.0", "--version", "-v", help="Semantic version for the dataset."),
    description: Optional[str] = typer.Option(None, "--description", "-d", help="Human-readable description."),
    out_dir: Path = typer.Option(..., "--out-dir", "-o", help="Output directory for CSV/FASTA/metadata.json"),
    # Query parameters
    taxonomy_id: List[int] = typer.Option(None, "--taxonomy-id", help="Taxonomy ID(s). Repeat to OR-combine."),
    ec: List[str] = typer.Option(None, "--ec", help="EC filter(s), e.g., 1.*  or  2.7.11.*. Repeat to OR-combine."),
    length_min: Optional[int] = typer.Option(100, "--length-min", help="Minimum sequence length."),
    length_max: Optional[int] = typer.Option(800, "--length-max", help="Maximum sequence length."),
    fragment: Optional[bool] = typer.Option(False, "--fragment/--no-fragment", help="Include fragments (default: no-fragment)."),
    reviewed: Optional[bool] = typer.Option(None, "--reviewed/--no-reviewed", help="Restrict to reviewed entries (optional)."),
    include: List[str] = typer.Option(None, "--include", help="Additional positive filter(s), AND-combined. Repeat as needed."),
    exclude: List[str] = typer.Option(None, "--exclude", help="Negative filter(s) (NOT ...). Repeat as needed."),
    # Fetch control
    page_size: int = typer.Option(500, "--page-size", min=1, max=500, help="UniProt page size."),
    max_total: int = typer.Option(1000, "--max-total", min=1, help="Maximum total records to fetch."),
    timeout: int = typer.Option(60, "--timeout", help="HTTP timeout seconds."),
    # Metadata extras
    database_release: Optional[str] = typer.Option(None, "--database-release", help="Source DB release label (e.g., 2025_03)."),
    retrieved_at: Optional[str] = typer.Option(None, "--retrieved-at", help="ISO 8601 timestamp; default is now (UTC)."),
    license_id: str = typer.Option("LicenseRef-UniProt-Terms", "--license-id", help="SPDX/custom license identifier."),
    license_url: str = typer.Option("https://www.uniprot.org/help/license", "--license-url", help="License URL."),
    # Output formats / filenames
    csv_name: Optional[str] = typer.Option(None, "--csv-name", help="Override CSV filename; default: <name>.csv"),
    fasta_name: Optional[str] = typer.Option(None, "--fasta-name", help="Override FASTA filename; default: <name>.fasta"),
    # Debug & logging
    dry_run: bool = typer.Option(False, "--dry-run", help="Print the query and exit without fetching."),
    count_only: bool = typer.Option(False, "--count-only", help="Only count results; do not write artifacts."),
    log_level: str = typer.Option("INFO", "--log-level", help="DEBUG|INFO|WARNING|ERROR|CRITICAL"),
    log_dir: Optional[Path] = typer.Option(None, "--log-dir", help="Directory for the log file."),
    log_file: str = typer.Option("bioclust_catalogs.log", "--log-file", help="Log filename."),
    no_rotation: bool = typer.Option(True, "--no-rotation/--rotation", help="Disable daily rotation for logs."),
) -> None:
    """
    Run the UniProt query, export CSV/FASTA, and write metadata.json.
    """
    # Configure logging
    configure_logging(
        level=getattr(logging, log_level.upper(), logging.INFO),
        log_dir=str(log_dir) if log_dir else None,
        use_rotation=not no_rotation,
        filename=log_file,
    )

    # Build query params
    qp = QueryParams(
        taxonomy_ids=list(taxonomy_id or []),
        ec_filters=[s.strip() for s in (ec or []) if s and s.strip()],
        length_min=length_min,
        length_max=length_max,
        fragment=fragment,
        reviewed=reviewed,
        include_terms=[s for s in (include or []) if s.strip()],
        exclude_terms=[s for s in (exclude or []) if s.strip()],
    )
    query = build_query(qp)

    # Dry-run mode to preview the query
    log.info("UniProt query: %s", query)
    if dry_run:
        typer.echo(f"[DRY-RUN] Query: {query}")
        raise typer.Exit(code=0)

    # Fetch
    log.info("Fetching up to %d entries (page size %d)...", max_total, page_size)
    items = fetch_uniprot(
        query=query,
        fields=FIELDS,
        page_size=page_size,
        max_total=max_total,
        timeout=timeout,
    )
    log.info("Fetched %d entries.", len(items))

    if count_only:
        typer.echo(f"COUNT: {len(items)}")
        raise typer.Exit(0)

    # Basic validation (non-fatal)
    accessions = [it.get("accession") for it in items if it.get("accession")]
    if len(accessions) != len(set(accessions)):
        log.warning("Duplicate accessions detected in response.")
    bad_seq = sum(1 for it in items if not isinstance(it.get("sequence"), str) or not AA_PATTERN.match(it["sequence"]))
    if bad_seq > 0:
        log.warning("Found %d entries with non-canonical or missing sequences.", bad_seq)

    # Prepare outputs
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / (csv_name or f"{name}.csv")
    fasta_path = out_dir / (fasta_name or f"{name}.fasta")

    # Export
    write_csv(items, csv_path)
    write_fasta(items, fasta_path)

    # Metadata
    retrieved_iso = retrieved_at or _now_iso()
    meta: Dict[str, Any] = {
        "name": name,
        "slug": name,
        "version": version,
        "description": description or f"UniProtKB-derived dataset for '{name}'.",
        "task": "multi_label_classification",
        "labels": ["ec"],
        "label_sep": ";",
        "query": query,
        "filters": {
            "taxonomy_id": qp.taxonomy_ids or None,
            "ec_list": qp.ec_filters or None,
            "fragment": qp.fragment,
            "reviewed": qp.reviewed,
            "length_range": [qp.length_min or 0, qp.length_max or 100000],
            "include": qp.include_terms or None,
            "exclude": qp.exclude_terms or None,
        },
        "count": len(items),
        "class_stats": _class_stats(items),
        "license": "UniProt Terms of Use",
        "license_id": license_id,
        "license_url": license_url,
        "source": [
            {
                "database": "UniProtKB",
                "endpoint": UNIPROT_ENDPOINT,
                "database_release": database_release,
                "retrieved_at": retrieved_iso,
            }
        ],
        "schema": _build_schema(),
        "artifacts": [
            {
                "kind": "table",
                "format": "csv",
                "path": csv_path.name,
                "primary": True,
                "size": csv_path.stat().st_size,
                "sha256": sha256sum(csv_path),
                "media_type": "text/csv",
            },
            {
                "kind": "sequences",
                "format": "fasta",
                "path": fasta_path.name,
                "size": fasta_path.stat().st_size,
                "sha256": sha256sum(fasta_path),
                "media_type": "chemical/seq+fasta",
            },
        ],
    }
    meta["size_total"] = sum(int(a.get("size", 0)) for a in meta["artifacts"])

    (out_dir / "metadata.json").write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    log.info("Wrote artifacts:\n- %s\n- %s\n- %s", csv_path, fasta_path, out_dir / "metadata.json")
    typer.echo(f"Done. Dataset saved under: {out_dir}")


if __name__ == "__main__":
    app()
