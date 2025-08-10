#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


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
    if name.startswith(("metadata.", "table.")):
        if name.endswith(".csv"):
            return "table", "csv"
        if name.endswith(".parquet"):
            return "table", "parquet"
    if "sequence" in name or name.endswith((".fa", ".fasta", ".fa.gz", ".fasta.gz")):
        return "sequences", "fasta"  # we keep logical format (even if gzipped)
    if name.endswith((".npy",)):
        return "embeddings", "npy"
    if name.endswith((".npz",)):
        return "embeddings", "npz"
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
    # fallback
    return "artifact", ext.lstrip(".")


def pick_primary(artifacts: List[Dict]) -> int:
    """
    Pick the index of the most likely 'primary' artifact in absence of metadata.json.
    Preference: table.(parquet|csv) > sequences.fasta > anything else.
    """
    for i, a in enumerate(artifacts):
        if a["kind"] == "table" and a["format"] in ("parquet", "csv"):
            return i
    for i, a in enumerate(artifacts):
        if a["kind"] == "sequences":
            return i
    return 0 if artifacts else -1


# ------------------------- Data models -------------------------

@dataclass
class Artifact:
    path: str
    sha256: str
    size: int
    optional: bool = False
    # Optional extras from your internal metadata.json (not required by catalog spec)
    kind: Optional[str] = None
    format: Optional[str] = None
    primary: Optional[bool] = None

    def to_catalog(self) -> Dict:
        d = {
            "path": self.path,
            "sha256": self.sha256,
            "size": self.size,
        }
        if self.optional:
            d["optional"] = True
        # Only include kind/format/primary if you run an extended catalog;
        # Otherwise comment these out to keep the minimal spec.
        if self.kind:
            d["kind"] = self.kind
        if self.format:
            d["format"] = self.format
        if self.primary is not None:
            d["primary"] = self.primary
        return d


@dataclass
class DatasetEntry:
    name: str
    title: str
    version: str
    license: str
    description: str
    tasks: List[str]
    artifacts: List[Artifact]
    homepage: Optional[str] = None
    extras: Optional[Dict[str, object]] = None

    def to_catalog(self) -> Dict:
        return {
            "name": self.name,
            "title": self.title or self.name,
            "version": self.version,
            "license": self.license or "NA",
            "description": self.description or "",
            "tasks": self.tasks or [],
            "artifacts": [a.to_catalog() for a in self.artifacts],
            **({"homepage": self.homepage} if self.homepage else {}),
            **({"extras": self.extras} if self.extras else {}),
        }


# ------------------------- Builders -------------------------

def load_internal_metadata(meta_path: Path) -> Dict:
    """
    Load your per-version metadata.json if present.
    Expected (flexible) structure example:
    {
      "name": "...",
      "version": "...",
      "description": "...",
      "license": "...",
      "tasks": [...],
      "homepage": "...",
      "extras": {...},
      "artifacts": [
        {"path": "metadata.csv", "kind": "table", "format": "csv", "primary": true},
        {"path": "sequences.fasta.gz", "kind": "sequences", "format": "fasta", "optional": true}
      ]
    }
    """
    return json.loads(meta_path.read_text(encoding="utf-8"))


def build_entry_from_version_dir(
    ds_name: str,
    version_dir: Path,
    include_optional: bool,
    fail_on_missing: bool
) -> DatasetEntry:
    meta_json = version_dir / "metadata.json"
    files = [p for p in version_dir.iterdir() if p.is_file() and not p.name.endswith(".meta.json")]
    files_map = {p.name: p for p in files}

    if meta_json.exists():
        meta = load_internal_metadata(meta_json)

        name = meta.get("name", ds_name)
        version = meta.get("version", version_dir.name)
        title = meta.get("title", name)
        license_ = meta.get("license", "NA")
        description = meta.get("description", "")
        tasks = meta.get("tasks", [])
        homepage = meta.get("homepage")
        extras = meta.get("extras")

        artifacts: List[Artifact] = []
        for a in meta.get("artifacts", []):
            rel = a.get("path")
            if not rel:
                continue
            p = files_map.get(rel)
            if not p:
                if fail_on_missing:
                    raise FileNotFoundError(f"Listed artifact not found on disk: {version_dir / rel}")
                else:
                    # skip silently if not found
                    continue
            artifacts.append(
                Artifact(
                    path=rel,
                    sha256=sha256sum(p),
                    size=p.stat().st_size,
                    optional=bool(a.get("optional", False)),
                    kind=a.get("kind"),
                    format=a.get("format"),
                    primary=a.get("primary"),
                )
            )

        # If metadata.json forgot to list some files, we can choose to add them too:
        listed = set(a.path for a in artifacts)
        for p in files:
            if p.name == "metadata.json" or p.name in listed:
                continue
            k, f = guess_kind_format(p)
            artifacts.append(
                Artifact(
                    path=p.name,
                    sha256=sha256sum(p),
                    size=p.stat().st_size,
                    optional=True,  # treat unlisted files as optional
                    kind=k,
                    format=f,
                    primary=False,
                )
            )

        # Ensure exactly one primary (if any were marked)
        primaries = [i for i, a in enumerate(artifacts) if a.primary]
        if len(primaries) == 0 and artifacts:
            # pick a reasonable primary if none set
            idx = pick_primary([asdict(a) for a in artifacts])
            if idx >= 0:
                for i, a in enumerate(artifacts):
                    a.primary = (i == idx)
        elif len(primaries) > 1:
            # demote extras, keep the first as primary
            first = primaries[0]
            for i, a in enumerate(artifacts):
                a.primary = (i == first)

        return DatasetEntry(
            name=name,
            title=title,
            version=version,
            license=license_,
            description=description,
            tasks=tasks,
            artifacts=artifacts,
            homepage=homepage,
            extras=extras,
        )

    # --------- No metadata.json: infer everything ----------
    artifacts: List[Artifact] = []
    for p in files:
        if p.name == "metadata.json":
            continue
        k, f = guess_kind_format(p)
        artifacts.append(
            Artifact(
                path=p.name,
                sha256=sha256sum(p),
                size=p.stat().st_size,
                optional=False,  # set later based on primary
                kind=k,
                format=f,
                primary=False,
            )
        )

    if artifacts:
        idx = pick_primary([asdict(a) for a in artifacts])
        if idx >= 0:
            for i, a in enumerate(artifacts):
                a.primary = (i == idx)
                if i != idx and include_optional:
                    a.optional = True

    return DatasetEntry(
        name=ds_name,
        title=ds_name,
        version=version_dir.name,
        license="NA",
        description="",
        tasks=[],
        artifacts=artifacts,
        homepage=None,
        extras=None,
    )


def scan_datasets(
    datasets_root: Path,
    include_optional: bool,
    fail_on_missing: bool
) -> List[DatasetEntry]:
    entries: List[DatasetEntry] = []
    for ds_dir in sorted([p for p in datasets_root.iterdir() if p.is_dir()]):
        # each subfolder = dataset name
        name = ds_dir.name
        # each version is a subfolder
        versions = [p for p in ds_dir.iterdir() if p.is_dir()]
        if not versions:
            continue
        for vdir in sorted(versions, key=lambda p: p.name):
            entry = build_entry_from_version_dir(
                ds_name=name,
                version_dir=vdir,
                include_optional=include_optional,
                fail_on_missing=fail_on_missing
            )
            entries.append(entry)
    return entries


# ------------------------- Catalog merge/write -------------------------

def load_existing_catalog(path: Path) -> Dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {"datasets": []}


def merge_catalog(existing: Dict, new_entries: List[DatasetEntry], base_url: Optional[str]) -> Dict:
    """
    Merge strategy:
      - Keep datasets not present on disk (from existing).
      - Replace or append entries for (name, version) pairs found on disk.
      - Preserve top-level base_url if provided via arg; else keep existing.
    """
    by_key = {}
    for d in existing.get("datasets", []):
        key = (d.get("name"), d.get("version"))
        by_key[key] = d

    for e in new_entries:
        by_key[(e.name, e.version)] = e.to_catalog()

    datasets = [by_key[k] for k in sorted(by_key.keys())]
    out = {"datasets": datasets}
    if base_url:
        out["base_url"] = base_url
    elif "base_url" in existing:
        out["base_url"] = existing["base_url"]
    return out


def build_catalog(
    datasets_root: Path,
    output: Optional[Path],
    base_url: Optional[str],
    merge: bool,
    include_optional: bool,
    fail_on_missing: bool,
    dry_run: bool
) -> Dict:
    new_entries = scan_datasets(datasets_root, include_optional, fail_on_missing)
    if merge and output and output.exists():
        existing = load_existing_catalog(output)
    else:
        existing = {"datasets": []}
        if output and output.exists() and not merge:
            # read to preserve base_url if user forgets to pass it again
            try:
                existing = load_existing_catalog(output)
            except Exception:
                pass
    catalog = merge_catalog(existing, new_entries, base_url)
    if dry_run:
        print(json.dumps(catalog, indent=2))
    elif output:
        output.write_text(json.dumps(catalog, indent=2), encoding="utf-8")
        print(f"[OK] Wrote catalog: {output}")
    return catalog


# ------------------------- CLI -------------------------

def parse_args(argv: List[str]) -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Build or update catalog.json from a datasets folder.")
    ap.add_argument("--datasets-root", required=True, type=Path, help="Root folder containing datasets/<name>/<version> directories.")
    ap.add_argument("--output", required=True, type=Path, help="Output path for catalog.json.")
    ap.add_argument("--base-url", default=None, help="Base URL prefix for artifact paths (e.g., GitHub raw URL).")
    ap.add_argument("--merge", action="store_true", help="Merge with an existing catalog.json if present.")
    ap.add_argument("--include-optional", action="store_true", help="When inferring (no metadata.json), mark non-primary artifacts as optional.")
    ap.add_argument("--fail-on-missing", action="store_true", help="Fail if metadata.json lists missing files.")
    ap.add_argument("--dry-run", action="store_true", help="Print catalog to stdout instead of writing the file.")
    return ap.parse_args(argv)


def main(argv: List[str]) -> None:
    args = parse_args(argv)
    root = args.datasets_root
    if not root.exists():
        raise SystemExit(f"Datasets root not found: {root}")
    build_catalog(
        datasets_root=root,
        output=args.output,
        base_url=args.base_url,
        merge=args.merge,
        include_optional=args.include_optional,
        fail_on_missing=args.fail_on_missing,
        dry_run=args.dry_run
    )


if __name__ == "__main__":
    main(sys.argv[1:])
