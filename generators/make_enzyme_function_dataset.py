import csv, json, hashlib
from pathlib import Path
from typing import List, Dict, Optional
import requests

UNIPROT_ENDPOINT = "https://rest.uniprot.org/uniprotkb/search"

TAXON_ID = 451871  
TARGET_N = 1000
MIN_LEN, MAX_LEN = 100, 800
OUT_DIR = Path("dataset_transferases_eurotiomycetidae_v1")

FIELDS = [
    "accession", "id", "protein_name", "gene_primary", "organism_name",
    "length", "sequence", "ec", "go_id", "protein_existence",
    "xref_pdb", "xref_alphafolddb"
]

def sha256sum(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def fetch_uniprot(query: str, fields: List[str], size: int = 500, max_total: int = 1000) -> List[Dict]:
    params = {"query": query, "format": "json", "fields": ",".join(fields), "size": size}
    out, cursor = [], None
    while True:
        if cursor: params["cursor"] = cursor
        r = requests.get(UNIPROT_ENDPOINT, params=params, timeout=60)
        r.raise_for_status()
        data = r.json()
        for rec in data.get("results", []):
            try:
                out.append(flatten_uniprot(rec))
                if len(out) >= max_total:
                    return out
            except:
                continue
        cursor = data.get("next")
        if not cursor or not data.get("results"):
            break
    return out

def flatten_uniprot(rec: Dict) -> Dict:
    flat: Dict[str, Optional[str]] = {}

    # Accession & ID
    flat["accession"] = rec.get("primaryAccession") or rec.get("accession")
    flat["id"] = rec.get("uniProtkbId") or rec.get("id")

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
            alts = prot_desc.get("alternativeNames") or []
            if alts and isinstance(alts, list):
                alt0 = alts[0] or {}
                full = alt0.get("fullName")
                if isinstance(full, dict):
                    prot_name = full.get("value")
                elif isinstance(full, str):
                    prot_name = full
    else:
        prot_name = rec.get("proteinName") or rec.get("protein_name")
    flat["protein_name"] = prot_name

    # Gene primary
    genes = rec.get("genes")
    if isinstance(genes, list) and genes:
        gene_name = genes[0].get("geneName", {})
        flat["gene_primary"] = gene_name.get("value") if isinstance(gene_name, dict) else gene_name or None
    else:
        flat["gene_primary"] = rec.get("gene_primary") or None

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

    # EC numbers
    ecs = []
    if isinstance(prot_desc, dict):
        rn = prot_desc.get("recommendedName") or {}
        for e in rn.get("ecNumbers") or []:
            v = e.get("value") if isinstance(e, dict) else e
            if v: ecs.append(v)
        for alt in prot_desc.get("alternativeNames") or []:
            for e in alt.get("ecNumbers") or []:
                v = e.get("value") if isinstance(e, dict) else e
                if v: ecs.append(v)

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

    pe = rec.get("proteinExistence")
    if isinstance(pe, dict):
        flat["protein_existence"] = pe.get("evidenceCode") or pe.get("value")
    else:
        flat["protein_existence"] = pe or rec.get("protein_existence")

    return flat

def write_csv(items: List[Dict], path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = FIELDS
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for it in items:
            row = {k: it.get(k) for k in keys}
            w.writerow(row)

def write_fasta(items: List[Dict], path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        for it in items:
            acc = it.get("accession"); seq = it.get("sequence") or ""
            org = it.get("organism_name") or ""
            if not acc or not seq: continue
            header = f">{acc} {org}"
            f.write(header + "\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + "\n")

def main():
    query = (
        f"taxonomy_id:{TAXON_ID} AND ec:2.* AND fragment:false AND length:[{MIN_LEN} TO {MAX_LEN}]"
    )
    print("UniProt query:", query)
    items = fetch_uniprot(query, fields=FIELDS, size=500, max_total=TARGET_N)
    print(f"Fetched: {len(items)} entries")

    # Export
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "transferases_Eurotiomycetidae.csv"
    fasta_path = OUT_DIR / "transferases_Eurotiomycetidae.fasta"
    write_csv(items, csv_path)
    write_fasta(items, fasta_path)

    # Dataset metadata
    meta = {
        "name": "transferases_Eurotiomycetidae",
        "version": "1.0.0",
        "description": "UniProtKB entries filtered by Eurotiomycetidae (taxon 451871), EC 2.*, non-fragments, length 100–800 aa.",
        "query": query,
        "count": len(items),
        "license": "UniProt Terms of Use",
        "source": [{"database": "UniProtKB", "endpoint": UNIPROT_ENDPOINT}],
        "artifacts": [
            {
                "path": csv_path.name,
                "size": csv_path.stat().st_size,
                "sha256": sha256sum(csv_path),
                "kind": "table",
                "format": "csv",
                "primary": True
            },
            {
                "path": fasta_path.name,
                "size": fasta_path.stat().st_size,
                "sha256": sha256sum(fasta_path),
                "kind": "sequences",
                "format": "fasta"
            }
        ]
    }
    (OUT_DIR / "metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print("Wrote:", csv_path, fasta_path, OUT_DIR / "metadata.json")

if __name__ == "__main__":
    main()
