#!/usr/bin/env python3
"""Split an existing AA/PDB dataset into organism directories without API calls."""
import argparse
import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

DIRECTORIES = {"Arabidopsis": "Arabidopsis", "E. coli": "E_coli", "Human": "Human", "Mouse": "Mouse", "Yeast": "Yeast"}
ALL_FILE = "all_AA_PDB_complexes.tsv"
CONTROL_FILE = "nonredundant_AA_protein_controls.tsv"
SUMMARY_FILE = "AA_organism_summary.tsv"
STRICT_FILE = "strict_WT_single_protein_AA_controls.tsv"
ULTRA_FILE = "ultra_strict_WT_monomer_AA_controls.tsv"


def read_tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def write_tsv(path, fields, rows):
    temporary = path.with_suffix(".tsv.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def assigned_groups(row, controls=False):
    if controls:
        # The builder permits unavailable source metadata for a single-source
        # entity; retain the chain organism as a fallback in that case.
        value = row.get("UniProt_organism_group") or row.get("organism_group", "")
    else:
        # Unresolved contact rows have no assigned chain. Preserve them under
        # each organism present in the entry, retaining their original status.
        value = row.get("organism_group") or row.get("entry_organism_groups", "")
    return set(value.split(";")) - {""}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=Path("."), help="Directory containing the combined TSV files")
    args = parser.parse_args(argv)
    root = args.out.resolve()
    tables = {name: read_tsv(root / name) for name in (ALL_FILE, CONTROL_FILE, SUMMARY_FILE)}
    for name in (STRICT_FILE, ULTRA_FILE):
        if (root / name).exists():
            tables[name] = read_tsv(root / name)
    source_hashes = {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in tables}
    for name in [n for n in tables if n != SUMMARY_FILE]:
        for row in tables[name][1]:
            groups = assigned_groups(row, controls=name != ALL_FILE)
            if not groups or not groups.issubset(DIRECTORIES):
                raise ValueError(f"Unrecognized or missing organism assignment in {name}: {row['PDB_ID']}: {groups}")
            for field in ("raw_response_file", "strict_metadata_response_file"):
                if row.get(field) and not (root / row[field]).is_file():
                    raise FileNotFoundError(row[field])
    for organism, directory in DIRECTORIES.items():
        destination = root / directory
        destination.mkdir(exist_ok=True)
        counts = {}
        for name, (fields, rows) in tables.items():
            if name == SUMMARY_FILE:
                selected = [dict(row) for row in rows if row["organism_group"] == organism]
            else:
                selected = [dict(row) for row in rows if organism in assigned_groups(row, controls=name != ALL_FILE)]
            for row in selected:
                for field in ("raw_response_file", "strict_metadata_response_file"):
                    if row.get(field):
                        row[field] = os.path.relpath(root / row[field], destination).replace(os.sep, "/")
            write_tsv(destination / name, fields, selected)
            counts[name] = {"rows": len(selected)}
            if name != SUMMARY_FILE:
                counts[name]["structures"] = len({row["PDB_ID"] for row in selected})
                counts[name]["rows_without_confirmed_target_contact"] = sum(row["status"] != "target_protein_contact" for row in selected)
        manifest = {
            "organism": organism,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "source_dataset_manifest": "../dataset_manifest.json",
            "source_file_sha256": source_hashes,
            "counts": counts,
            "output_sha256": {name: hashlib.sha256((destination / name).read_bytes()).hexdigest() for name in tables},
            "assignment": "Complexes use contacting-chain organism; unresolved rows use entry organisms. Controls use UniProt organism, with chain organism fallback when unavailable.",
            "shared_entries": "Rows belonging to multiple organisms appear in each applicable directory.",
            "raw_data": "Shared raw archive remains in the parent directory; raw_response_file paths are relative to this directory.",
        }
        (destination / "organism_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        print(f"{directory}: {counts[ALL_FILE]['rows']} complex rows; {counts[CONTROL_FILE]['rows']} controls")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
