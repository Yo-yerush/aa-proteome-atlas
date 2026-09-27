#!/usr/bin/env python3

import argparse
import csv
import re
import statistics
from pathlib import Path

PRED_DIR = Path("pockets_p2rank")
PDB_DIR  = Path("structures_clean")
MANIFEST = Path("metadata/structure_manifest.tsv")
OUTPUT   = Path("results/master_pockets.tsv")


# ---------------------------------------------------------
# Arguments
# ---------------------------------------------------------

parser = argparse.ArgumentParser(
    description="Build a master pocket table with a standardized gene_id column."
)
parser.add_argument(
    "--gene-id-column",
    required=True,
    help="Gene identifier column in metadata/structure_manifest.tsv",
)
args = parser.parse_args()


# ---------------------------------------------------------
# UniProt -> gene ID mapping
# ---------------------------------------------------------

uniprot_to_gene_id = {}

if MANIFEST.exists():
    with MANIFEST.open() as f:
        reader = csv.DictReader(f, delimiter="\t")

        required_columns = {"uniprot_id", args.gene_id_column}
        missing_columns = required_columns.difference(reader.fieldnames or [])
        if missing_columns:
            missing = ", ".join(sorted(missing_columns))
            raise ValueError(f"Missing manifest column(s): {missing}")

        for row in reader:
            uid = row.get("uniprot_id", "").strip()
            gene_id = row.get(args.gene_id_column, "").strip()
            if uid:
                uniprot_to_gene_id[uid] = gene_id


# ---------------------------------------------------------
# Read AlphaFold pLDDT from PDB B-factor field
# One value per residue, using CA atom
# ---------------------------------------------------------

def read_plddt(pdb_file):

    plddt = {}

    with open(pdb_file) as f:

        for line in f:

            if not line.startswith("ATOM"):
                continue

            atom = line[12:16].strip()

            if atom != "CA":
                continue

            chain = line[21].strip()
            if not chain:
                chain = "A"

            try:
                residue_number = int(line[22:26])
                bfactor = float(line[60:66])
            except ValueError:
                continue

            key = f"{chain}_{residue_number}"

            plddt[key] = bfactor

    return plddt


# ---------------------------------------------------------
# Output
# ---------------------------------------------------------

OUTPUT.parent.mkdir(parents=True, exist_ok=True)

columns = [
    "uniprot_id",
    "gene_id",
    "protein",
    "pocket",
    "rank",
    "score",
    "probability",
    "center_x",
    "center_y",
    "center_z",
    "sas_points",
    "surf_atoms",
    "n_pocket_residues",
    "n_plddt_matched",
    "mean_pocket_plddt",
    "min_pocket_plddt",
    "fraction_plddt_ge70",
    "fraction_plddt_ge90",
    "residue_ids"
]


prediction_files = sorted(PRED_DIR.rglob("*_predictions.csv"))

print(f"Found {len(prediction_files)} prediction files")


with OUTPUT.open("w", newline="") as out:

    writer = csv.DictWriter(
        out,
        fieldnames=columns,
        delimiter="\t"
    )

    writer.writeheader()

    for i, pred_file in enumerate(prediction_files, 1):

        # Example:
        # AF-ACCESSION-F2-model_v6.pdb_predictions.csv
        protein_file = pred_file.name.replace("_predictions.csv", "")
        protein_name = Path(protein_file).stem

        m = re.search(r"^AF-(.+?)-F\d+(?:-|$)", protein_name)

        if m:
            uniprot = m.group(1)
        else:
            uniprot = ""

        gene_id = uniprot_to_gene_id.get(uniprot, "")

        pdb_file = PDB_DIR / protein_file

        if not pdb_file.exists():
            pdb_file = PDB_DIR / f"{protein_name}.pdb"

        if not pdb_file.exists():
            print(f"WARNING: missing PDB: {protein_file}")
            continue

        plddt = read_plddt(pdb_file)

        with pred_file.open() as f:

            reader = csv.DictReader(f)

            for raw in reader:

                # P2Rank headers contain spaces
                row = {k.strip(): v.strip() for k, v in raw.items()}

                residues = row.get("residue_ids", "").split()

                pocket_plddt = [
                    plddt[r]
                    for r in residues
                    if r in plddt
                ]

                if pocket_plddt:

                    mean_plddt = statistics.mean(pocket_plddt)
                    min_plddt = min(pocket_plddt)

                    frac70 = sum(
                        x >= 70 for x in pocket_plddt
                    ) / len(pocket_plddt)

                    frac90 = sum(
                        x >= 90 for x in pocket_plddt
                    ) / len(pocket_plddt)

                else:

                    mean_plddt = ""
                    min_plddt = ""
                    frac70 = ""
                    frac90 = ""

                writer.writerow({

                    "uniprot_id": uniprot,
                    "gene_id": gene_id,
                    "protein": protein_name,

                    "pocket": row.get("name", ""),
                    "rank": row.get("rank", ""),
                    "score": row.get("score", ""),
                    "probability": row.get("probability", ""),

                    "center_x": row.get("center_x", ""),
                    "center_y": row.get("center_y", ""),
                    "center_z": row.get("center_z", ""),

                    "sas_points": row.get("sas_points", ""),
                    "surf_atoms": row.get("surf_atoms", ""),

                    "n_pocket_residues": len(residues),
                    "n_plddt_matched": len(pocket_plddt),

                    "mean_pocket_plddt":
                        f"{mean_plddt:.2f}" if pocket_plddt else "",

                    "min_pocket_plddt":
                        f"{min_plddt:.2f}" if pocket_plddt else "",

                    "fraction_plddt_ge70":
                        f"{frac70:.3f}" if pocket_plddt else "",

                    "fraction_plddt_ge90":
                        f"{frac90:.3f}" if pocket_plddt else "",

                    "residue_ids": " ".join(residues)
                })


        if i % 1000 == 0:
            print(f"Processed {i}/{len(prediction_files)} proteins")


print()
print(f"Finished.")
print(f"Master table: {OUTPUT}")
