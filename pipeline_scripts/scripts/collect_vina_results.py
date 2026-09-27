#!/usr/bin/env python3

import argparse
import re
from pathlib import Path

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("aa", help="Amino-acid code, e.g. ALA or DALA")
parser.add_argument("--d", action="store_true", help="Use the separate D-AA folders")
args = parser.parse_args()

AA = args.aa.upper()
aa = AA.lower()
affinity = f"vina_{aa}_affinity"
DOCKING_DIR = Path("docking_vina_d" if args.d else "docking_vina") / aa
LOG_DIR = Path(f"logs/vina_d/{aa}" if args.d else f"logs/vina_{aa}")
RESULTS = Path("results/d_amino_acids" if args.d else "results")
RESULTS.mkdir(parents=True, exist_ok=True)

pockets = pd.read_csv(
    "results/pockets_for_docking.tsv",
    sep="\t"
)

results = []
nonempty_outputs = 0

for _, row in pockets.iterrows():

    protein = row["protein"]
    pocket = row["pocket"]

    out = DOCKING_DIR / f"{protein}_{pocket}_{AA}.pdbqt"
    log = LOG_DIR / f"{protein}_{pocket}_{AA}.log"

    if out.is_file() and out.stat().st_size > 0:
        nonempty_outputs += 1

    score = None
    status = "missing"

    if log.exists():

        text = log.read_text(errors="ignore")

        if "MISSING RECEPTOR" in text:
            status = "missing_receptor"

        elif "MISSING LIGAND" in text:
            status = "missing_ligand"

        else:
            m = re.search(
                r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+",
                text,
                re.MULTILINE
            )

            if m:
                score = float(m.group(1))
                status = "success"
            else:
                status = "failed"

    results.append({
        "protein": protein,
        "pocket": pocket,
        affinity: score,
        "vina_status": status
    })


vina = pd.DataFrame(results)

merged = pockets.merge(
    vina,
    on=["protein", "pocket"],
    how="left"
)

merged.to_csv(
    RESULTS / f"vina_{aa}_all_pockets.tsv",
    sep="\t",
    index=False
)

success = merged[
    merged["vina_status"] == "success"
].copy()

# More negative Vina score = better.
if success.empty:
    best = success
else:
    best_idx = success.groupby("protein")[affinity].idxmin()
    best = success.loc[best_idx].sort_values(affinity)

best.to_csv(
    RESULTS / f"vina_{aa}_best_per_protein.tsv",
    sep="\t",
    index=False
)

print("Total expected dockings:", len(merged))
print("Non-empty output files:", nonempty_outputs)
print("Successful:", (merged["vina_status"] == "success").sum())
print("Failed:", (merged["vina_status"] == "failed").sum())
print("Missing receptor:", (merged["vina_status"] == "missing_receptor").sum())
print("Missing ligand:", (merged["vina_status"] == "missing_ligand").sum())
print("Missing/no log:", (merged["vina_status"] == "missing").sum())

print("\nProteins with successful docking:", best["protein"].nunique())

print(f"\nBest 20 {AA} scores:")
print(
    best[
        [
            "gene_id",
            "uniprot_id",
            "protein",
            "pocket",
            "probability",
            "mean_pocket_plddt",
            affinity
        ]
    ].head(20).to_string(index=False)
)
