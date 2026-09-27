#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("aa", help="Amino-acid code, e.g. ALA or DALA")
parser.add_argument("--d", action="store_true", help="Use the separate D-AA folders")
args = parser.parse_args()

AA = args.aa.upper()
aa = AA.lower()
POCKETS = "results/pockets_for_docking.tsv"
SFCT_DIR = Path("docking_sfct_d" if args.d else "docking_sfct") / aa
DOCKING_DIR = Path("docking_vina_d" if args.d else "docking_vina") / aa
LOG_DIR = Path(f"logs/sfct_d/{aa}" if args.d else f"logs/sfct_{aa}")
RESULTS = Path("results/d_amino_acids" if args.d else "results")
RESULTS.mkdir(parents=True, exist_ok=True)

pockets = pd.read_csv(POCKETS, sep="\t")

records = []
vina_outputs = [
    f for f in sorted(DOCKING_DIR.glob(f"*_{AA}.pdbqt"))
    if f.stat().st_size > 0
]

for ligand in vina_outputs:

    # Example:
    # AF-Q8LPL7-F1-model_v6_pocket1_ALA.pdbqt

    stem = ligand.stem
    prefix = stem[:-(len(AA) + 1)]

    try:
        protein, pocket = prefix.rsplit("_", 1)
    except ValueError:
        print("Could not parse:", ligand)
        continue

    f = SFCT_DIR / f"{stem}.dat"
    log = LOG_DIR / f"{stem}.log"

    if not f.is_file() or f.stat().st_size == 0:
        status = "missing"

        if log.is_file():
            text = log.read_text(errors="ignore")
            status = (
                "missing_receptor"
                if "MISSING RECEPTOR" in text
                else "failed"
            )

        records.append({
            "protein": protein,
            "pocket": pocket,
            "sfct_status": status
        })
        continue

    try:
        df = pd.read_csv(
            f,
            sep=r"\s+",
            comment="#",
            header=None,
            names=[
                "pose_name",
                "pose_index",
                "origin_score",
                "combined_score",
                "sfct_score"
            ]
        )

        if df.empty:
            records.append({
                "protein": protein,
                "pocket": pocket,
                "sfct_status": "failed"
            })
            continue

        # Best SFCT-corrected pose:
        # lower combined score = better
        best = df.loc[df["combined_score"].idxmin()]

        records.append({
            "protein": protein,
            "pocket": pocket,
            "sfct_best_pose": int(best["pose_index"]),
            "sfct_vina_score": float(best["origin_score"]),
            "sfct_score": float(best["sfct_score"]),
            "vina_sfct_combined": float(best["combined_score"]),
            "sfct_n_poses": len(df),
            "sfct_status": "success"
        })

    except Exception as e:
        print("ERROR:", f, e)
        records.append({
            "protein": protein,
            "pocket": pocket,
            "sfct_status": "failed"
        })


sfct = pd.DataFrame(records, columns=[
    "protein", "pocket", "sfct_best_pose", "sfct_vina_score",
    "sfct_score", "vina_sfct_combined", "sfct_n_poses", "sfct_status"
])

print(f"\n=== SFCT {AA} summary ===")
print("Expected jobs:", len(vina_outputs))
print("Successful:", (sfct["sfct_status"] == "success").sum())
print("Missing receptor:", (sfct["sfct_status"] == "missing_receptor").sum())
print("Failed:", (sfct["sfct_status"] == "failed").sum())
print("Missing/no log:", (sfct["sfct_status"] == "missing").sum())


# ---------------------------------------------------------
# Merge with P2Rank / pLDDT / original docking information
# ---------------------------------------------------------

merged = pockets.merge(
    sfct,
    on=["protein", "pocket"],
    how="left"
)

merged["sfct_status"] = merged["sfct_status"].fillna("missing")

merged.to_csv(
    RESULTS / f"sfct_{aa}_all_pockets.tsv",
    sep="\t",
    index=False
)


# ---------------------------------------------------------
# Best corrected pocket per protein
# ---------------------------------------------------------

ok = merged[
    merged["sfct_status"] == "success"
].copy()

if ok.empty:
    best = ok
else:
    best_idx = ok.groupby("protein")["vina_sfct_combined"].idxmin()
    best = ok.loc[best_idx].sort_values("vina_sfct_combined")

best.to_csv(
    RESULTS / f"sfct_{aa}_best_per_protein.tsv",
    sep="\t",
    index=False
)


print("Successful SFCT pockets:", len(ok))
print("Proteins with SFCT result:", best["protein"].nunique())

print("\nCombined score distribution:")
print(
    best["vina_sfct_combined"].describe(
        percentiles=[.01,.05,.10,.25,.50,.75,.90,.95,.99]
    )
)

print("\nTop 20 after SFCT:")
print(
    best[
        [
            "gene_id",
            "uniprot_id",
            "pocket",
            "probability",
            "mean_pocket_plddt",
            "sfct_vina_score",
            "sfct_score",
            "vina_sfct_combined"
        ]
    ].head(20).to_string(index=False)
)
