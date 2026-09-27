#!/usr/bin/env python3

import argparse
import math
from pathlib import Path

import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--d", action="store_true", help="Use the separate D-AA folders")
args = parser.parse_args()

RESULTS = Path("results/d_amino_acids" if args.d else "results")
LIGAND_DIR = Path("ligands/d_amino_acids_pdbqt" if args.d else "ligands/amino_acids_pdbqt")
OUT = RESULTS / "filtered"
(OUT / "vina").mkdir(parents=True, exist_ok=True)
(OUT / "sfct").mkdir(parents=True, exist_ok=True)


def rank_scores(df, score_column):
    ranked = df.copy()
    ranked[score_column] = pd.to_numeric(ranked[score_column], errors="coerce")
    ranked = ranked.dropna(subset=[score_column]).sort_values(score_column)

    n = len(ranked)
    ranked["rank"] = range(1, n + 1)
    ranked["percentile"] = ranked["rank"] / n * 100 if n else None

    if not n:
        ranked["hit_tier"] = ""
        return ranked, 0, 0

    top_1 = max(1, math.ceil(n * 0.01))
    top_5 = max(1, math.ceil(n * 0.05))
    ranked["hit_tier"] = ""
    ranked.loc[ranked["rank"] <= top_5, "hit_tier"] = "top_5pct"
    ranked.loc[ranked["rank"] <= top_1, "hit_tier"] = "top_1pct"

    return ranked, top_1, top_5


thresholds = []
all_sfct_hits = []

ligands = sorted(LIGAND_DIR.glob("*.pdbqt"))

for ligand in ligands:
    AA = ligand.stem.upper()
    if not args.d and AA == "DMET":
        continue  # Ignore the legacy DMET control in the regular-AA folder.
    aa = AA.lower()
    vina_file = RESULTS / f"vina_{aa}_best_per_protein.tsv"
    sfct_file = RESULTS / f"sfct_{aa}_best_per_protein.tsv"

    if not vina_file.exists() or not sfct_file.exists():
        print(f"Skipping {AA}: missing Vina or SFCT result table")
        continue

    vina_score = f"vina_{aa}_affinity"
    vina, vina_top1, vina_top5 = rank_scores(
        pd.read_csv(vina_file, sep="\t"), vina_score
    )
    sfct, sfct_top1, sfct_top5 = rank_scores(
        pd.read_csv(sfct_file, sep="\t"), "vina_sfct_combined"
    )

    for score_type, ranked, score, n1, n5 in [
        ("vina", vina, vina_score, vina_top1, vina_top5),
        ("sfct", sfct, "vina_sfct_combined", sfct_top1, sfct_top5),
    ]:
        thresholds.append({
            "aa": AA,
            "score_type": score_type,
            "score_column": score,
            "successful_proteins": len(ranked),
            "median": ranked[score].median(),
            "top_5pct_cutoff": ranked.iloc[n5 - 1][score] if n5 else None,
            "top_1pct_cutoff": ranked.iloc[n1 - 1][score] if n1 else None,
        })

    vina.insert(0, "aa", AA)
    sfct.insert(0, "aa", AA)
    vina["score_type"] = "vina"
    sfct["score_type"] = "sfct"
    vina["score"] = vina[vina_score]
    sfct["score"] = sfct["vina_sfct_combined"]

    vina_hits = vina[vina["hit_tier"] != ""].copy()
    vina_hits.to_csv(
        OUT / "vina" / f"vina_{aa}_top_hits.tsv", sep="\t", index=False
    )

    vina_support = vina[
        ["protein", "pocket", "score", "rank", "percentile", "hit_tier"]
    ].rename(columns={
        "pocket": "vina_best_pocket",
        "score": "vina_best_score",
        "rank": "vina_rank",
        "percentile": "vina_percentile",
        "hit_tier": "vina_hit_tier",
    })

    sfct = sfct.merge(vina_support, on="protein", how="left")
    sfct["consensus_top_5pct"] = sfct["vina_hit_tier"] != ""
    sfct_hits = sfct[sfct["hit_tier"] != ""].copy()
    sfct_hits.to_csv(
        OUT / "sfct" / f"sfct_{aa}_top_hits.tsv", sep="\t", index=False
    )
    all_sfct_hits.append(sfct_hits)

pd.DataFrame(thresholds).to_csv(
    OUT / "aa_thresholds.tsv", sep="\t", index=False
)

combined = pd.concat(all_sfct_hits, ignore_index=True) if all_sfct_hits else pd.DataFrame()
combined.to_csv(OUT / "all_aa_top_hits.tsv", sep="\t", index=False)

print(f"Thresholds: {OUT / 'aa_thresholds.tsv'}")
print(f"Combined SFCT hits: {OUT / 'all_aa_top_hits.tsv'}")
