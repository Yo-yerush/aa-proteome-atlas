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
affinity = f"vina_{aa}_affinity"
RESULTS = Path("results/d_amino_acids" if args.d else "results")

vina = pd.read_csv(RESULTS / f"vina_{aa}_best_per_protein.tsv", sep="\t")
sfct = pd.read_csv(RESULTS / f"sfct_{aa}_best_per_protein.tsv", sep="\t")

vina = vina.sort_values(affinity)
vina["vina_rank"] = range(1, len(vina) + 1)
vina["vina_percentile"] = vina["vina_rank"] / len(vina) * 100

sfct = sfct.sort_values("vina_sfct_combined")
sfct["sfct_rank"] = range(1, len(sfct) + 1)
sfct["sfct_percentile"] = sfct["sfct_rank"] / len(sfct) * 100

comp = vina[
    ["protein", "gene_id", "uniprot_id", affinity,
     "vina_rank", "vina_percentile"]
].merge(
    sfct[
        ["protein", "pocket", "probability", "mean_pocket_plddt",
         "sfct_vina_score", "sfct_score", "vina_sfct_combined",
         "sfct_rank", "sfct_percentile"]
    ],
    on="protein",
    how="inner"
)

comp["rank_change"] = comp["vina_rank"] - comp["sfct_rank"]
comp.to_csv(RESULTS / f"vina_vs_sfct_{aa}.tsv", sep="\t", index=False)

print("Proteins compared:", len(comp))
print("\nCorrelation:")
print(comp[[affinity, "vina_sfct_combined"]].corr())
print("\nLargest improvements after SFCT:")
print(comp.sort_values("rank_change", ascending=False).head(20).to_string(index=False))
