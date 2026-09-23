<img src="assets/favicon.svg" width="64" height="64" align="right" alt="AA Proteome Atlas icon" />

# AA Proteome Interaction Atlas

Explore predicted amino-acid binding across the proteome.

**[Open the web app →](https://yo-yerush.github.io/aa-proteome-atlas/)**

AA Proteome Interaction Atlas is an interactive web app for finding, comparing and investigating candidate proteins that may bind free amino acids. It brings together proteome-wide docking scores, predicted binding pockets, 3D ligand poses, functional annotations and experimental positive controls in one place.

The current release contains **Arabidopsis thaliana** and **Escherichia coli** results, each with **20 canonical amino acids** and **19 D-amino-acid datasets** for stereochemical comparisons. Glycine has no distinct L/D pair.

Choose an organism in the header. E. coli is the default for links without an organism parameter; `?organism=arabidopsis` opens Arabidopsis directly. Switching reloads the atlas with only that organism's data, preserves shared filters and the current tab, and clears protein searches, selected proteins/pockets and all analysis/viewer caches. Rankings, percentiles, normalization, GO backgrounds and experimental-control QC are computed separately within each organism. Optional files load on demand; missing datasets are shown as **Unavailable** without substituting another organism's data.

E. coli uses the supplied gene symbols and `b`/`JW` locus identifiers, with product descriptions keyed by `b` locus ID. UniProt links remain available for both organisms. Download filenames begin with `arabidopsis_` or `ecoli_`; protein exports use `tair_id` for Arabidopsis and `gene_id` for E. coli.

Use the atlas to move from a ranked protein list to a specific pocket, compare amino-acid preferences, explore enriched biological functions, and assess how the scoring performs against known controls.

> The atlas supports hypothesis generation. Docking scores indicate predicted compatibility—not experimentally confirmed binding, binding probability or measured affinity.

[Explore the app](#explore-the-app) · [Analysis settings](#analysis-settings) · [Scoring and interpretation](#scoring-and-interpretation) · [Run the app](#run-the-app) · [Technical guide](docs/technical-guide.md)

## Explore the app

| Tab | What you can do |
| --- | --- |
| **Explorer** | Search proteins or genes; rank candidates for a target AA; filter by pocket quality, proteome tier, selectivity and L/D preference; switch between best and all pockets. |
| **Protein profile** | Inspect the full 20-AA profile, compare scores or AA-normalized values, and view the protein, pocket residues, P2Rank points and docked ligand in Mol*. |
| **AA matrix** | Compare proteome-wide AA score patterns with a 20 × 20 Spearman correlation matrix. |
| **Compare AAs** | Compare any two available L/D AAs, inspect paired scores and identify proteins predicted to favor one over the other. |
| **Top-hit overlap** | Find shared or unique top-ranked proteins across canonical AAs using intersections, Jaccard similarity and Venn/UpSet views. |
| **GO analysis** | Test for overrepresented Biological Process, Molecular Function and Cellular Component terms in a selected protein set. |
| **Control QC** | Measure recovery of known AA–protein controls, inspect experimental-site overlap and compare L/D preferences. |
| **Statistics** | Review score distributions, pocket-quality metrics, retained-pocket counts and score correlations. |
| **Methods** | Read the calculation rules, parameter definitions, data-handling conventions and interpretation limits. |

Protein rows link to UniProt, provide gene descriptions and open the corresponding profile. Download filtered tables in the formats offered by each view, including TSV/CSV, and export the protein-profile and GO plots as SVG or PNG.

### A typical workflow

1. **Find candidates.** Choose a target AA in Explorer, apply the quality filters, and narrow the list to a Top 1%, 5% or 10% tier.
2. **Check specificity.** Review the within-protein AA rank, competitive AAs, Δ to median and optional L>D filter.
3. **Inspect the pocket.** Open a protein profile, inspect individual pockets, compare the 20-AA scores and examine the exported ligand poses in 3D.
4. **Put the results in context.** Explore functional enrichment in GO analysis, compare hit lists across AAs, and review Control QC before prioritizing experimental follow-up.

## Analysis settings

These are the default app settings; they can be changed in the relevant tab.

| Setting | Default | Meaning |
| --- | --- | --- |
| Target amino acid | Alanine (ALA) in Explorer and GO | The AA used to select and rank candidate proteins. |
| Ranking score | Combined 50% | Lower scores rank better. Vina, SFCT and Combined 80% are also available. |
| P2Rank probability | ≥ 0.7 | Retain pockets meeting the predicted-pocket probability threshold. |
| Pocket mean pLDDT | ≥ 90 | Retain pockets meeting the local structure-confidence threshold. |
| Proteome tier | All ranks in Explorer; Top 5% in GO | Restrict the list using the target AA’s QC-qualified protein ranking. |
| Pocket display | Best pocket in Explorer | One best qualifying pocket per protein; All pockets shows individual qualifying pocket rows. |
| Competitive AAs maximum | 19 in Explorer and GO | Initially unrestricted; lower values narrow the selected protein set. |
| L>D | Off | Require the L-AA to score better than its matching D-AA. |
| Protein-profile plot | Combined 50%, best-to-worst score order; normalization off | Change Value type, ordering or AA-normalized Z-score independently of the Explorer ranking score. |

Quality filters apply **before best-pocket selection and ranking**. Protein profile and AA matrix share Explorer’s quality thresholds; Compare AAs, Top-hit overlap, GO, Control QC and Statistics have their own analysis settings.

The input dataset’s upstream pocket-retention settings, as documented in Methods, are less restrictive: mean pocket pLDDT ≥70, at least 80% of pocket residues with pLDDT ≥70, and P2Rank probability ≥0.20 **or** rank ≤3 with probability ≥0.05, with at most 10 retained pockets per protein. These are distinct from the app defaults; relaxing an app filter cannot restore pockets absent from the exported data.

## Scoring and interpretation

### Available scores

All four ranking scores use **lower is better**.

| Score | Definition |
| --- | --- |
| **Vina** | The independently selected best AutoDock Vina score for the pocket, reported in kcal/mol. |
| **SFCT** | The OnionNet-SFCT score of the saved SFCT/Combined-selected pose. |
| **Combined 80%** | `0.2 × origin_score + 0.8 × sfct_score`, as supplied by the pipeline. |
| **Combined 50%** | `0.5 × origin_score + 0.5 × sfct_score`, calculated from the saved pose’s values. |

Here, `origin_score` is the saved pose’s Vina score (`sfct_vina_score`), **not necessarily the independent Vina minimum**. SFCT and both Combined values describe the pose selected upstream using Combined 80%. Combined 50% reweights that pose and can change the best pocket; it does not search again for a different best pose.

### Ranking, selectivity and missing results

- **Proteome percentile:** rank within the same AA and score after quality filtering: `100 × (rank − 1) / (N − 1)`, with a single protein assigned 0. Lower is better. Ties are ordered by UniProt ID.
- **Within-protein AA rank:** rank of the target score among the available successful canonical AA scores for that protein. Inspecting a pocket restricts comparisons to the same pocket geometry.
- **Competitive AAs:** other canonical AAs with `other score ≤ target score + tolerance`. Tolerances are 0.10 for either Combined score, 0.25 for Vina and 0.15 for SFCT. This is a score-based heuristic, not a count of experimentally established binders.
- **Δ to median:** `median(other AA scores) − target score`. Positive values favor the target.
- **Selectivity Z:** `[mean(other AA scores) − target score] / SD(other AA scores)`. Positive values favor the target.
- **L>D:** requires `L score < D score`, equivalently a positive `D − L` difference. In Explorer and GO, L and D use their independently best QC-passing pockets. Ties and missing pairs do not pass; D-AA results are comparisons, not proven nonbinding controls.

The **Top-hit overlap** table uses a different competitive-AA count: `x/20` is membership in the top-tier lists for all 20 canonical AAs under that tab’s settings.

Only successful, finite scores contribute to rankings and statistics. Missing or failed results for one protein × AA do not remove successful results for another AA. Missing values are never replaced with zero. D-AA datasets do not enter the canonical 20-AA ranks, selectivity statistics or correlation matrix.

### AA-normalized profiles

Raw docking-score baselines differ between AAs. The protein-profile plot offers an **AA-normalized Z-score** checkbox for each raw score type:

`Z = (score − AA-specific mean) / AA-specific population SD`

The reference is one best QC-passing pocket per protein for that AA and score. **Negative is favorable here**, unlike the within-protein Selectivity Z above. The plot marks `Z = 0` and each AA’s observed Top 5% boundary; the latter is not a fixed normal-distribution significance cutoff. Normalization changes the plot, not the raw score table or ranking calculations.

### GO enrichment

GO analysis performs a **one-sided Fisher exact test** (hypergeometric upper tail), followed by **Benjamini–Hochberg FDR correction** across all background GO terms from BP, MF and CC, including terms with zero selected hits.

- **Defaults:** ALA, Top 5%, Combined 50%, P2Rank ≥0.7, pLDDT ≥90, Competitive AAs ≤19, L>D off, Biological Process, and FDR ≤0.05. The significant-only view is enabled by default.
- **Background:** Parameter-filtered by default—all proteins with a successful target-AA score passing the same quality filters, without the tier, competitor or L>D restrictions. Total (target AA) removes the quality filters from the background, not from the selected hits.
- **Counting:** selected `k/n` and background `K/N` use proteins annotated in the corresponding GO aspect. Fold enrichment is `(k/n) / (K/N)`; the background includes the selected proteins.
- **Annotations:** only the supplied UniProt GO mappings are analyzed. There is no additional ancestor propagation, evidence-code filtering or live enrichment API call.

FDR refers to the **GO enrichment tests**, not to the probability that a docking candidate is a false-positive binder. Changing the GO aspect or display threshold does not change the testing family; changing the selected/background protein populations recalculates enrichment.

### Experimental-control QC

Control QC compares Vina, SFCT and both Combined scores against known AA–protein controls. Curves show cumulative recovery versus top proteome percentile, with Top 1%, 5% and 10% summaries and an illustrative random-reference diagonal.

The default is **Best-scoring pocket**, with P2Rank ≥0.7 and pLDDT ≥90. **Control-matched pocket** instead uses the predicted pocket with the largest overlap with mapped experimental-site residues; the default requirement is at least 50% coverage. The default curve denominator includes controls evaluable for all four scores; **All known controls** also includes missing or QC-excluded cases.

Site-overlap recovery is reported separately from score recovery. Control-matched analysis is conditioned on the known site, so passing the site threshold is not independent validation. Only retained pockets are available: this is not complete P2Rank sensitivity or ligand-pose accuracy.

## Pocket visualization

The Mol* viewer combines a version-matched AlphaFold structure with the selected pocket and exported ligand coordinates.

- Switch the protein between **cartoon and surface**, with **pLDDT, residue charge, residue polarity or uniform** coloring.
- Toggle the green pocket-residue highlight and optional sticks independently; sticks are off by default.
- Show the selected pocket’s **P2Rank points**, colored gold or by precomputed electrostatic potential on a fixed **−5 → 0 → +5 kT/e** red–white–blue scale.
- Use **Ligand pose: Auto** to follow the profile’s Value type: Vina uses MODEL 1; SFCT/Combined use the saved selection. Manual Vina-best and saved-pose choices change the displayed ligand only—not the scores or plots.

Explorer’s **qφ (kT)** is the precomputed sum of ligand charge × receptor potential for **Vina MODEL 1** in the displayed pocket. Negative values are favorable and positive values unfavorable in the fixed receptor field. It is not binding free energy and does not follow manual pose switching. Failed or unavailable values remain missing, never zero.

## Run the app

The atlas is a static web app: no build step, package installation or backend is required. From the app folder:

```sh
python -m http.server 8765 --bind 127.0.0.1
```

Open **http://127.0.0.1:8765/** in an up-to-date browser. The 3D viewer requires an internet connection for Mol* and AlphaFold structures.

For GitHub Pages, publish `index.html`, `js/`, `css/`, `assets/`, `At_results/`, `Ec_results/`, `annotations/arabidopsis/` and `annotations/ecoli/` together, preserving filename capitalization. See the [deployment instructions](docs/technical-guide.md#github-pages) for the required files and configuration. Confirm data-redistribution permissions and attribution before publishing.

## Troubleshooting and development

- **The app will not load when opening the HTML file directly:** use the local server above or HTTPS hosting; `file://` is not supported.
- **Loading is slow:** all available score datasets for the selected organism are expanded in browser memory at startup. The progress bar tracks completed steps, not remaining time; the first 3D pocket view also loads additional coordinates. The organism selector remains usable during loading.
- **A structure or overlay is unavailable:** check the panel’s message, network access and matching data exports. Use the offered retry control; missing optional overlays do not invalidate the docking scores.
- **A checksum fails after updating files:** publish the complete matching bundle and manifest. Do not mix exports or recompress electrostatics files without updating their checksums.

Keep `tests/` in the repository for regression checks; the browser does not load it and deployment does not require it. Node.js is only needed to [run the tests](docs/technical-guide.md#regression-checks).

The [technical guide](docs/technical-guide.md) preserves the full setup, data-format, export, caching, validation and implementation notes. The app’s **Methods** tab provides detailed scientific definitions and interpretation guidance.
