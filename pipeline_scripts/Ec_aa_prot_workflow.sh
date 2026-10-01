# Whole-proteome amino acids reverse docking workflow in Escherichia coli

mkdir -p /home/yoyerush/yo/whole_at_proteins_docking/
cd /home/yoyerush/yo/whole_at_proteins_docking/

################################################################

### 1. Keep UniProt IDs + E. coli locus tags associated with every structure
# download AlphaFold structures
# https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/

mkdir -p ecoli_docking/structures_raw
cd ecoli_docking/structures_raw

wget https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000000625_83333_ECOLI_v6.tar

tar -xf UP000000625_83333_ECOLI_v6.tar
gunzip *.gz

# extract the UniProt accession from every AlphaFold filename
cd ../
mkdir -p metadata structures_clean receptors_pdbqt logs

find structures_raw -type f \( -name "*.cif" -o -name "*.pdb" \) \
    | sed 's#.*/##' \
    | sed -E 's/^AF-([A-Z0-9]+)-F[0-9]+-model.*$/\1/' \
    | sort -u \
    > metadata/uniprot_ids.txt

# Download the E.coli UniProt annotation table:
wget -O metadata/ecoli_uniprot.tsv \
"https://rest.uniprot.org/uniprotkb/stream?query=%28organism_id%3A83333%29&format=tsv&fields=accession,id,gene_names,gene_oln,length"

# Download the E.coli UniProt GO annotation table:
wget -O metadata/ecoli_uniprot_go.tsv \
"https://rest.uniprot.org/uniprotkb/stream?query=%28organism_id%3A83333%29&format=tsv&fields=accession,go_p,go_f,go_c"
gzip metadata/ecoli_uniprot_go.tsv

# create a structure manifest:
python3 - <<'PY'
import os
import re
import csv

raw = "structures_raw"
annotation = "metadata/ecoli_uniprot.tsv"
output = "metadata/structure_manifest.tsv"

mapping = {}

with open(annotation) as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        acc = row["Entry"]

        locus_tag = ""
        for key in row:
            if key == "Gene Names (ordered locus)":
                locus_tag = row[key].strip()
                break

        genes = row.get("Gene Names", "")

        mapping[acc] = (locus_tag, genes)

with open(output, "w") as out:
    out.write("uniprot_id\tlocus_tag\tgene_names\tstructure_file\n")

    for filename in sorted(os.listdir(raw)):

        if not filename.endswith((".cif", ".pdb")):
            continue

        m = re.match(r"AF-([A-Z0-9]+)-F\d+-model", filename)

        if not m:
            continue

        accession = m.group(1)

        locus_tag, genes = mapping.get(accession, ("", ""))

        out.write(
            f"{accession}\t{locus_tag}\t{genes}\t{filename}\n"
        )
PY

# Check:
column -t -s $'\t' metadata/structure_manifest.tsv | head

################################################################

### 2. Preprocess all AlphaFold structures

### using Meeko + ProDy for receptor preparation

# conda create -n docking python=3.11 -y
conda activate docking
# pip install meeko prody gemmi
# conda install -c conda-forge rdkit parallel -y

# convert all AlphaFold CIF files to standardized PDB files:
python3 - <<'PY'
import os
import gemmi

input_dir = "structures_raw"
output_dir = "structures_clean"

os.makedirs(output_dir, exist_ok=True)

for fn in os.listdir(input_dir):

    if not fn.endswith(".cif"):
        continue

    infile = os.path.join(input_dir, fn)

    structure = gemmi.read_structure(infile)

    outfile = os.path.join(
        output_dir,
        fn.replace(".cif", ".pdb")
    )

    structure.write_pdb(outfile)

print("Conversion complete.")
PY

# Check the number of structures
find structures_clean -name "*.pdb" | wc -l

# prepare every receptor for Vina - using 80 cores in parallel
mkdir -p receptors_pdbqt logs logs/receptor_preparation

find structures_clean -name "*.pdb" | \
parallel -j 80 '
    base=$(basename {} .pdb)
    echo "Preparing $base"
    mk_prepare_receptor.py \
        --read_pdb "{}" \
        -o "receptors_pdbqt/$base" \
        -p \
        > "logs/receptor_preparation/${base}.log" 2>&1
'

# Check how many receptors succeeded:
find receptors_pdbqt -name "*.pdbqt" | wc -l
# total:
find logs/receptor_preparation -name "*.log" | wc -l

# remove the structures_raw extracted files (still keep the original tar file)
rm -f structures_raw/AF-*.pdb
rm -f structures_raw/AF-*.cif

################################################################

### 3. Predict potential binding pockets
# conda install -c conda-forge openjdk=17 fpocket -y

mkdir -p tools
cd tools

# download P2Rank (v2.5.1)
wget https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz
tar -xzf p2rank_2.5.1.tar.gz

cd /home/yoyerush/yo/whole_at_proteins_docking/ecoli_docking

## Run P2Rank on ALL proteins
# Create the dataset:
find structures_clean \
    -name "*.pdb" \
    -type f \
    | sort \
    > ecoli_structures.ds

mkdir -p pockets_p2rank logs

tools/p2rank_2.5.1/prank predict \
    -c alphafold \
    -threads 80 \
    -visualizations 1 \
    -o pockets_p2rank \
    ecoli_structures.ds \
    > logs/p2rank.log 2>&1

# check for 'center_x/y/z' coordinates
find pockets_p2rank -name "*_predictions.csv" | wc -l

# move the raw '*_points.pdb.gz' files
mkdir -p visualizations_p2rank
find pockets_p2rank/visualizations/data \
    -maxdepth 1 -type f -name '*_points.pdb.gz' \
    -exec mv -t visualizations_p2rank -- {} +

rm -r pockets_p2rank/visualizations

################################################################

### 4. Rank and filter pockets

mkdir -p results

# build the master-pocket table
python scripts/build_pocket_master.py --gene-id-column locus_tag


# Before filtering, inspect the distribution
python - <<'PY'
import pandas as pd

x = pd.read_csv(
    "results/master_pockets.tsv",
    sep="\t"
)

print("\nTotal pockets:")
print(len(x))

print("\nProteins:")
print(x["protein"].nunique())

print("\nPockets per protein:")
print(x.groupby("protein").size().describe())

print("\nP2Rank probability:")
print(x["probability"].describe(
    percentiles=[.10,.25,.50,.75,.90,.95,.99]
))

print("\nPocket mean pLDDT:")
print(x["mean_pocket_plddt"].describe(
    percentiles=[.10,.25,.50,.75,.90,.95,.99]
))
PY

# results # Total pockets:
# results # 18298
# results # 
# results # Proteins:
# results # 3589


## Run the pocket QC/filtering
# parameters:
# mean_pLDDT >= 70
# AND
# fraction_pLDDT_ge70 >= 0.80
# AND
# (
#     probability >= 0.20
#     OR
#     (rank <= 3 AND probability >= 0.05)
# )

python - <<'PY'
import pandas as pd

x = pd.read_csv("results/master_pockets.tsv", sep="\t")

# structural quality filter
good_structure = (
    (x["mean_pocket_plddt"] >= 70) &
    (x["fraction_plddt_ge70"] >= 0.80)
)

# pocket retention rule
keep = good_structure & (
    (x["probability"] >= 0.20) |
    (
        (x["rank"] <= 3) &
        (x["probability"] >= 0.05)
    )
)

y = x[keep].copy()

# Keep top 10 pockets per protein
y = (
    y.sort_values(["protein", "rank"])
     .groupby("protein", group_keys=False)
     .head(10)
     .copy()
)

print("Total predicted pockets:", len(x))
print("Retained pockets:", len(y))
print("Retained proteins:", y["protein"].nunique())

print("\nPockets per retained protein:")
print(y.groupby("protein").size().describe())

print("\nRetained probability distribution:")
print(y["probability"].describe(
    percentiles=[.1,.25,.5,.75,.9,.95,.99]
))

print("\nRetained mean pocket pLDDT:")
print(y["mean_pocket_plddt"].describe(
    percentiles=[.1,.25,.5,.75,.9,.95,.99]
))

y.to_csv(
    "results/pockets_for_docking.tsv",
    sep="\t",
    index=False
)

print("\nSaved to results/pockets_for_docking.tsv")
PY

################################################################

### 5. Dock amino acids to all pockets

mkdir -p ligands

# Generate AAs in its zwitterionic form, which is the appropriate primary form for free AAs near physiological pH:
python scripts/create_amino_acid_sdf.py

# Convert L-AAs to PDBQT, skipping any legacy DMET files.
mkdir -p ligands/amino_acids_pdbqt
for sdf in ligands/amino_acids_sdf/*.sdf
do
    [ -f "$sdf" ] || continue
    code=$(basename "$sdf" .sdf)
    [ "${code^^}" = "DMET" ] && continue
    echo "Preparing $code"
    mk_prepare_ligand.py \
        -i "$sdf" \
        -o "ligands/amino_acids_pdbqt/${code}.pdbqt"
done

# convert D-AAs to PDBQT:
mkdir -p ligands/d_amino_acids_pdbqt
for sdf in ligands/d_amino_acids_sdf/*.sdf
do
    [ -f "$sdf" ] || continue
    code=$(basename "$sdf" .sdf)
    echo "Preparing $code"
    mk_prepare_ligand.py \
        -i "$sdf" \
        -o "ligands/d_amino_acids_pdbqt/${code}.pdbqt"
done

## dock using VINA
# conda install -c conda-forge vina -y

# create one pocket-job table and reuse it for every amino acid
python - <<'PY'
import pandas as pd

x = pd.read_csv(
    "results/pockets_for_docking.tsv",
    sep="\t"
)

x[
    [
        "protein",
        "pocket",
        "center_x",
        "center_y",
        "center_z"
    ]
].to_csv(
    "results/vina_jobs.tsv",
    sep="\t",
    index=False,
    header=False
)

print("Jobs:", len(x))
PY

### run vina
# Run one amino acid at a time; its pockets run in parallel.
for ligand in ligands/amino_acids_pdbqt/*.pdbqt
do
    [ -f "$ligand" ] || continue
    aa=$(basename "$ligand" .pdbqt)
    [ "${aa^^}" = "DMET" ] && continue
    aa_lower=${aa,,}

    echo "Docking $aa"
    mkdir -p "docking_vina/${aa_lower}" "logs/vina_${aa_lower}"

    parallel -j 80 \
        --colsep '\t' \
        scripts/run_one_vina.sh {1} {2} {3} {4} {5} "$ligand" \
        :::: results/vina_jobs.tsv

    python scripts/collect_vina_results.py "$aa"
done

################################################################

################################################################

################################################################

################################################################

################################################################

################################################################

### 6. use OnionNet-SFCT for correction of docking scores
conda deactivate

# # a. Download OnionNet-SFCT
# cd /home/yoyerush/yo/whole_at_proteins_docking/ecoli_docking/tools
# git clone https://github.com/zhenglz/OnionNet-SFCT.git
# cd OnionNet-SFCT

# # b. Create a separate SFCT environment
# conda create -n sfct -c conda-forge python=3.8 openbabel=3.1.1 pip -y
conda activate sfct
# python -m pip install "pip<25"
# python -m pip install --only-binary=:all: numpy==1.19.5 scipy==1.5.4 pandas==1.1.5 scikit-learn==0.23.2
# 
# python -m pip install "Cython<3"
# python -m pip install --no-build-isolation mdtraj==1.9.7
# python -m pip install biopandas==0.2.9

# # c. Download the trained SFCT model
# mkdir -p data
# 
# pip install gdown
# gdown "https://drive.google.com/uc?id=1iiJvW4GBfg4D7LCuTRLKv9qnRYu5L2o5" -O data/sfct.model


# d. Run SFCT for every completed amino-acid docking
cd /home/yoyerush/yo/whole_at_proteins_docking/ecoli_docking

for ligand in ligands/amino_acids_pdbqt/*.pdbqt
do
    [ -f "$ligand" ] || continue
    aa=$(basename "$ligand" .pdbqt)
    [ "${aa^^}" = "DMET" ] && continue
    aa_lower=${aa,,}

    echo "Scoring $aa with SFCT"
    mkdir -p "docking_sfct/${aa_lower}" "logs/sfct_${aa_lower}"

    find "docking_vina/${aa_lower}" -name "*_${aa}.pdbqt" | \
    parallel -j 80 \
        --joblog "results/sfct_${aa_lower}_parallel.joblog" \
        'nice -n 10 scripts/run_one_sfct.sh {}'

    python scripts/collect_sfct_results.py "$aa"
    python scripts/compare_vina_sfct_results.py "$aa"
done

### 7. Create regular-AA thresholds and filtered hit tables
python scripts/filter_all_aa_results.py

# Export the L-AA bundle before starting the D-AA runs.
python scripts/prepare_atlas_data.py \
    --configuration L \
    --vina-dir results --sfct-dir results \
    --output results/compact/L --overwrite

################################################################

# Export Vina ligand pose coordinates
python scripts/export_vina_ligand_positions.py --bundle results/compact/L

# Export P2Rank pocket-point coordinates
python scripts/export_pocket_points.py \
    --bundle results/compact/L \
    --points-dir visualizations_p2rank
    
gzip results/compact/pocket_points.tsv

# not working # ## Export Vina pose diagnostics: H-bonds, salt bridges, clashes, and pocket proximity
# not working # # conda create -n pose_diag -c conda-forge python=3.11 rdkit "meeko>=0.6" "prolif>=2" numpy scipy -y
# not working # conda deactivate
# not working # conda activate pose_diag
# not working # python scripts/export_pose_diagnostics.py \
# not working #   --bundle results/compact/L \
# not working #   --uniprot-tsv metadata/ecoli_uniprot.tsv \
# not working #   --points-dir visualizations_p2rank \
# not working #   --jobs 8 2>&1 | tee logs/pose_diagnostics.log

## Export electrostatic potential values at the pocket points
# conda create -n electrostatics -c conda-forge python=3.11 numpy scipy rdkit "meeko>=0.6" pdb2pqr propka apbs ambertools -y
# python -m pip install requests
conda deactivate
conda activate electrostatics
python scripts/export_electrostatics.py \
  --bundle results/compact/L \
  --uniprot-tsv metadata/ecoli_uniprot.tsv \
  --points-dir visualizations_p2rank \
  --pdb2pqr pdb2pqr30 \
  --jobs 64 2>&1 | tee logs/electrostatics_potential.log

# before download to the app, also gzip the uniprot file
gzip metadata/ecoli_uniprot.tsv

################################################################

### 8. Run D-AA Vina docking and collect results
conda deactivate
conda activate docking

# Run D-AAs separately using the same pocket jobs and docking settings.
for ligand in ligands/d_amino_acids_pdbqt/*.pdbqt
do
    [ -f "$ligand" ] || continue
    aa=$(basename "$ligand" .pdbqt)
    aa_lower=${aa,,}

    echo "Docking $aa (D-AA set)"
    mkdir -p "docking_vina_d/${aa_lower}" "logs/vina_d/${aa_lower}"

    parallel -j 80 \
        --colsep '\t' \
        scripts/run_one_vina.sh {1} {2} {3} {4} {5} "$ligand" --d \
        :::: results/vina_jobs.tsv

    python scripts/collect_vina_results.py "$aa" --d
done

### 9. Run D-AA SFCT and collect results
conda deactivate
conda activate sfct

# Score the separate D-AA docking outputs.
mkdir -p results/d_amino_acids
for ligand in ligands/d_amino_acids_pdbqt/*.pdbqt
do
    [ -f "$ligand" ] || continue
    aa=$(basename "$ligand" .pdbqt)
    aa_lower=${aa,,}

    echo "Scoring $aa with SFCT (D-AA set)"
    mkdir -p "docking_sfct_d/${aa_lower}" "logs/sfct_d/${aa_lower}"

    find "docking_vina_d/${aa_lower}" -name "*_${aa}.pdbqt" | \
    parallel -j 80 \
        --joblog "results/d_amino_acids/sfct_${aa_lower}_parallel.joblog" \
        'nice -n 10 scripts/run_one_sfct.sh {} --d'

    python scripts/collect_sfct_results.py "$aa" --d
    python scripts/compare_vina_sfct_results.py "$aa" --d
done

################################################################

### 10. Create D-AA thresholds and filtered hit tables
python scripts/filter_all_aa_results.py --d

# Export only the D-AA results into their own bundle.
python scripts/prepare_atlas_data.py \
    --configuration D \
    --vina-d-dir results/d_amino_acids --sfct-d-dir results/d_amino_acids \
    --output results/compact/D --overwrite

conda deactivate
