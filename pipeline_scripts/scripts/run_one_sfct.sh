#!/bin/bash

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export BLIS_NUM_THREADS=1

lig="$1"

base=$(basename "$lig" .pdbqt)
aa=${base##*_}
aa_lower=${aa,,}

# Example:
# AF-Q8LPL7-F1-model_v6_pocket1_ALA
# ->
# AF-Q8LPL7-F1-model_v6
protein=${base%%_pocket*}

receptor="structures_clean/${protein}.pdb"
out="docking_sfct/${aa_lower}/${base}.dat"
log="logs/sfct_${aa_lower}/${base}.log"
model="tools/OnionNet-SFCT/data/sfct.model"

# Pass --d after the docked ligand path for the separate D-AA run.
if [ "${2:-}" = "--d" ]; then
    out="docking_sfct_d/${aa_lower}/${base}.dat"
    log="logs/sfct_d/${aa_lower}/${base}.log"
fi

# Resume-safe
if [ -s "$out" ]; then
    exit 0
fi

if [ ! -s "$receptor" ]; then
    echo "MISSING RECEPTOR: $receptor" > "$log"
    exit 1
fi

python tools/OnionNet-SFCT/scorer.py \
    -r "$receptor" \
    -l "$lig" \
    -o "$out" \
    --model "$model" \
    -w 0.8 \
    --ncpus 1 \
    > "$log" 2>&1
