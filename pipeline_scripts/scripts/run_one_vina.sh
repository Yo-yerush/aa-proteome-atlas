#!/bin/bash

protein="$1"
pocket="$2"
x="$3"
y="$4"
z="$5"
ligand="$6"

aa=$(basename "$ligand" .pdbqt)
aa_lower=${aa,,}

receptor="receptors_pdbqt/${protein}.pdbqt"
out="docking_vina/${aa_lower}/${protein}_${pocket}_${aa}.pdbqt"
log="logs/vina_${aa_lower}/${protein}_${pocket}_${aa}.log"

# Pass --d after the ligand path for the separate D-AA run.
if [ "${7:-}" = "--d" ]; then
    out="docking_vina_d/${aa_lower}/${protein}_${pocket}_${aa}.pdbqt"
    log="logs/vina_d/${aa_lower}/${protein}_${pocket}_${aa}.log"
fi

if [ ! -f "$receptor" ]; then
    echo "MISSING RECEPTOR: $receptor" > "$log"
    exit 1
fi

if [ ! -f "$ligand" ]; then
    echo "MISSING LIGAND: $ligand" > "$log"
    exit 1
fi

if [ -s "$out" ]; then
    exit 0
fi

vina \
    --receptor "$receptor" \
    --ligand "$ligand" \
    --center_x "$x" \
    --center_y "$y" \
    --center_z "$z" \
    --size_x 15 \
    --size_y 15 \
    --size_z 15 \
    --exhaustiveness 32 \
    --num_modes 10 \
    --cpu 1 \
    --out "$out" \
    > "$log" 2>&1
