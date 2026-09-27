#!/usr/bin/env python3

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
import csv

OUTDIR = Path("ligands/amino_acids_sdf")
OUTDIR.mkdir(parents=True, exist_ok=True)
D_OUTDIR = Path("ligands/d_amino_acids_sdf")
D_OUTDIR.mkdir(parents=True, exist_ok=True)

# L-amino acids, physiologically relevant zwitterionic forms
# His is represented with a neutral imidazole side chain.

AMINO_ACIDS = {
    "ALA": ("Alanine",       "A", "C[C@H]([NH3+])C(=O)[O-]"),
    "ARG": ("Arginine",      "R", "NC(=[NH2+])NCCC[C@H]([NH3+])C(=O)[O-]"),
    "ASN": ("Asparagine",    "N", "NC(=O)C[C@H]([NH3+])C(=O)[O-]"),
    "ASP": ("Aspartate",     "D", "[O-]C(=O)C[C@H]([NH3+])C(=O)[O-]"),
    "CYS": ("Cysteine",      "C", "SC[C@H]([NH3+])C(=O)[O-]"),
    "GLN": ("Glutamine",     "Q", "NC(=O)CC[C@H]([NH3+])C(=O)[O-]"),
    "GLU": ("Glutamate",     "E", "[O-]C(=O)CC[C@H]([NH3+])C(=O)[O-]"),
    "GLY": ("Glycine",       "G", "[NH3+]CC(=O)[O-]"),
    "HIS": ("Histidine",     "H", "[NH3+][C@@H](Cc1cnc[nH]1)C(=O)[O-]"),
    "ILE": ("Isoleucine",    "I", "CC[C@H](C)[C@H]([NH3+])C(=O)[O-]"),
    "LEU": ("Leucine",       "L", "CC(C)C[C@H]([NH3+])C(=O)[O-]"),
    "LYS": ("Lysine",        "K", "[NH3+]CCCC[C@H]([NH3+])C(=O)[O-]"),
    "MET": ("Methionine",    "M", "CSCC[C@H]([NH3+])C(=O)[O-]"),
    "PHE": ("Phenylalanine", "F", "c1ccccc1C[C@H]([NH3+])C(=O)[O-]"),
    "PRO": ("Proline",       "P", "[NH2+]1CCC[C@H]1C(=O)[O-]"),
    "SER": ("Serine",        "S", "OC[C@H]([NH3+])C(=O)[O-]"),
    "THR": ("Threonine",     "T", "C[C@@H](O)[C@H]([NH3+])C(=O)[O-]"),
    "TRP": ("Tryptophan",    "W", "c1ccc2[nH]cc(C[C@H]([NH3+])C(=O)[O-])c2c1"),
    "TYR": ("Tyrosine",      "Y", "Oc1ccc(C[C@H]([NH3+])C(=O)[O-])cc1"),
    "VAL": ("Valine",        "V", "CC(C)[C@H]([NH3+])C(=O)[O-]"),
}

D_AMINO_ACIDS = {
    "DALA": ("D-Alanine",       "dA", "C[C@@H]([NH3+])C(=O)[O-]"),
    "DARG": ("D-Arginine",      "dR", "NC(=[NH2+])NCCC[C@@H]([NH3+])C(=O)[O-]"),
    "DASN": ("D-Asparagine",    "dN", "NC(=O)C[C@@H]([NH3+])C(=O)[O-]"),
    "DASP": ("D-Aspartate",     "dD", "[O-]C(=O)C[C@@H]([NH3+])C(=O)[O-]"),
    "DCYS": ("D-Cysteine",      "dC", "SC[C@@H]([NH3+])C(=O)[O-]"),
    "DGLN": ("D-Glutamine",     "dQ", "NC(=O)CC[C@@H]([NH3+])C(=O)[O-]"),
    "DGLU": ("D-Glutamate",     "dE", "[O-]C(=O)CC[C@@H]([NH3+])C(=O)[O-]"),
    "DHIS": ("D-Histidine",     "dH", "[NH3+][C@H](Cc1cnc[nH]1)C(=O)[O-]"),
    "DILE": ("D-Isoleucine",    "dI", "CC[C@@H](C)[C@@H]([NH3+])C(=O)[O-]"),
    "DLEU": ("D-Leucine",       "dL", "CC(C)C[C@@H]([NH3+])C(=O)[O-]"),
    "DLYS": ("D-Lysine",        "dK", "[NH3+]CCCC[C@@H]([NH3+])C(=O)[O-]"),
    "DMET": ("D-Methionine",    "dM", "CSCC[C@@H]([NH3+])C(=O)[O-]"),
    "DPHE": ("D-Phenylalanine", "dF", "c1ccccc1C[C@@H]([NH3+])C(=O)[O-]"),
    "DPRO": ("D-Proline",       "dP", "[NH2+]1CCC[C@@H]1C(=O)[O-]"),
    "DSER": ("D-Serine",        "dS", "OC[C@@H]([NH3+])C(=O)[O-]"),
    "DTHR": ("D-Threonine",     "dT", "C[C@H](O)[C@@H]([NH3+])C(=O)[O-]"),
    "DTRP": ("D-Tryptophan",    "dW", "c1ccc2[nH]cc(C[C@@H]([NH3+])C(=O)[O-])c2c1"),
    "DTYR": ("D-Tyrosine",     "dY", "Oc1ccc(C[C@@H]([NH3+])C(=O)[O-])cc1"),
    "DVAL": ("D-Valine",        "dV", "CC(C)[C@@H]([NH3+])C(=O)[O-]"),
}

manifests = {OUTDIR: [], D_OUTDIR: []}

ligand_entries = (
    (code, details, outdir)
    for amino_acids, outdir in ((AMINO_ACIDS, OUTDIR), (D_AMINO_ACIDS, D_OUTDIR))
    for code, details in amino_acids.items()
)

for code, (name, one_letter, smiles), outdir in ligand_entries:

    print(f"Preparing {code}: {name}")

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise RuntimeError(f"Could not parse SMILES for {code}: {smiles}")

    mol = Chem.AddHs(mol)

    # Generate several conformers so the starting geometry is not arbitrary
    params = AllChem.ETKDGv3()
    params.randomSeed = 42

    conf_ids = AllChem.EmbedMultipleConfs(
        mol,
        numConfs=20,
        params=params
    )

    if not conf_ids:
        raise RuntimeError(f"3D embedding failed for {code}")

    energies = []

    # MMFF where possible, otherwise UFF
    if AllChem.MMFFHasAllMoleculeParams(mol):

        props = AllChem.MMFFGetMoleculeProperties(mol)

        for cid in conf_ids:
            ff = AllChem.MMFFGetMoleculeForceField(
                mol,
                props,
                confId=cid
            )

            ff.Minimize(maxIts=1000)
            energies.append((ff.CalcEnergy(), cid))

        method = "MMFF94"

    else:

        for cid in conf_ids:
            ff = AllChem.UFFGetMoleculeForceField(
                mol,
                confId=cid
            )

            ff.Minimize(maxIts=1000)
            energies.append((ff.CalcEnergy(), cid))

        method = "UFF"

    # Keep lowest-energy conformer
    best_energy, best_cid = min(energies)

    outfile = outdir / f"{code}.sdf"

    mol.SetProp("_Name", name)
    mol.SetProp("three_letter_code", code)
    mol.SetProp("one_letter_code", one_letter)
    mol.SetProp("input_smiles", smiles)
    mol.SetProp("formal_charge", str(Chem.GetFormalCharge(mol)))
    mol.SetProp("optimization_method", method)
    mol.SetProp("optimized_energy", f"{best_energy:.6f}")

    writer = Chem.SDWriter(str(outfile))
    writer.write(mol, confId=best_cid)
    writer.close()

    manifests[outdir].append({
        "code": code,
        "name": name,
        "one_letter": one_letter,
        "formal_charge": Chem.GetFormalCharge(mol),
        "smiles": smiles,
        "optimization": method,
        "energy": best_energy,
        "sdf_file": str(outfile)
    })


# Save a separate metadata table for each ligand folder
for outdir, manifest in manifests.items():
    manifest_file = outdir / "amino_acid_manifest.tsv"

    with manifest_file.open("w", newline="") as f:

        writer = csv.DictWriter(
            f,
            fieldnames=manifest[0].keys(),
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(manifest)


    print()
    print("Finished.")
    print(f"Created {len(manifest)} amino-acid SDF files.")
    print(f"Directory: {outdir}")
    print(f"Manifest: {manifest_file}")
