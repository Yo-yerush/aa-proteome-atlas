#!/usr/bin/env python3
"""Export Vina MODEL 1 diagnostics for an existing compact atlas bundle.

Requires Python >=3.10, RDKit, Meeko >=0.6, ProLIF >=2, NumPy and SciPy.
Run from the species project directory, for example:

    python scripts/export_pose_diagnostics.py --bundle results/compact/L \
        --points-dir visualizations_p2rank --jobs 8

Run separately with --bundle results/compact/D for D-AAs. Output defaults to
<bundle>/pose_diagnostics/pose_diagnostics_<aa>.tsv.gz.
L/D runs write into their respective bundle directories.
Existing per-AA output requires --overwrite; results are replaced, not merged.
Only the requested AA tables are written; no docking, minimization or SFCT is run.
Progress counts completed pose rows (including failures) across all workers and
updates about once per second while poses finish, with decimal percentages.

Terminal metadata (read locally; no downloads):
  --uniprot-tsv selects a UniProt table with Entry and Length columns. If omitted,
  use the sole metadata/*_uniprot.tsv(.gz) file. Optional Sequence and Fragment
  columns are also checked. Isoform accessions require their own length entry.
  Read the modeled range from UniProt DBREF/DBREF1/DBREF2 records in
  --structure-dir/<protein>.pdb (default structures_clean). Validate it against
  receptor residue numbering and the full UniProt length. Without DBREF, use an
  exact unique Sequence match, if available, or a complete chain of full length.
  For partial models lacking these records, --model-ranges accepts a TSV with:
      protein  chain  uniprot_start  uniprot_end
  One row per chain; inclusive UniProt positions; '.' for a blank chain ID.
  The range must cover the entire observed chain, with no internal gaps.
  Missing, ambiguous or inconsistent coverage gives coverage_error, not guessed
  charges. Fragment numbers such as F1/F2 are never used to infer coverage.

Chemistry:
  * Extract only a complete MODEL 1, including Meeko SMILES/index/H-parent
    remarks. Meeko restores the docked coordinates, including polar hydrogens.
    Transfer these onto the original SDF molecule through a full, chiral atom
    match, preserving its bonds, formal charges and stereochemistry.
  * Reconstruct the prepared receptor PDBQT with Meeko residue templates.
    Preserve existing heavy-atom coordinates and restore nonpolar Hs. Complete
    biological termini (UniProt positions 1/full length) as NH3+ (NH2+ for
    proline) and COO-. Artificial ends receive neutral acetyl/N-methylamide
    caps. Only terminal H/OXT atoms may be removed to establish these groups;
    all added terminal/cap coordinates are modeled, without minimization.
    Cap atoms are excluded from both ProLIF contacts and steric clashes.
    If any ligand heavy atom or retained pocket point is within
    --boundary-cutoff (default 6.0 A) of an artificial terminal residue's heavy
    atoms or cap, write unreliable_fragment_boundary with blank diagnostics.
    Reject non-terminal missing polar Hs and all remaining valence/radical errors.
    Evaluate ligand and receptor together in their unchanged coordinate frame
    with ProLIF; do not write intermediate complex files.

Definitions (distances in angstroms):
  h_bonds: distinct directed donor/acceptor heavy-atom pairs found by ProLIF
    HBDonor/HBAcceptor, default donor-acceptor distance <=3.5 and DHA 130-180 deg.
  salt_bridges: distinct oppositely charged group pairs found by ProLIF
    Cationic/Anionic, default distance <=4.5. Carboxylate O atoms and
    guanidinium/amidinium N atoms belonging to one group count together.
    A contact can count as both an H-bond and a salt bridge.
  clashes: ligand/receptor heavy-atom pairs whose RDKit van der Waals radii
    overlap by more than --clash-overlap (default 0.4). A geometric diagnostic,
    not a force-field energy; H atoms are excluded.
  pocket_proximity: fraction (0-1) of ligand heavy atoms within --point-cutoff
    (default 2.0) of any point assigned to this pocket's P2Rank rank.
  pocket_mean_distance: mean ligand-heavy-atom distance to its nearest pocket
    point. Neither proximity metric is a binding probability.

Raw *_points.pdb(.gz) files are matched by protein and pocket rank, avoiding
assumptions that pocket_id values match between L and D bundles. Rank 0 and
other pockets' points are excluded. Both aa/ and docking_aa/ Vina folders work.
Failed/unreliable rows have blank vina_pose and diagnostics, plus a status/error.
Exit codes: 0 all rows succeeded; 2 tables contain failed/unreliable rows;
1 export failed. Complete tables of failures are retained for diagnosis.

API references:
https://prolif.readthedocs.io/en/stable/notebooks/docking.html
https://prolif.readthedocs.io/en/stable/source/modules/interaction-fingerprint.html
https://meeko.readthedocs.io/en/develop/py_rec_prep.html
https://www.wwpdb.org/documentation/file-format-content/format33/sect3.html
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from contextlib import ExitStack
import csv
import gzip
from importlib.metadata import PackageNotFoundError, version
import io
import json
import math
from multiprocessing import get_context
import os
from pathlib import Path
import re
import sys
import tempfile
import time


L_AAS = "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
D_AAS = {"D" + aa for aa in L_AAS if aa != "GLY"}
RESIDUE_CODES = dict(zip(L_AAS, "ARNDCQEGHILKMFPSTWYV"))
FIELDS = (
    "pocket_id", "aa", "vina_pose", "h_bonds", "salt_bridges", "clashes",
    "pocket_proximity", "pocket_mean_distance", "status", "error",
)
# These are AutoDock atom types, not the PDB element column.
AD_ELEMENTS = {
    "H": "H", "HD": "H", "HS": "H", "C": "C", "A": "C",
    "N": "N", "NA": "N", "NS": "N", "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S", "P": "P", "F": "F", "Cl": "Cl", "Br": "Br",
    "I": "I", "Mg": "Mg", "Ca": "Ca", "Mn": "Mn", "Fe": "Fe", "Zn": "Zn",
    "Cu": "Cu", "Si": "Si", "B": "B",
}


class DiagnosticError(ValueError):
    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


def load_dependencies():
    # Delayed imports keep --help usable without the chemistry environment.
    global np, Chem, plf, cKDTree, Polymer, ResidueChemTemplates
    global PDBQTMolecule, RDKitMolCreate, find_inter_mols_bonds
    try:
        import numpy as np
        from rdkit import Chem
        import prolif as plf
        from scipy.spatial import cKDTree
        from meeko import Polymer, ResidueChemTemplates, PDBQTMolecule, RDKitMolCreate
        from meeko.polymer import find_inter_mols_bonds
    except ImportError as exc:
        raise ValueError(
            "Missing chemistry dependency. Use Python >=3.10 with rdkit, "
            "meeko>=0.6, prolif>=2, numpy and scipy installed. " + str(exc)
        ) from exc


def table_rows(path, required):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        headers = reader.fieldnames or []
        if not required.issubset(headers) or len(headers) != len(set(headers)):
            raise ValueError(f"{path}: missing required columns or duplicate headers")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}, line {reader.line_num}: invalid TSV row")
            if any(not row[key].strip() for key in required):
                raise ValueError(f"{path}, line {reader.line_num}: empty required field")
            yield row


def bundle_file(bundle, name):
    path = (bundle / name).resolve()
    if path.parent != bundle:
        raise ValueError(f"Expected a file directly inside the bundle: {name}")
    return path


def filename_part(value):
    if value in (".", "..") or any(c in value for c in "/\\\t\r\n\0"):
        raise ValueError(f"Invalid protein/pocket filename component: {value!r}")
    return value


def load_jobs(args):
    manifest = json.loads((args.bundle / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("format") != "aa-proteome-atlas-compact":
        raise ValueError("--bundle must be an existing compact atlas bundle")
    pockets = {}
    for row in table_rows(bundle_file(args.bundle, manifest["pockets"]["file"]),
                          {"pocket_id", "protein", "pocket", "rank"}):
        pid = row["pocket_id"]
        rank = float(row["rank"])
        if pid in pockets or not math.isfinite(rank) or rank < 1 or not rank.is_integer():
            raise ValueError(f"Duplicate pocket_id or invalid pocket rank: {pid}")
        pockets[pid] = (filename_part(row["protein"]), filename_part(row["pocket"]), int(rank))
    selected = set(args.aa or [])
    codes, groups = [], {}
    for entry in manifest["ligands"]:
        aa = entry["code"]
        if aa not in set(L_AAS) | D_AAS or aa in codes:
            raise ValueError(f"Unsupported or duplicate AA code: {aa}")
        codes.append(aa)
        if selected and aa not in selected:
            continue
        seen = set()
        for row in table_rows(bundle_file(args.bundle, entry["file"]), {"pocket_id"}):
            pid = row["pocket_id"]
            if pid not in pockets or pid in seen:
                raise ValueError(f"{aa}: unknown or repeated pocket_id {pid}")
            seen.add(pid)
            protein, pocket, rank = pockets[pid]
            groups.setdefault(protein, []).append((pid, aa, pocket, rank))
        if len(seen) != entry["rows"]:
            raise ValueError(f"{aa}: score-table row count differs from manifest")
    if selected - set(codes):
        raise ValueError(f"AA codes absent from this bundle: {sorted(selected - set(codes))}")
    if not groups:
        raise ValueError("No pocket/AA rows to process")
    return groups, [aa for aa in codes if not selected or aa in selected]


def alphafold_accession(protein):
    match = re.fullmatch(r"AF-([A-Z0-9]+(?:-\d+)?)-F\d+-model_v\d+", protein)
    return match[1] if match else None


def load_uniprot_metadata(args):
    source = args.uniprot_tsv
    if source is None:
        candidates = sorted(Path("metadata").glob("*_uniprot.tsv"))
        candidates += sorted(Path("metadata").glob("*_uniprot.tsv.gz"))
        if len(candidates) != 1:
            raise ValueError("Use --uniprot-tsv to select the table containing Entry and Length")
        source = candidates[0]
    records = {}
    for row in table_rows(source, {"Entry", "Length"}):
        accession, length = row["Entry"].strip(), int(row["Length"])
        sequence = row.get("Sequence", "").strip().upper()
        if length < 1 or (sequence and (len(sequence) != length
                                       or not re.fullmatch(r"[A-Z]+", sequence))):
            raise ValueError(f"Invalid UniProt length/sequence: {accession}")
        fragment = row.get("Fragment", "").strip().lower() not in ("", "no", "false", "0")
        record = (length, sequence, fragment)
        if accession in records and records[accession] != record:
            raise ValueError(f"Conflicting UniProt metadata: {accession}")
        records[accession] = record
    print(f"UniProt metadata: {source}", flush=True)
    return records


def load_model_ranges(path):
    """Optional explicit AlphaFold coverage if the clean PDB lacks DBREF records."""
    ranges = {}
    if path is not None:
        for row in table_rows(path, {"protein", "chain", "uniprot_start", "uniprot_end"}):
            key = (row["protein"], row["chain"] if row["chain"] != "." else "")
            start, end = int(row["uniprot_start"]), int(row["uniprot_end"])
            if key[1] in ranges.get(key[0], {}) or not 1 <= start <= end:
                raise ValueError(f"Duplicate or invalid modeled range: {key}")
            ranges.setdefault(key[0], {})[key[1]] = (start, end)
    return ranges


def point_files(root, proteins):
    if not root.is_dir():
        raise ValueError(f"Point directory does not exist: {root}")
    files = {}
    for path in root.rglob("*_points.pdb*"):
        suffix = next((s for s in ("_points.pdb.gz", "_points.pdb") if path.name.endswith(s)), None)
        if suffix is None or not path.is_file():
            continue
        stem = path.name[:-len(suffix)]
        protein = stem if stem in proteins else Path(stem).stem
        if protein not in proteins:
            continue
        if protein in files:
            raise ValueError(f"Multiple point files for {protein}: {files[protein]} and {path}")
        files[protein] = path
    if not files:
        raise ValueError(f"No point files matching bundle proteins under {root}")
    return files


def read_points(path, ranks):
    if path is None:
        raise DiagnosticError("missing_points", "No matching P2Rank point file")
    points = {rank: [] for rank in ranks}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            rank = int(line[22:26])
            if rank not in points:
                continue
            xyz = [float(line[start:start + 8]) for start in (30, 38, 46)]
            if not all(math.isfinite(value) for value in xyz):
                raise ValueError(f"Nonfinite point coordinates in {path}")
            points[rank].append(xyz)
    return {rank: cKDTree(xyz) for rank, xyz in points.items() if xyz}


def model_one(path):
    if not path.is_file():
        raise DiagnosticError("missing_pose", str(path))
    preamble, lines = [], []
    active, seen_model = False, False
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                seen_model = True
                if active:
                    break
                active = int(line.split()[1]) == 1
                if active:
                    lines = preamble + [line]
            elif active:
                lines.append(line)
                if line.startswith("ENDMDL"):
                    block = "".join(lines)
                    if "BEGIN_RES" in block:
                        raise ValueError("Flexible-receptor poses are not supported")
                    if not any(s.startswith(("ATOM", "HETATM")) for s in lines):
                        break
                    if "REMARK SMILES IDX " not in block or "REMARK SMILES " not in block:
                        raise ValueError("Missing Meeko SMILES/atom-mapping remarks")
                    return block
            elif line.startswith("REMARK") and not seen_model:
                preamble.append(line)
    raise DiagnosticError("incomplete_pose", f"No complete MODEL 1 in {path}")


def load_ligand(path):
    supplier = Chem.SDMolSupplier(str(path), removeHs=False)
    if len(supplier) != 1 or supplier[0] is None:
        raise ValueError(f"Expected one valid original AA molecule: {path}")
    mol = Chem.AddHs(supplier[0])
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError(f"Disconnected original ligand: {path}")
    return mol


def reconstruct_ligand(path, original):
    block = model_one(path)
    molecules = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(block, skip_typing=True))
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError("Meeko could not reconstruct exactly one ligand")
    docked = molecules[0]
    if docked.GetNumConformers() != 1:
        raise ValueError("Expected exactly one MODEL 1 conformer")
    if Chem.MolToSmiles(Chem.RemoveHs(docked)) != Chem.MolToSmiles(Chem.RemoveHs(original)):
        raise ValueError("Docked ligand chemistry/stereochemistry differs from original SDF")
    mapping = docked.GetSubstructMatch(original, useChirality=True)
    if len(mapping) != original.GetNumAtoms() or docked.GetNumAtoms() != original.GetNumAtoms():
        raise ValueError("Cannot map every original SDF atom to the docked molecule")
    ligand = Chem.Mol(original)
    ligand.RemoveAllConformers()
    conf = Chem.Conformer(ligand.GetNumAtoms())
    conf.Set3D(True)
    pose_conf = docked.GetConformer()
    for index, docked_index in enumerate(mapping):
        conf.SetAtomPosition(index, pose_conf.GetAtomPosition(docked_index))
    ligand.AddConformer(conf)
    if not np.isfinite(conf.GetPositions()).all():
        raise ValueError("Nonfinite ligand coordinates")
    # Prove that heavy atoms and polar Hs came from the pose, not generated Hs.
    source_coords = Counter()
    for line in block.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            element = AD_ELEMENTS[line[77:].strip()]
            xyz = tuple(round(float(line[start:start + 8]), 3) for start in (30, 38, 46))
            source_coords[(element, *xyz)] += 1
    restored_coords, required_coords = Counter(), Counter()
    for atom in ligand.GetAtoms():
        xyz = tuple(round(value, 3) for value in conf.GetAtomPosition(atom.GetIdx()))
        key = (atom.GetSymbol(), *xyz)
        restored_coords[key] += 1
        if atom.GetAtomicNum() != 1 or any(n.GetAtomicNum() != 6 for n in atom.GetNeighbors()):
            required_coords[key] += 1
    if source_coords - restored_coords or required_coords - source_coords:
        raise ValueError("MODEL 1 atom coordinates or polar hydrogen mapping are incomplete")
    return ligand


def receptor_pdb(path):
    if not path.is_file():
        raise DiagnosticError("missing_receptor", str(path))
    lines, seen, count = [], set(), 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(("MODEL", "BEGIN_RES")):
                raise ValueError("Expected a single rigid receptor PDBQT")
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_type = line[77:].strip()
            if atom_type not in AD_ELEMENTS:
                raise ValueError(f"Unsupported receptor AutoDock atom type: {atom_type}")
            key = (line[21:27], line[12:16])
            if key in seen or line[16:17].strip():
                raise ValueError("Duplicate receptor atom or unresolved alternate location")
            seen.add(key)
            xyz = [float(line[start:start + 8]) for start in (30, 38, 46)]
            if not all(math.isfinite(value) for value in xyz):
                raise ValueError("Nonfinite receptor coordinates")
            lines.append(line[:66].ljust(76) + f"{AD_ELEMENTS[atom_type]:>2}\n")
            count += 1
    if not count:
        raise ValueError("Receptor has no atoms")
    return "".join(lines) + "END\n", count


def pdb_reference_ranges(path, accession):
    """Read original model-to-UniProt mappings, including long-accession DBREF1/2."""
    references, pending = {}, {}
    if not path.is_file():
        return references
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            kind = line[:6].strip()
            if kind in ("ATOM", "HETATM", "MODEL"):
                break
            if kind not in ("DBREF", "DBREF1", "DBREF2"):
                continue
            line = line.rstrip("\r\n").ljust(80)
            chain = line[12].strip()
            if kind in ("DBREF", "DBREF1"):
                if line[26:32].strip() != "UNP":
                    continue
                if line[18].strip() or line[24].strip():
                    raise DiagnosticError("coverage_error", f"{path}: insertion codes in DBREF")
                bounds = (int(line[14:18]), int(line[20:24]))
                if kind == "DBREF1":
                    if chain in pending:
                        raise DiagnosticError("coverage_error", f"{path}: ambiguous DBREF1")
                    pending[chain] = bounds
                    continue
                entry = line[33:41].strip()
                start, end = int(line[55:60]), int(line[62:67])
            else:
                if chain not in pending:
                    continue
                bounds = pending.pop(chain)
                entry = line[18:40].strip()
                start, end = int(line[45:55]), int(line[57:67])
            if entry != accession or chain in references:
                raise DiagnosticError("coverage_error", f"{path}: mismatched/ambiguous UniProt DBREF")
            references[chain] = (*bounds, start, end)
    if pending:
        raise DiagnosticError("coverage_error", f"{path}: incomplete DBREF1/DBREF2")
    return references


def terminal_backbone_atoms(monomers, inter_bonds, protein, metadata, model_ranges):
    """Classify free ends only after establishing this model's UniProt coverage."""
    accession = alphafold_accession(protein)
    if accession is None or metadata is None:
        raise DiagnosticError("coverage_error", f"{protein}: missing exact UniProt accession/Length")
    full_length, full_sequence, is_fragment = metadata
    if is_fragment:
        raise DiagnosticError("coverage_error", f"{accession}: UniProt entry is itself a fragment")
    references = pdb_reference_ranges(OPTIONS.structure_dir / f"{protein}.pdb", accession)
    linked = set()
    connected = set()
    for (left, right), bonds in inter_bonds.items():
        for i, j in bonds:
            a, b = (left, monomers[left].mapidx_from_raw[i]), (right, monomers[right].mapidx_from_raw[j])
            linked.update((a, b))
            connected.add(frozenset((a, b)))
    chains = {}
    for resid, monomer in monomers.items():
        if not {"N", "CA", "C", "O"}.issubset(monomer.atom_names):
            continue
        chain, number = resid.rsplit(":", 1)
        chain = chain.strip()
        insertion = number[-1] if number[-1].isalpha() else ""
        sequence = int(number[:-1] if insertion else number)
        chains.setdefault(chain, []).append((sequence, insertion, resid))
    if not chains:
        raise DiagnosticError("coverage_error", f"{protein}: no protein chains")
    n_termini, c_termini = {}, {}
    for chain, residues in chains.items():
        residues.sort()
        first, last = residues[0][2], residues[-1][2]
        numbers = [number for number, insertion, resid in residues]
        if (any(insertion for number, insertion, resid in residues)
                or numbers != list(range(numbers[0], numbers[-1] + 1))):
            raise DiagnosticError("coverage_error", f"{protein}, chain {chain}: noncontiguous residues")
        # Internal missing peptide links must remain failures, never new termini.
        for left, right in zip(residues, residues[1:]):
            a = (left[2], monomers[left[2]].atom_names.index("C"))
            b = (right[2], monomers[right[2]].atom_names.index("N"))
            if frozenset((a, b)) not in connected:
                raise ValueError(f"Internal peptide break between {left[2]} and {right[2]}")
        modeled_sequence = "".join(RESIDUE_CODES.get(monomers[r[2]].input_resname, "?") for r in residues)
        bounds = model_ranges.get(chain)
        if chain in references:
            pdb_start, pdb_end, start, end = references[chain]
            if (pdb_start, pdb_end) != (numbers[0], numbers[-1]):
                raise DiagnosticError("coverage_error", f"{protein}, chain {chain}: DBREF/receptor mismatch")
            if bounds is not None and bounds != (start, end):
                raise DiagnosticError("coverage_error", f"{protein}, chain {chain}: conflicting modeled ranges")
            bounds = (start, end)
        if bounds is None and full_sequence:
            offset = full_sequence.find(modeled_sequence)
            if offset >= 0 and full_sequence.find(modeled_sequence, offset + 1) < 0:
                bounds = (offset + 1, offset + len(residues))
        if bounds is None and len(residues) == full_length:
            # A complete contiguous chain covers the full entry regardless of PDB numbering.
            bounds = (1, full_length)
        if bounds is None:
            raise DiagnosticError("coverage_error", f"{protein}, chain {chain}: missing modeled range; "
                                  "supply --model-ranges or a clean PDB with UniProt DBREF")
        start, end = bounds
        if not 1 <= start <= end <= full_length or end - start + 1 != len(residues):
            raise DiagnosticError("coverage_error", f"{protein}, chain {chain}: range/Length mismatch")
        if full_sequence and full_sequence[start - 1:end] != modeled_sequence:
            raise DiagnosticError("coverage_error", f"{protein}, chain {chain}: UniProt sequence mismatch")
        for resid, name, termini, biological in ((first, "N", n_termini, start == 1),
                                                 (last, "C", c_termini, end == full_length)):
            key = (resid, monomers[resid].atom_names.index(name))
            if key in linked:
                raise ValueError(f"Unexpected inter-residue bond at chain end {resid} {name}")
            termini[key] = biological
    return n_termini, c_termini


def unit_vector(vector):
    length = np.linalg.norm(vector)
    if not math.isfinite(float(length)) or length < 1e-6:
        raise ValueError("Degenerate terminal geometry")
    return vector / length


def remaining_bond_direction(center, neighbors):
    vectors = np.array([tuple(point) for point in neighbors]) - center
    return unit_vector(-sum(unit_vector(vector) for vector in vectors))


def planar_branches(center, bonded, plane_point):
    axis = unit_vector(np.asarray(tuple(bonded)) - center)
    reference = np.asarray(tuple(plane_point)) - center
    perpendicular = unit_vector(reference - np.dot(reference, axis) * axis)
    return (-0.5 * axis + math.sqrt(0.75) * perpendicular,
            -0.5 * axis - math.sqrt(0.75) * perpendicular)


def add_cap_atom(receptor, coords, parent, element, name, xyz, hydrogens=0,
                 order=None):
    """Keep caps in the same residue for ProLIF typing, but tag all cap contacts out."""
    atom = Chem.Atom(element)
    atom.SetNoImplicit(True)
    atom.SetNumExplicitHs(hydrogens)
    atom.SetBoolProp("diagnostics_cap", True)
    atom.SetBoolProp("diagnostics_boundary", True)
    if hydrogens:
        atom.SetBoolProp("diagnostics_add_h", True)
    info = receptor.GetAtomWithIdx(parent).GetPDBResidueInfo()
    atom.SetPDBResidueInfo(Chem.AtomPDBResidueInfo(
        atomName=name, residueName=info.GetResidueName(),
        residueNumber=info.GetResidueNumber(), chainId=info.GetChainId(),
    ))
    index = receptor.AddAtom(atom)
    receptor.AddBond(parent, index, order if order is not None else Chem.BondType.SINGLE)
    coords.append(tuple(xyz))
    return index


def complete_receptor_termini(receptor, coords, monomers, offsets, n_termini, c_termini):
    """Charge biological ends; acetylate/amidate artificial ends before sanitization."""
    remove_atoms = set()
    for (resid, local_index), biological in sorted(n_termini.items()):
        index = offsets[resid] + local_index
        nitrogen = receptor.GetAtomWithIdx(index)
        heavy = [a for a in nitrogen.GetNeighbors() if a.GetAtomicNum() > 1]
        proline = monomers[resid].input_resname == "PRO"
        expected = {"CA", "CD"} if proline else {"CA"}
        names = {a.GetPDBResidueInfo().GetName().strip() for a in heavy}
        if (nitrogen.GetAtomicNum() != 7 or names != expected or len(heavy) != len(expected)
                or any(a.GetAtomicNum() != 6 for a in heavy)
                or any(b.GetBondType() != Chem.BondType.SINGLE for b in nitrogen.GetBonds())
                or nitrogen.GetDegree() > 4):
            raise ValueError(f"Unsupported N-terminal connectivity in {resid}")
        nitrogen.SetFormalCharge(1 if biological else 0)
        nitrogen.SetNumRadicalElectrons(0)
        nitrogen.SetNoImplicit(True)
        nitrogen.SetBoolProp("diagnostics_add_h", True)
        if biological:
            # Existing explicit H atoms are kept; AddHs places only missing Hs.
            nitrogen.SetNumExplicitHs(4 - nitrogen.GetDegree())
            continue
        # An N-acetyl cap gives a neutral peptide N, not a neutral amine.
        hydrogens = sorted(a.GetIdx() for a in nitrogen.GetNeighbors() if a.GetAtomicNum() == 1)
        keep = hydrogens[:0 if proline else 1]
        remove_atoms.update(set(hydrogens) - set(keep))
        nitrogen.SetNumExplicitHs((0 if proline else 1) - len(keep))
        center = np.asarray(tuple(coords[index]))
        ca = offsets[resid] + monomers[resid].atom_names.index("CA")
        neighbors = [coords[a.GetIdx()] for a in heavy] + [coords[i] for i in keep]
        if len(neighbors) == 2:
            direction = remaining_bond_direction(center, neighbors)
        else:
            carbonyl = offsets[resid] + monomers[resid].atom_names.index("C")
            direction = planar_branches(center, coords[ca], coords[carbonyl])[0]
        cap_xyz = center + 1.33 * direction
        carbon = add_cap_atom(receptor, coords, index, 6, "AC", cap_xyz)
        oxygen_dir, methyl_dir = planar_branches(cap_xyz, center, coords[ca])
        add_cap_atom(receptor, coords, carbon, 8, "AO", cap_xyz + 1.23 * oxygen_dir,
                     order=Chem.BondType.DOUBLE)
        add_cap_atom(receptor, coords, carbon, 6, "AM", cap_xyz + 1.50 * methyl_dir,
                     hydrogens=3)

    for (resid, local_index), biological in sorted(c_termini.items()):
        index = offsets[resid] + local_index
        carbon = receptor.GetAtomWithIdx(index)
        ca = offsets[resid] + monomers[resid].atom_names.index("CA")
        oxygens = [a.GetIdx() for a in carbon.GetNeighbors() if a.GetAtomicNum() == 8]
        neighbors = {a.GetIdx() for a in carbon.GetNeighbors()}
        if (carbon.GetAtomicNum() != 6 or len(oxygens) not in (1, 2)
                or neighbors != {ca, *oxygens}
                or receptor.GetBondBetweenAtoms(index, ca).GetBondType() != Chem.BondType.SINGLE):
            raise ValueError(f"Unsupported C-terminal connectivity in {resid}")
        oxygen = offsets[resid] + monomers[resid].atom_names.index("O")
        if oxygen not in oxygens:
            raise ValueError(f"Missing backbone carbonyl O in {resid}")
        for oi in oxygens:
            for neighbor in receptor.GetAtomWithIdx(oi).GetNeighbors():
                if neighbor.GetIdx() == index:
                    continue
                if neighbor.GetAtomicNum() != 1 or neighbor.GetDegree() != 1:
                    raise ValueError(f"Unsupported terminal oxygen connectivity in {resid}")
                remove_atoms.add(neighbor.GetIdx())
        center = np.asarray(tuple(coords[index]))
        oxt = next((oi for oi in oxygens if oi != oxygen), None)
        if biological:
            if oxt is None:
                direction = remaining_bond_direction(center, [coords[ca], coords[oxygen]])
                new_oxygen = Chem.Atom(8)
                info = carbon.GetPDBResidueInfo()
                new_oxygen.SetPDBResidueInfo(Chem.AtomPDBResidueInfo(
                    atomName="OXT", residueName=info.GetResidueName(),
                    residueNumber=info.GetResidueNumber(), chainId=info.GetChainId(),
                ))
                oxt = receptor.AddAtom(new_oxygen)
                receptor.AddBond(index, oxt, Chem.BondType.SINGLE)
                coords.append(tuple(center + 1.25 * direction))
            oxygen_groups = ((oxygen, 0, Chem.BondType.DOUBLE), (oxt, -1, Chem.BondType.SINGLE))
        else:
            # Replace terminal OXT (if present) with a neutral N-methylamide cap.
            if oxt is not None:
                remove_atoms.add(oxt)
            direction = remaining_bond_direction(center, [coords[ca], coords[oxygen]])
            cap_xyz = center + 1.33 * direction
            cap_n = add_cap_atom(receptor, coords, index, 7, "NM", cap_xyz, hydrogens=1)
            methyl_dir, _ = planar_branches(cap_xyz, center, coords[oxygen])
            add_cap_atom(receptor, coords, cap_n, 6, "CM", cap_xyz + 1.45 * methyl_dir,
                         hydrogens=3)
            oxygen_groups = ((oxygen, 0, Chem.BondType.DOUBLE),)
        for oi, charge, order in oxygen_groups:
            atom = receptor.GetAtomWithIdx(oi)
            atom.SetFormalCharge(charge)
            atom.SetNumRadicalElectrons(0)
            atom.SetNumExplicitHs(0)
            atom.SetNoImplicit(True)
            receptor.GetBondBetweenAtoms(index, oi).SetBondType(order)
        carbon = receptor.GetAtomWithIdx(index)
        carbon.SetFormalCharge(0)
        carbon.SetNumRadicalElectrons(0)
        carbon.SetNumExplicitHs(0)
        carbon.SetNoImplicit(True)

    # Only terminal atoms are removed. Properties survive the index changes.
    for index in sorted(remove_atoms, reverse=True):
        receptor.RemoveAtom(index)
        coords.pop(index)
    return [atom.GetIdx() for atom in receptor.GetAtoms() if atom.HasProp("diagnostics_add_h")]


def reconstruct_receptor(path, protein, metadata, model_ranges):
    pdb, input_count = receptor_pdb(path)
    polymer = Polymer.from_pdb_string(pdb, CHEM_TEMPLATES, None, allow_bad_res=False)
    monomers = polymer.get_valid_monomers()
    if not monomers or len(monomers) != len(polymer.monomers):
        raise ValueError("Not every receptor residue has valid chemistry")
    raw_mols = {key: (m.raw_rdkit_mol, m.input_resname) for key, m in monomers.items()}
    inter_bonds = find_inter_mols_bonds(raw_mols)
    n_termini, c_termini = terminal_backbone_atoms(
        monomers, inter_bonds, protein, metadata, model_ranges,
    )
    artificial_residues = {resid for termini in (n_termini, c_termini)
                           for (resid, index), biological in termini.items() if not biological}
    receptor = Chem.RWMol()
    coords, offsets, matched = [], {}, 0
    for residue_number, (resid, monomer) in enumerate(monomers.items(), 1):
        mol = monomer.rdkit_mol
        offsets[resid] = receptor.GetNumAtoms()
        original_indices = set(monomer.mapidx_to_raw.values())
        if original_indices != set(range(monomer.raw_rdkit_mol.GetNumAtoms())):
            raise ValueError(f"Receptor template discards input atoms in {resid}")
        matched += len(original_indices)
        conf = mol.GetConformer()
        raw_conf = monomer.raw_rdkit_mol.GetConformer()
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()
            if idx not in monomer.mapidx_to_raw:
                terminal_h = (atom.GetAtomicNum() == 1 and atom.GetDegree() == 1
                              and (resid, atom.GetNeighbors()[0].GetIdx()) in n_termini)
                if not terminal_h and (atom.GetAtomicNum() != 1
                                       or any(n.GetAtomicNum() != 6 for n in atom.GetNeighbors())):
                    raise ValueError(f"Missing heavy atom or polar H in receptor residue {resid}")
            else:
                # Preserve every atom actually present in the prepared receptor.
                conf.SetAtomPosition(idx, raw_conf.GetAtomPosition(monomer.mapidx_to_raw[idx]))
            copied = Chem.Atom(atom)
            copied.SetNoImplicit(True)  # Residue-link placeholders are replaced below.
            if resid in artificial_residues:
                copied.SetBoolProp("diagnostics_boundary", True)
            copied.SetPDBResidueInfo(Chem.AtomPDBResidueInfo(
                atomName=monomer.atom_names[idx], residueName=monomer.input_resname,
                residueNumber=residue_number, chainId="R",
            ))
            receptor.AddAtom(copied)
            coords.append(conf.GetAtomPosition(idx))
        for bond in mol.GetBonds():
            receptor.AddBond(offsets[resid] + bond.GetBeginAtomIdx(),
                             offsets[resid] + bond.GetEndAtomIdx(), bond.GetBondType())
    if matched != input_count:
        raise ValueError("Receptor reconstruction did not preserve every input atom")
    # Reuse Meeko's inter-residue connectivity, including peptide/disulfide bonds.
    for (left, right), bonds in inter_bonds.items():
        for i, j in bonds:
            a = offsets[left] + monomers[left].mapidx_from_raw[i]
            b = offsets[right] + monomers[right].mapidx_from_raw[j]
            receptor.AddBond(a, b, Chem.BondType.SINGLE)
    add_hydrogens = complete_receptor_termini(
        receptor, coords, monomers, offsets, n_termini, c_termini,
    )
    mol = receptor.GetMol()
    Chem.SanitizeMol(mol)
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for idx, xyz in enumerate(coords):
        conf.SetAtomPosition(idx, xyz)
    mol.AddConformer(conf)
    if add_hydrogens:
        original_count = mol.GetNumAtoms()
        mol = Chem.AddHs(mol, explicitOnly=True, addCoords=True,
                         onlyOnAtoms=add_hydrogens, addResidueInfo=True)
        for index in range(original_count, mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(index)
            if any(n.HasProp("diagnostics_boundary") for n in atom.GetNeighbors()):
                atom.SetBoolProp("diagnostics_cap", True)
        Chem.SanitizeMol(mol)

    # if any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
    #     raise ValueError("Unresolved receptor valence after template reconstruction")
    radicals = []
    for atom in mol.GetAtoms():
        if atom.GetNumRadicalElectrons():
            info = atom.GetPDBResidueInfo()
            radicals.append((
                atom.GetIdx(),
                atom.GetSymbol(),
                info.GetResidueName() if info else "?",
                info.GetName().strip() if info else "?",
                info.GetResidueNumber() if info else -1,
                atom.GetDegree(),
                atom.GetFormalCharge(),
                atom.GetNumRadicalElectrons(),
            ))

    if radicals:
        print("RADICAL ATOMS:", radicals[:30], file=sys.stderr)
        raise ValueError("Unresolved receptor valence after template reconstruction")
    
    if not np.isfinite(mol.GetConformer().GetPositions()).all():
        raise ValueError("Nonfinite reconstructed receptor coordinates")
    return mol


def charged_groups(mol):
    """Group resonance-equivalent O/N atoms; other atoms remain separate."""
    groups = list(range(mol.GetNumAtoms()))
    for center in mol.GetAtoms():
        if center.GetAtomicNum() not in (6, 15, 16):
            continue
        for element, sign in ((8, -1), (7, 1)):
            atoms = [a for a in center.GetNeighbors() if a.GetAtomicNum() == element]
            if not any(a.GetFormalCharge() * sign > 0 for a in atoms):
                continue
            if not any(mol.GetBondBetweenAtoms(center.GetIdx(), a.GetIdx()).GetBondType()
                       == Chem.BondType.DOUBLE for a in atoms):
                continue
            group = min(a.GetIdx() for a in atoms)
            for atom in atoms:
                groups[atom.GetIdx()] = group
    return groups


def heavy_geometry(mol):
    atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() > 1 and not a.HasProp("diagnostics_cap")]
    coords = mol.GetConformer().GetPositions()[[a.GetIdx() for a in atoms]]
    table = Chem.GetPeriodicTable()
    radii = np.array([table.GetRvdw(a.GetAtomicNum()) for a in atoms])
    if not len(atoms) or not np.isfinite(coords).all() or not (radii > 0).all():
        raise ValueError("Invalid heavy-atom coordinates or van der Waals radii")
    return coords, radii


def count_interactions(ligand, receptor, receptor_groups, cap_atoms):
    # ProLIF uses a ligand/receptor pair in one coordinate frame as the complex.
    lig = plf.Molecule.from_rdkit(ligand)
    contacts = FINGERPRINT.generate(lig, receptor, metadata=True)
    ligand_groups = charged_groups(ligand)
    hbonds, salts = set(), set()
    for interactions in contacts.values():
        for kind, matches in interactions.items():
            for match in matches:
                indices = match["parent_indices"]
                if cap_atoms.intersection(indices["protein"]):
                    continue
                la = next(i for i in indices["ligand"] if ligand.GetAtomWithIdx(i).GetAtomicNum() > 1)
                ra = next(i for i in indices["protein"] if receptor.GetAtomWithIdx(i).GetAtomicNum() > 1)
                if kind in ("HBDonor", "HBAcceptor"):
                    hbonds.add((kind, la, ra))
                else:
                    salts.add((kind, ligand_groups[la], receptor_groups[ra]))
    return len(hbonds), len(salts)


def count_clashes(lig_xyz, lig_radii, rec_tree, rec_radii, overlap):
    total = 0
    for xyz, radius in zip(lig_xyz, lig_radii):
        search_radius = float(radius + rec_radii.max() - overlap)
        if search_radius <= 0:
            continue
        indices = rec_tree.query_ball_point(xyz, search_radius)
        if indices:
            distances = np.linalg.norm(rec_tree.data[indices] - xyz, axis=1)
            total += int(np.count_nonzero(radius + rec_radii[indices] - distances > overlap))
    return total


class PoseProgress:
    """Count finished rows independently of whole-protein result delivery."""

    def __init__(self, total):
        self.total = total
        self.context = get_context()
        self.counter = self.context.Value("Q", 0)
        self.last_completed = 0
        self.last_update = time.monotonic()

    def report(self, force=False):
        now = time.monotonic()
        if not force and now - self.last_update < 1.0:
            return
        with self.counter.get_lock():
            completed = self.counter.value
        if completed == self.last_completed:
            return
        # Truncate so rounding cannot display 100% before the final row finishes.
        percent = (completed * 100_000 // self.total) / 1000
        print(f"Diagnostics: {percent:.3f}% ({completed:,}/{self.total:,})", flush=True)
        self.last_completed = completed
        self.last_update = now


def record_pose_progress(count, progress):
    with COMPLETED_POSES.get_lock():
        COMPLETED_POSES.value += count
    if progress is not None:  # The single-worker path reports from the main process.
        progress.report()


def initialize_worker(options, codes, completed_poses=None):
    global OPTIONS, LIGANDS, CHEM_TEMPLATES, FINGERPRINT, COMPLETED_POSES
    load_dependencies()
    OPTIONS = options
    COMPLETED_POSES = completed_poses
    LIGANDS = {}
    for aa in codes:
        root = options.ligand_d_dir if aa in D_AAS else options.ligand_dir
        LIGANDS[aa] = load_ligand(root / f"{aa}.sdf")
    CHEM_TEMPLATES = ResidueChemTemplates.create_from_defaults()
    FINGERPRINT = plf.Fingerprint(
        ["HBDonor", "HBAcceptor", "Cationic", "Anionic"], count=True,
        parameters={
            "HBDonor": {"distance": options.hbond_distance, "DHA_angle": (options.hbond_angle, 180)},
            "HBAcceptor": {"distance": options.hbond_distance, "DHA_angle": (options.hbond_angle, 180)},
            "Cationic": {"distance": options.salt_distance},
            "Anionic": {"distance": options.salt_distance},
        },
        vicinity_cutoff=max(6.0, options.hbond_distance + 1, options.salt_distance + 1),
    )


def error_row(job, exc, stage):
    status = exc.status if isinstance(exc, DiagnosticError) else stage + "_error"
    message = " ".join(str(exc).split())[:300]
    return (job[0], job[1], "", "", "", "", "", "", status, message)


def process_protein(task, progress=None):
    protein, jobs, points_path, metadata, model_ranges = task
    stage = "receptor"
    try:
        rec = reconstruct_receptor(OPTIONS.receptor_dir / f"{protein}.pdbqt",
                                   protein, metadata, model_ranges)
        cap_atoms = {a.GetIdx() for a in rec.GetAtoms() if a.HasProp("diagnostics_cap")}
        boundary_atoms = [a.GetIdx() for a in rec.GetAtoms()
                          if a.GetAtomicNum() > 1 and a.HasProp("diagnostics_boundary")]
        boundary_tree = (cKDTree(rec.GetConformer().GetPositions()[boundary_atoms])
                         if boundary_atoms else None)
        rec_xyz, rec_radii = heavy_geometry(rec)
        rec_tree = cKDTree(rec_xyz)
        rec_groups = charged_groups(rec)
        receptor = plf.Molecule.from_rdkit(rec)
        stage = "points"
        points = read_points(points_path, {job[3] for job in jobs})
        boundary_distances = ({rank: float(boundary_tree.query(tree.data)[0].min())
                               for rank, tree in points.items()} if boundary_tree is not None else {})
    except Exception as exc:
        rows = [error_row(job, exc, stage) for job in jobs]
        record_pose_progress(len(rows), progress)
        return rows
    rows = []
    for job in jobs:
        pid, aa, pocket, rank = job
        stage = "points"
        try:
            if rank not in points:
                raise DiagnosticError("missing_points", f"{protein}: no points for {pocket}, rank {rank}")
            stage = "ligand"
            root = OPTIONS.vina_d_dir if aa in D_AAS else OPTIONS.vina_dir
            name = f"{protein}_{pocket}_{aa}.pdbqt"
            path = root / aa.lower() / name
            if not path.is_file():
                path = root / f"docking_{aa.lower()}" / name
            ligand = reconstruct_ligand(path, LIGANDS[aa])
            lig_xyz, lig_radii = heavy_geometry(ligand)
            if boundary_tree is not None:
                ligand_distance = float(boundary_tree.query(lig_xyz)[0].min())
                pocket_distance = boundary_distances[rank]
                if min(ligand_distance, pocket_distance) <= OPTIONS.boundary_cutoff:
                    raise DiagnosticError(
                        "unreliable_fragment_boundary",
                        f"Artificial chain end/cap within {OPTIONS.boundary_cutoff:g} A: "
                        f"ligand={ligand_distance:.2f} A, pocket points={pocket_distance:.2f} A",
                    )
            stage = "interaction"
            hbonds, salts = count_interactions(ligand, receptor, rec_groups, cap_atoms)
            clashes = count_clashes(lig_xyz, lig_radii, rec_tree, rec_radii, OPTIONS.clash_overlap)
            distances, _ = points[rank].query(lig_xyz)
            proximity = float(np.mean(distances <= OPTIONS.point_cutoff))
            rows.append((pid, aa, 1, hbonds, salts, clashes, f"{proximity:.4f}",
                         f"{float(distances.mean()):.4f}", "success", ""))
        except Exception as exc:
            rows.append(error_row(job, exc, stage))
        record_pose_progress(1, progress)
    return rows


def protein_results(tasks, args, codes, progress):
    if args.jobs == 1:
        global COMPLETED_POSES
        COMPLETED_POSES = progress.counter
        for task in tasks:
            yield process_protein(task, progress)
        progress.report(force=True)
        return

    with ProcessPoolExecutor(
        max_workers=args.jobs,
        mp_context=progress.context,
        initializer=initialize_worker,
        initargs=(args, codes, progress.counter),
    ) as pool:

        tasks = iter(tasks)
        pending = set()

        # Keep up to 2 x jobs tasks queued
        for _ in range(2 * args.jobs):
            task = next(tasks, None)
            if task is None:
                break
            pending.add(pool.submit(process_protein, task))

        while pending:
            done, pending = wait(
                pending,
                timeout=1.0,
                return_when=FIRST_COMPLETED,
            )
            progress.report()

            for future in done:
                yield future.result()

                task = next(tasks, None)
                if task is not None:
                    pending.add(
                        pool.submit(process_protein, task)
                    )
    progress.report(force=True)


def export(args):
    args.bundle = args.bundle.resolve()
    output_dir = (args.output_dir or args.bundle / "pose_diagnostics").resolve()
    if output_dir.exists() and not output_dir.is_dir():
        raise ValueError(f"Output directory is not a directory: {output_dir}")
    if not output_dir.parent.is_dir():
        raise ValueError(f"Output parent must already exist: {output_dir.parent}")
    groups, codes = load_jobs(args)
    outputs = {aa: output_dir / f"pose_diagnostics_{aa.lower()}.tsv.gz" for aa in codes}
    for output in outputs.values():
        if output.exists() and not args.overwrite:
            raise ValueError(f"Output exists; use --overwrite to replace it: {output}")
    for aa in codes:
        root = args.vina_d_dir if aa in D_AAS else args.vina_dir
        if not any((root / name).is_dir() for name in (aa.lower(), f"docking_{aa.lower()}")):
            raise ValueError(f"{aa}: no {aa.lower()}/ or docking_{aa.lower()}/ directory under {root}")
    if not args.receptor_dir.is_dir():
        raise ValueError(f"Receptor directory does not exist: {args.receptor_dir}")
    files = point_files(args.points_dir, groups)
    metadata = load_uniprot_metadata(args)
    model_ranges = load_model_ranges(args.model_ranges)
    initialize_worker(args, codes)  # Validate imports, templates and parameters before starting.
    versions = []
    for name in ("rdkit", "meeko", "prolif", "numpy", "scipy"):
        try:
            installed = version(name)
        except PackageNotFoundError:
            installed = "unknown"
        versions.append(f"{name}={installed}")
    print("Versions: " + ", ".join(versions), flush=True)
    print(f"Bundle: {args.bundle}\nOutput directory: {output_dir}\nVina pose: MODEL 1", flush=True)
    print(f"H-bonds: <= {args.hbond_distance} A, DHA >= {args.hbond_angle} deg; "
          f"salt bridges: <= {args.salt_distance} A; clashes: overlap > {args.clash_overlap} A; "
          f"pocket proximity: <= {args.point_cutoff} A; "
          f"artificial-boundary exclusion: <= {args.boundary_cutoff} A", flush=True)
    total = sum(map(len, groups.values()))
    tasks = ((protein, jobs, files.get(protein), metadata.get(alphafold_accession(protein)),
              model_ranges.get(protein, {})) for protein, jobs in groups.items())
    completed = 0
    progress = PoseProgress(total)
    counts, reported = Counter(), Counter()
    print(f"Diagnostics: 0.000% (0/{total:,}); {len(groups):,} proteins; {args.jobs} workers", flush=True)
    output_dir.mkdir(exist_ok=True)
    aa_counts = Counter()
    with tempfile.TemporaryDirectory(prefix=".pose-diagnostics-", dir=output_dir) as temporary:
        stages = {aa: Path(temporary) / path.name for aa, path in outputs.items()}
        with ExitStack() as handles:
            writers = {}
            for aa, stage in stages.items():
                raw = handles.enter_context(stage.open("wb"))
                compressed = handles.enter_context(
                    gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0))
                handle = handles.enter_context(io.TextIOWrapper(compressed, encoding="utf-8", newline=""))
                writers[aa] = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writers[aa].writerow(FIELDS)
            for rows in protein_results(tasks, args, codes, progress):
                for row in rows:
                    writers[row[1]].writerow(row)
                    aa_counts[row[1]] += 1
                    status = row[-2]
                    counts[status] += 1
                    if status != "success" and reported[status] < 3:
                        print(f"{status}: pocket_id={row[0]}, aa={row[1]}: {row[-1]}",
                              file=sys.stderr, flush=True)
                        reported[status] += 1
                completed += len(rows)
        print("Rows by status: " + ", ".join(f"{key}={value:,}" for key, value in counts.items()), flush=True)
        if completed != total:
            raise ValueError("Incomplete row count; no output was replaced")
        for aa, output in outputs.items():
            os.replace(stages[aa], output)
            print(f"Wrote {aa_counts[aa]:,} rows to {output} "
                  f"({output.stat().st_size / 1_000_000:.2f} MB)", flush=True)
    return 0 if counts["success"] == total else 2


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bundle", type=Path, default=Path("results/compact/L"))
    parser.add_argument("--points-dir", type=Path, default=Path("visualizations_p2rank"),
                        help="Recursively read original P2Rank *_points.pdb(.gz) files")
    parser.add_argument("--receptor-dir", type=Path, default=Path("receptors_pdbqt"))
    parser.add_argument("--structure-dir", type=Path, default=Path("structures_clean"),
                        help="Original clean AlphaFold PDB files with UniProt DBREF coverage")
    parser.add_argument("--uniprot-tsv", type=Path,
                        help="Entry/Length table; default: sole metadata/*_uniprot.tsv(.gz)")
    parser.add_argument("--model-ranges", type=Path,
                        help="Optional TSV: protein, chain, uniprot_start, uniprot_end")
    parser.add_argument("--vina-dir", type=Path, default=Path("docking_vina"))
    parser.add_argument("--vina-d-dir", type=Path, default=Path("docking_vina_d"))
    parser.add_argument("--ligand-dir", type=Path, default=Path("ligands/amino_acids_sdf"))
    parser.add_argument("--ligand-d-dir", type=Path, default=Path("ligands/d_amino_acids_sdf"))
    parser.add_argument("--output-dir", type=Path,
                        help="Default: <bundle>/pose_diagnostics; one gzip table per AA")
    parser.add_argument("--aa", nargs="+", type=str.upper, help="Optional subset, e.g. --aa MET ALA")
    parser.add_argument("--jobs", type=int, default=1, help="Protein workers (default: 1)")
    parser.add_argument("--point-cutoff", type=float, default=2.0)
    parser.add_argument("--clash-overlap", type=float, default=0.4)
    parser.add_argument("--hbond-distance", type=float, default=3.5)
    parser.add_argument("--hbond-angle", type=float, default=130.0)
    parser.add_argument("--salt-distance", type=float, default=4.5)
    parser.add_argument("--boundary-cutoff", type=float, default=6.0,
                        help="Mark poses near artificial ends/caps unreliable (angstroms; default 6)")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    for name in ("point_cutoff", "clash_overlap", "hbond_distance", "salt_distance", "boundary_cutoff"):
        value = getattr(args, name)
        if not math.isfinite(value) or value <= 0:
            parser.error(f"--{name.replace('_', '-')} must be finite and positive")
    if args.jobs < 1 or not math.isfinite(args.hbond_angle) or not 0 <= args.hbond_angle < 180:
        parser.error("--jobs must be positive; --hbond-angle must be in [0, 180)")
    if args.boundary_cutoff < max(args.hbond_distance, args.salt_distance):
        parser.error("--boundary-cutoff must cover the H-bond and salt-bridge distance cutoffs")
    try:
        return export(args)
    except (OSError, ValueError, RuntimeError, KeyError, csv.Error, EOFError) as exc:
        print(f"Export failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
