#!/usr/bin/env python3
"""Sample receptor electrostatic potentials at retained pockets and Vina MODEL 1.

Run from a species project directory (Python >=3.10):
    python scripts/export_electrostatics.py \
        --bundle results/compact/L results/compact/D \
        --points-dir visualizations_p2rank --jobs 4

Requires numpy, scipy, pdb2pqr (with PROPKA), and the APBS executable.
Pose exports additionally require rdkit, meeko>=0.6 and AmberTools antechamber.
--pockets-only omits ligand preparation/poses and their dependencies.
No other scripts are imported; original structures, docking files and atlas
manifests are never modified. This script does not run docking or SFCT.

Outputs:
  <bundle-parent>/electrostatics/manifest.json
  <bundle-parent>/electrostatics/pocket_summary.tsv.gz
  <bundle-parent>/electrostatics/points_<protein>.tsv.gz
  <bundle>/pose_electrostatics/electrostatics_<aa>.tsv.gz

Select L and D together to share one receptor preparation/coarse APBS solution.
A fine solve per pocket covers all requested poses and pocket points. Grids,
PQRs and external-program scratch files are temporary. --keep-grids requires
an explicit --protein selection and retains local DX files only for those
candidates. --output-dir changes the shared output directory, not pose roots.

Chemistry and interpretation:
  PDB2PQR AMBER charges/radii, PROPKA at --ph; no heavy-atom debumping or
  hydrogen-network optimization. Every input heavy atom must be preserved;
  only a biological terminal OXT may be added. Unknown residues/parameters,
  internal chain breaks and inconsistent coverage fail explicitly.
  UniProt Entry/Length and original AlphaFold DBREF coverage identify termini.
  Artificial ends receive modeled ACE/NME caps BEFORE PDB2PQR/PROPKA. All
  prepared residues/caps are checked against the installed AMBER parameter
  table. Caps contribute to the receptor field. A pocket/pose within
  --boundary-cutoff of an artificial terminal residue or cap is unreliable.
  Other fragment results have status fragment_only: absent residues' long-range
  contributions are not restored by capping.

  Ligand chemistry and stereochemistry come from the original AA SDF. Meeko
  restores MODEL 1 heavy atoms and polar Hs; nonpolar Hs may be modeled. AM1-BCC
  charges are computed once at the original SDF geometry (no gas-phase geometry
  minimization), then mapped by the full chemical graph. Charge assignment
  never changes the docked coordinates or the ligand protonation state.
  qphi_kT = sum(q_i/e * phi_i/(kT/e)): a fixed-receptor-field descriptor,
  NOT binding free energy, desolvation energy, or a surface-complementarity
  correlation. All atoms, including hydrogens, enter qphi and pose statistics.
  Heavy-atom overlaps > --clash-overlap give unreliable_clash.

  APBS solves nonlinear Poisson-Boltzmann with AMBER radii, protein/solvent
  dielectric constants, and symmetric monovalent salt. Samples use trilinear
  interpolation of the unrounded potential grid (kT/e). Point visualization
  values alone are rounded to two decimals; fractions/statistics use the full
  precision samples. No zero-fill, extrapolation or partial-pocket averages.

Point contract for Mol*:
  point_index is ZERO-BASED in ATOM/HETATM encounter order within pocket rank,
  matching export_pocket_points.py. Coordinates are not duplicated.
  points_sha256 hashes ASCII "x,y,z;x,y,z;..." with each coordinate formatted
  to exactly 3 decimals, negative zero normalized to zero, no trailing newline.
  Check the fingerprint and n_points against the displayed point geometry.
  Join by protein + pocket_rank; pocket_id belongs to its own L/D bundle.
  Use a per-point color callback with a fixed red/white/blue kT/e scale.
  This script exports data only; it does not modify the web application.

Metadata:
  --uniprot-tsv: local Entry/Length TSV(.gz); optional Sequence/Fragment.
  Default: the sole metadata/*_uniprot.tsv(.gz). No network downloads.
  --model-ranges: optional TSV protein,chain,uniprot_start,uniprot_end;
  inclusive UniProt positions, '.' for blank chain. Without a DBREF/range,
  an exact unique Sequence match or a complete chain of full Length is needed.

Output publication is staged. Existing files need --overwrite; include every
previously exported bundle when updating shared data. Exit 0: all full-protein
results succeeded; 2: tables include failures or fragment_only results;
1: export aborted. Missing/failed poses have blank vina_pose and scores.
External tools have --timeout seconds per invocation. --max-grid-points and
--max-coarse-spacing bound per-worker grids; oversize domains fail rather than
silently reducing the requested resolution.

References:
https://apbs.readthedocs.io/en/latest/using/input/old/elec/bcfl.html
https://apbs.readthedocs.io/en/latest/using/input/old/elec/usemap.html
https://apbs.readthedocs.io/en/latest/formats/opendx.html
https://pdb2pqr.readthedocs.io/en/latest/using/index.html
https://ambermd.org/tutorials/basic/tutorial5/index.php
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import csv
import gzip
import hashlib
from importlib.metadata import PackageNotFoundError, version
from importlib.resources import files as package_files
import io
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

FORMAT = "aa-pocket-electrostatics"
L_AAS = "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
D_AAS = {"D" + aa for aa in L_AAS if aa != "GLY"}
RESIDUE_CODES = dict(zip(L_AAS, "ARNDCQEGHILKMFPSTWYV"))
POINT_FIELDS = ("pocket_rank", "point_index", "phi_kT_per_e")
SUMMARY_FIELDS = ("protein", "pocket_rank", "n_points", "points_sha256", "phi_mean",
                  "phi_min", "phi_max", "fraction_positive", "fraction_negative",
                  "status", "error")
POSE_FIELDS = ("pocket_id", "vina_pose", "n_atoms", "phi_mean", "phi_min", "phi_max",
               "qphi_kT", "status", "error")
GOOD = {"success", "fragment_only"}
AD_ELEMENTS = {"H": "H", "HD": "H", "HS": "H", "C": "C", "A": "C", "N": "N",
               "NA": "N", "NS": "N", "O": "O", "OA": "O", "OS": "O",
               "S": "S", "SA": "S", "P": "P", "F": "F", "Cl": "Cl", "Br": "Br",
               "I": "I", "Mg": "Mg", "Ca": "Ca", "Mn": "Mn", "Fe": "Fe", "Zn": "Zn"}

class DiagnosticError(ValueError):
    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


def dependencies(poses=True):
    global np, cKDTree, Chem, PDBQTMolecule, RDKitMolCreate
    try:
        import numpy as np
        from scipy.spatial import cKDTree
        if poses:
            from rdkit import Chem
            from meeko import PDBQTMolecule, RDKitMolCreate
    except ImportError as exc:
        raise ValueError("Missing dependency: use Python >=3.10 with numpy, scipy, "
                         "pdb2pqr; poses also need rdkit and meeko>=0.6. " + str(exc)) from exc


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def error_info(exc, fallback="processing_error"):
    return getattr(exc, "status", fallback), " ".join(str(exc).split())[:400]


@contextmanager
def gzip_writer(path, fields):
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                writer = csv.writer(text, delimiter="\t", lineterminator="\n")
                writer.writerow(fields)
                yield writer


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def run_program(command, cwd, timeout, name):
    """No shell; bounded logs are read only after completion, full logs stay in scratch."""
    log = cwd / f"{name}.log"
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with log.open("wb") as handle:
        try:
            completed = subprocess.run(command, cwd=cwd, stdout=handle,
                                       stderr=subprocess.STDOUT, env=env, timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            raise DiagnosticError(name + "_error", f"{name} exceeded {timeout:g} seconds") from exc
    with log.open("rb") as handle:
        handle.seek(max(0, log.stat().st_size - 12000))
        tail = handle.read().decode("utf-8", errors="replace")
    if completed.returncode:
        raise DiagnosticError(name + "_error", f"{name} exit {completed.returncode}: {tail[-1800:]}")
    return tail


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
    args.uniprot_tsv = source.resolve()
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
    mol = Chem.AddHs(supplier[0], addCoords=True)
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


def load_selection(args):
    keep = {(row["protein"], row["pocket"])
            for row in table_rows(args.retained, {"protein", "pocket"})}
    bundles, proteins, labels, seen_bundles = [], {}, {}, set()
    selected = set(args.protein or [])
    for bundle_index, bundle in enumerate(args.bundle):
        bundle = bundle.resolve()
        if bundle in seen_bundles:
            raise ValueError(f"Repeated bundle: {bundle}")
        seen_bundles.add(bundle)
        manifest_path = bundle / "manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest.get("format") != "aa-proteome-atlas-compact":
            raise ValueError(f"Not a compact atlas bundle: {bundle}")
        meta_path = bundle_file(bundle, manifest["pockets"]["file"])
        pockets, all_ids = {}, set()
        for row in table_rows(meta_path, {"pocket_id", "protein", "pocket", "rank"}):
            pid, protein, pocket = row["pocket_id"], filename_part(row["protein"]), filename_part(row["pocket"])
            rank = float(row["rank"])
            if pid in all_ids or not math.isfinite(rank) or rank < 1 or not rank.is_integer():
                raise ValueError(f"{meta_path}: duplicate pocket ID or invalid rank")
            all_ids.add(pid)
            if (protein, pocket) not in keep or (selected and protein not in selected):
                continue
            rank = int(rank)
            if (protein, pocket) in labels and labels[(protein, pocket)] != rank:
                raise ValueError(f"Conflicting pocket ranks across bundles: {protein} {pocket}")
            labels[(protein, pocket)] = rank
            groups = proteins.setdefault(protein, {})
            group = groups.setdefault(rank, {"pocket": pocket, "poses": []})
            if group["pocket"] != pocket:
                raise ValueError(f"Conflicting pocket names: {protein}, rank {rank}")
            pockets[pid] = (protein, rank)
        info = {"path": str(bundle), "manifest_sha256": fingerprint(manifest_path),
                "pockets_sha256": fingerprint(meta_path), "scores": {}, "codes": [], "rows": {}}
        for entry in manifest["ligands"]:
            aa = entry["code"]
            if aa not in set(L_AAS) | D_AAS or aa in info["codes"]:
                raise ValueError(f"{bundle}: unsupported or repeated AA {aa}")
            info["codes"].append(aa)
            info["rows"][aa] = 0
            if args.pockets_only:
                continue
            path = bundle_file(bundle, entry["file"])
            seen = set()
            for row in table_rows(path, {"pocket_id"}):
                pid = row["pocket_id"]
                if pid not in all_ids or pid in seen:
                    raise ValueError(f"{path}: unknown or repeated pocket_id {pid}")
                seen.add(pid)
                if pid in pockets:
                    protein, rank = pockets[pid]
                    proteins[protein][rank]["poses"].append((bundle_index, aa, pid))
                    info["rows"][aa] += 1
            if len(seen) != entry["rows"]:
                raise ValueError(f"{path}: row count differs from manifest")
            info["scores"][aa] = fingerprint(path)
        bundles.append(info)
    if selected - proteins.keys():
        raise ValueError(f"No retained pockets for requested proteins: {sorted(selected - proteins.keys())}")
    if not proteins:
        raise ValueError("No retained bundle pockets to export")
    parents = {Path(info["path"]).parent for info in bundles}
    if len(parents) != 1:
        raise ValueError("Selected L/D bundles must share the same compact parent directory")
    return bundles, proteins, next(iter(parents))


def point_files(root, proteins):
    if not root.is_dir():
        raise ValueError(f"Missing points directory: {root}")
    found = {}
    for path in root.rglob("*_points.pdb*"):
        suffix = next((s for s in ("_points.pdb.gz", "_points.pdb") if path.name.endswith(s)), None)
        if suffix is None or not path.is_file():
            continue
        stem = path.name[:-len(suffix)]
        protein = stem if stem in proteins else Path(stem).stem
        if protein in proteins:
            if protein in found:
                raise ValueError(f"Multiple point files for {protein}")
            found[protein] = path
    return found


def read_points(path, ranks):
    if path is None:
        raise DiagnosticError("missing_points", "No matching P2Rank point file")
    arrays = {rank: [] for rank in ranks}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            rank = int(line[22:26])
            if rank in arrays:
                xyz = tuple(float(line[start:start + 8]) for start in (30, 38, 46))
                if not all(math.isfinite(v) for v in xyz):
                    raise DiagnosticError("points_error", "Nonfinite P2Rank point")
                arrays[rank].append(xyz)
    return {rank: np.array(values, dtype=float).reshape((-1, 3)) for rank, values in arrays.items()}


def point_hash(xyz):
    def number(value):
        text = f"{value:.3f}"
        return "0.000" if text == "-0.000" else text
    text = ";".join(",".join(number(value) for value in point) for point in xyz)
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def load_amber_parameters():
    resource = package_files("pdb2pqr").joinpath("dat", "AMBER.DAT")
    contents = resource.read_text(encoding="utf-8")
    params = {}
    for line in contents.splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if len(fields) < 4:
            raise ValueError("Malformed AMBER parameter table")
        residue, atom, charge, radius = fields[:4]
        params.setdefault(residue, {})[atom] = (float(charge), float(radius))
    for residue in ("ACE", "NME", "ALA", "NALA", "CALA"):
        if residue not in params:
            raise ValueError(f"Missing {residue} in AMBER parameter table")
    return params, hashlib.sha256(contents.encode("utf-8")).hexdigest()


def mol2_charges(path, original):
    """
    Read AM1-BCC charges from Antechamber Mol2 and map them back to the
    original SDF atoms using element + unchanged 3D coordinates.

    Mol2/Sybyl resonance representations are NOT used to define ligand
    chemistry or stereochemistry. The original SDF remains authoritative.
    """
    from scipy.optimize import linear_sum_assignment

    text = path.read_text(encoding="utf-8")

    rows = []
    bonds = []
    section = None

    element_map = {
        "H": "H",
        "C": "C",
        "N": "N",
        "O": "O",
        "S": "S",
        "P": "P",
        "F": "F",
        "CL": "Cl",
        "BR": "Br",
        "I": "I",
    }

    for line in text.splitlines():

        if line.startswith("@<TRIPOS>"):
            section = line.strip()
            continue

        if not line.strip():
            continue

        if section == "@<TRIPOS>ATOM":
            fields = line.split()

            if len(fields) < 9:
                raise ValueError(
                    "Unsupported Mol2 atom row or missing partial charge"
                )

            atom_id = int(fields[0])

            if atom_id != len(rows) + 1:
                raise ValueError(
                    "Unsupported Mol2 atom numbering/order"
                )

            atom_type = fields[5].split(".")[0].upper()

            if atom_type not in element_map:
                raise ValueError(
                    f"Unsupported Mol2 atom type: {fields[5]}"
                )

            element = element_map[atom_type]

            xyz = np.array(
                [float(fields[2]), float(fields[3]), float(fields[4])],
                dtype=float
            )

            charge = float(fields[8])

            if not np.isfinite(xyz).all() or not math.isfinite(charge):
                raise ValueError(
                    "Nonfinite Mol2 coordinates or charge"
                )

            rows.append({
                "element": element,
                "xyz": xyz,
                "charge": charge
            })

        elif section == "@<TRIPOS>BOND":
            fields = line.split()

            if len(fields) < 4:
                raise ValueError("Malformed Mol2 bond row")

            a = int(fields[1]) - 1
            b = int(fields[2]) - 1

            bonds.append((a, b))

    if len(rows) != original.GetNumAtoms():
        raise ValueError(
            f"AM1-BCC Mol2 atom count differs from original SDF: "
            f"{len(rows)} vs {original.GetNumAtoms()}"
        )

    if original.GetNumConformers() != 1:
        raise ValueError(
            "Original ligand must contain exactly one conformer"
        )

    original_conf = original.GetConformer()
    original_xyz = np.asarray(
        original_conf.GetPositions(),
        dtype=float
    )

    original_elements = [
        atom.GetSymbol()
        for atom in original.GetAtoms()
    ]

    mol2_elements = [
        row["element"]
        for row in rows
    ]

    # Confirm that the elemental composition is identical.
    if Counter(original_elements) != Counter(mol2_elements):
        raise ValueError(
            "AM1-BCC Mol2 elemental composition differs from original SDF"
        )

    mol2_xyz = np.array(
        [row["xyz"] for row in rows],
        dtype=float
    )

    # Map Mol2 atoms to original atoms by element and 3D coordinates.
    # Antechamber is run with maxcyc=0, so coordinates should be unchanged
    # apart from tiny file-format rounding.
    mol2_to_original = {}
    maximum_displacement = 0.0

    for element in sorted(set(original_elements)):

        original_indices = [
            i for i, symbol in enumerate(original_elements)
            if symbol == element
        ]

        mol2_indices = [
            i for i, symbol in enumerate(mol2_elements)
            if symbol == element
        ]

        original_group = original_xyz[original_indices]
        mol2_group = mol2_xyz[mol2_indices]

        distance_matrix = np.linalg.norm(
            original_group[:, None, :] -
            mol2_group[None, :, :],
            axis=2
        )

        original_assignment, mol2_assignment = (
            linear_sum_assignment(distance_matrix)
        )

        for oi, mi in zip(
            original_assignment,
            mol2_assignment
        ):
            distance = float(distance_matrix[oi, mi])

            maximum_displacement = max(
                maximum_displacement,
                distance
            )

            original_index = original_indices[oi]
            mol2_index = mol2_indices[mi]

            mol2_to_original[mol2_index] = original_index

    # maxcyc=0 should preserve coordinates essentially exactly.
    # 0.02 A allows harmless Mol2 coordinate rounding.
    if maximum_displacement > 0.02:
        raise ValueError(
            "AM1-BCC output moved/reordered atoms beyond coordinate "
            f"matching tolerance: max displacement "
            f"{maximum_displacement:.4f} A"
        )

    if len(mol2_to_original) != original.GetNumAtoms():
        raise ValueError(
            "Could not uniquely map every AM1-BCC atom to original SDF"
        )

    # Validate atom CONNECTIVITY while deliberately ignoring Mol2 bond order.
    # Sybyl resonance/bond-order conventions differ from the original SDF.
    mol2_edges = set()

    for a, b in bonds:

        if (
            a not in mol2_to_original
            or b not in mol2_to_original
        ):
            raise ValueError(
                "Mol2 bond refers to an unmapped atom"
            )

        oa = mol2_to_original[a]
        ob = mol2_to_original[b]

        mol2_edges.add(
            tuple(sorted((oa, ob)))
        )

    original_edges = {
        tuple(sorted((
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx()
        )))
        for bond in original.GetBonds()
    }

    if mol2_edges != original_edges:
        missing = original_edges - mol2_edges
        added = mol2_edges - original_edges

        raise ValueError(
            "AM1-BCC output changed atom connectivity: "
            f"missing={sorted(missing)}, "
            f"added={sorted(added)}"
        )

    # Transfer Mol2 charges into ORIGINAL SDF atom order.
    charges = np.empty(
        original.GetNumAtoms(),
        dtype=float
    )

    for mol2_index, row in enumerate(rows):
        original_index = mol2_to_original[mol2_index]
        charges[original_index] = row["charge"]

    formal_charge = Chem.GetFormalCharge(original)

    if not np.isfinite(charges).all():
        raise ValueError(
            "Nonfinite AM1-BCC partial charges"
        )

    charge_sum = float(charges.sum())

    if abs(charge_sum - formal_charge) > 0.0025:
        raise ValueError(
            "AM1-BCC partial charges do not sum to the "
            f"original formal charge: "
            f"sum={charge_sum:.6f}, expected={formal_charge}"
        )

    return charges.tolist()


def prepare_ligands(args, codes, scratch):
    records = {}
    for number, aa in enumerate(codes, 1):
        source = (args.ligand_d_dir if aa in D_AAS else args.ligand_dir) / f"{aa}.sdf"
        original = load_ligand(source)
        if original.GetNumConformers() != 1 or not original.GetConformer().Is3D():
            raise ValueError(f"{source}: expected one original 3D conformer")
        if (not np.isfinite(original.GetConformer().GetPositions()).all()
                or any(atom.GetNumRadicalElectrons() for atom in original.GetAtoms())):
            raise ValueError(f"{source}: nonfinite coordinates or radical ligand")
        folder = scratch / aa
        folder.mkdir()
        Chem.MolToMolFile(original, str(folder / "original.mol"))
        print(f"AM1-BCC charges: {aa} ({number}/{len(codes)})", flush=True)
        # Sybyl output types let RDKit verify the charged molecule's full graph.
        run_program([args.antechamber, "-i", "original.mol", "-fi", "mdl",
                     "-o", "charged.mol2", "-fo", "mol2", "-c", "bcc",
                     "-nc", str(Chem.GetFormalCharge(original)), "-at", "sybyl", "-pf", "y",
                     "-ek", "qm_theory='AM1', scfconv=1.d-10, maxcyc=0, ndiis_attempts=700,"],
                    folder, args.timeout, "antechamber")
        charges = mol2_charges(folder / "charged.mol2", original)
        records[aa] = {"sdf": str(source), "sdf_sha256": fingerprint(source),
                       "molblock": Chem.MolToMolBlock(original), "charges": charges,
                       "formal_charge": Chem.GetFormalCharge(original),
                       "charge_model": "AM1-BCC (AmberTools antechamber; fixed original SDF conformer, maxcyc=0)"}
    return records


@dataclass
class PqrAtom:
    name: str
    residue: str
    chain: str
    number: int
    xyz: object
    charge: float = 0.0
    radius: float = 0.0
    cap: bool = False

    @property
    def element(self):
        return self.name.lstrip("0123456789")[0]


def read_structure(path):
    if not path.is_file():
        raise DiagnosticError("missing_structure", str(path))
    chains, models = {}, 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            kind = line[:6].strip()
            if kind == "MODEL":
                models += 1
                if models > 1:
                    raise ValueError("Multiple structure models are not supported")
            if kind not in ("ATOM", "HETATM"):
                continue
            if line[16:17].strip() or line[26:27].strip():
                raise ValueError("Alternate locations/insertion codes require an unambiguous structure")
            residue = line[17:20].strip()
            if residue not in L_AAS:
                raise ValueError(f"Unsupported receptor residue {residue}; no atoms are silently dropped")
            element = line[76:78].strip() or line[12:16].strip().lstrip("0123456789")[0]
            if element == "H":
                continue  # PDB2PQR assigns a documented protonation state from heavy atoms.
            if element not in ("C", "N", "O", "S"):
                raise ValueError(f"Unsupported receptor element: {element}")
            chain, number, name = line[21:22].strip(), int(line[22:26]), line[12:16].strip()
            xyz = np.array([float(line[start:start + 8]) for start in (30, 38, 46)])
            if not np.isfinite(xyz).all():
                raise ValueError("Nonfinite receptor coordinates")
            record = chains.setdefault(chain, {}).setdefault(number, {"residue": residue, "atoms": {}})
            if record["residue"] != residue or name in record["atoms"]:
                raise ValueError(f"Duplicate or inconsistent receptor residue: {chain}:{number}")
            record["atoms"][name] = xyz
    if not chains:
        raise ValueError("Structure has no protein atoms")
    for chain, residues in chains.items():
        numbers = sorted(residues)
        if len(numbers) < 2 or numbers != list(range(numbers[0], numbers[-1] + 1)):
            raise DiagnosticError("coverage_error", f"{chain}: incomplete/unsupported chain numbering")
        for number in numbers:
            if not {"N", "CA", "C", "O"}.issubset(residues[number]["atoms"]):
                raise ValueError(f"Missing backbone atoms: {chain}:{number}")
        for a, b in zip(numbers, numbers[1:]):
            distance = np.linalg.norm(residues[a]["atoms"]["C"] - residues[b]["atoms"]["N"])
            if not 0.9 <= distance <= 1.8:
                raise ValueError(f"Internal peptide break: {chain}:{a}-{b}, C-N={distance:.3f} A")
    return chains


def classify_coverage(protein, path, chains, metadata, explicit):
    accession = alphafold_accession(protein)
    if accession is None or metadata is None:
        raise DiagnosticError("coverage_error", f"{protein}: missing exact UniProt Entry/Length")
    length, sequence, fragment = metadata
    if fragment:
        raise DiagnosticError("coverage_error", f"{accession}: UniProt sequence is itself a fragment")
    references = pdb_reference_ranges(path, accession)
    coverage = {}
    for chain, residues in chains.items():
        numbers = sorted(residues)
        model_sequence = "".join(RESIDUE_CODES[residues[n]["residue"]] for n in numbers)
        bounds = explicit.get(chain)
        if chain in references:
            first, last, start, end = references[chain]
            if (first, last) != (numbers[0], numbers[-1]):
                raise DiagnosticError("coverage_error", f"{protein}:{chain}: DBREF/structure mismatch")
            if bounds is not None and bounds != (start, end):
                raise DiagnosticError("coverage_error", f"{protein}:{chain}: conflicting modeled ranges")
            bounds = (start, end)
        if bounds is None and sequence:
            start = sequence.find(model_sequence)
            if start >= 0 and sequence.find(model_sequence, start + 1) < 0:
                bounds = (start + 1, start + len(numbers))
        if bounds is None and len(numbers) == length:
            bounds = (1, length)
        if bounds is None:
            raise DiagnosticError("coverage_error", f"{protein}:{chain}: no verified modeled range; "
                                  "provide --model-ranges or original DBREF coverage")
        start, end = bounds
        if not 1 <= start <= end <= length or end - start + 1 != len(numbers):
            raise DiagnosticError("coverage_error", f"{protein}:{chain}: modeled range/Length mismatch")
        if sequence and sequence[start - 1:end] != model_sequence:
            raise DiagnosticError("coverage_error", f"{protein}:{chain}: sequence mismatch")
        coverage[chain] = {"uniprot_start": start, "uniprot_end": end, "full_length": length,
                           "biological_n": start == 1, "biological_c": end == length}
    return coverage


def methyl_hydrogens(center, bonded, plane):
    axis = unit_vector(bonded - center)
    reference = plane - center
    perpendicular = unit_vector(reference - np.dot(reference, axis) * axis)
    other = np.cross(axis, perpendicular)
    return [center + 1.09 * (-axis / 3 + math.sqrt(8 / 9) *
            (math.cos(angle) * perpendicular + math.sin(angle) * other))
            for angle in (0, 2 * math.pi / 3, 4 * math.pi / 3)]


def make_cap(residue, source, chain, number):
    if residue == "ACE":
        center, ca = source["N"], source["CA"]
        if "CD" in source and source.get("_proline", False):
            direction = remaining_bond_direction(center, [ca, source["CD"]])
        else:
            direction = planar_branches(center, ca, source["C"])[0]
        carbon = center + 1.33 * direction
        odir, mdir = planar_branches(carbon, center, ca)
        oxygen, methyl = carbon + 1.23 * odir, carbon + 1.50 * mdir
        coords = {"C": carbon, "O": oxygen, "CH3": methyl}
        hydrogens = methyl_hydrogens(methyl, carbon, oxygen)
    else:
        center, ca, oxygen = source["C"], source["CA"], source["O"]
        nitrogen = center + 1.33 * remaining_bond_direction(center, [ca, oxygen])
        mdir, _ = planar_branches(nitrogen, center, oxygen)
        methyl = nitrogen + 1.45 * mdir
        hydrogen = nitrogen + 1.01 * remaining_bond_direction(nitrogen, [center, methyl])
        coords = {"N": nitrogen, "H": hydrogen, "CH3": methyl}
        hydrogens = methyl_hydrogens(methyl, nitrogen, center)
    coords.update({f"HH3{i}": xyz for i, xyz in enumerate(hydrogens, 1)})
    return [PqrAtom(name, residue, chain, number, xyz, cap=True) for name, xyz in coords.items()]


def prepared_structure(chains, coverage, target):
    atoms, roles, boundary = [], {}, []
    occupied = set(chains) - {""}
    blank_chain = next((chr(i) for i in range(65, 91) if chr(i) not in occupied), None)
    for chain, residues in chains.items():
        output_chain = chain or blank_chain
        if output_chain is None or len(output_chain) != 1:
            raise ValueError("No supported PDB chain ID available")
        info, numbers = coverage[chain], sorted(residues)
        current = 1
        if not info["biological_n"]:
            source = dict(residues[numbers[0]]["atoms"])
            source["_proline"] = residues[numbers[0]]["residue"] == "PRO"
            cap = make_cap("ACE", source, output_chain, current)
            atoms.extend(cap)
            roles[(output_chain, current)] = "cap"
            boundary.extend(a.xyz for a in cap if a.element != "H")
            current += 1
        for number in numbers:
            record, coordinates = residues[number], dict(residues[number]["atoms"])
            artificial = ((number == numbers[0] and not info["biological_n"])
                          or (number == numbers[-1] and not info["biological_c"]))
            if number != numbers[-1] and "OXT" in coordinates:
                raise ValueError(f"Internal OXT at {chain}:{number}")
            if number == numbers[-1]:
                if not info["biological_c"]:
                    coordinates.pop("OXT", None)
                elif "OXT" not in coordinates:
                    coordinates["OXT"] = coordinates["C"] + 1.25 * remaining_bond_direction(
                        coordinates["C"], [coordinates["CA"], coordinates["O"]])
            if artificial:
                boundary.extend(coordinates.values())
            role = ("N" if number == numbers[0] and info["biological_n"] else
                    "C" if number == numbers[-1] and info["biological_c"] else "")
            roles[(output_chain, current)] = role
            atoms.extend(PqrAtom(name, record["residue"], output_chain, current, xyz)
                         for name, xyz in coordinates.items())
            current += 1
        if not info["biological_c"]:
            cap = make_cap("NME", residues[numbers[-1]]["atoms"], output_chain, current)
            atoms.extend(cap)
            roles[(output_chain, current)] = "cap"
            boundary.extend(a.xyz for a in cap if a.element != "H")
    if len(atoms) > 99999 or any(a.number > 9999 for a in atoms):
        raise ValueError("Prepared receptor exceeds PDB atom/residue numbering limits")
    with target.open("w", encoding="ascii", newline="\n") as handle:
        previous = None
        for serial, atom in enumerate(atoms, 1):
            if previous is not None and atom.chain != previous:
                handle.write("TER\n")
            name = atom.name if len(atom.name) == 4 else f" {atom.name:<3}"
            x, y, z = atom.xyz
            if any(len(f"{v:8.3f}") != 8 for v in atom.xyz):
                raise ValueError("Coordinate exceeds the PDB field width")
            handle.write(f"ATOM  {serial:5d} {name} {atom.residue:>3} {atom.chain}"
                         f"{atom.number:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
                         f"  1.00  0.00          {atom.element:>2}\n")
            previous = atom.chain
        handle.write("TER\nEND\n")
    return atoms, roles, np.array(boundary).reshape((-1, 3))


def read_pqr(path):
    atoms, names = [], set()
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if not fields or fields[0] not in ("ATOM", "HETATM"):
            continue
        if len(fields) != 11:
            raise ValueError("Expected whitespace PQR with retained chain IDs")
        _, serial, name, residue, chain, number, x, y, z, charge, radius = fields
        key = (chain, int(number), name)
        if key in names:
            raise ValueError("Duplicate PQR atom")
        names.add(key)
        atom = PqrAtom(name, residue, chain, int(number), np.array([float(x), float(y), float(z)]),
                       float(charge), float(radius))
        if (not np.isfinite([*atom.xyz, atom.charge, atom.radius]).all()
                or atom.radius < 0 or (atom.element != "H" and atom.radius == 0)):
            raise ValueError("Invalid PQR charge, radius or coordinate")
        atoms.append(atom)
    if not atoms:
        raise ValueError("PDB2PQR did not produce atoms")
    return atoms

AMBER_VARIANTS = {
    "HIS": ("HID", "HIE", "HIP"),
    "ASP": ("ASP", "ASH"),
    "GLU": ("GLU", "GLH"),
    "CYS": ("CYS", "CYM", "CYX"),
    "LYS": ("LYS", "LYN"),
}

def resolve_amber_template(prefix, residue, actual):
    """
    Resolve PDB2PQR residue names to the exact AMBER parameter template.

    PDB2PQR may retain a generic residue name such as HIS while the
    assigned AMBER parameters correspond to HID/HIE/HIP.  Select the
    template using the actual atom set and assigned charge/radius values.
    """
    variants = AMBER_VARIANTS.get(residue, (residue,))
    matches = []

    for variant in variants:
        template_name = prefix + variant
        template = AMBER.get(template_name)

        if template is None:
            continue

        if set(actual) != set(template):
            continue

        parameters_match = True

        for name, atom in actual.items():
            charge, radius = template[name]

            if (
                abs(atom.charge - charge) > 0.00011
                or abs(atom.radius - radius) > 0.00011
            ):
                parameters_match = False
                break

        if parameters_match:
            matches.append((template_name, template))

    if len(matches) == 1:
        return matches[0]

    candidates = ", ".join(prefix + v for v in variants)

    if not matches:
        raise ValueError(
            f"No matching AMBER template for {prefix + residue}; "
            f"tested [{candidates}]"
        )

    raise ValueError(
        f"Ambiguous AMBER template for {prefix + residue}: "
        + ", ".join(name for name, _ in matches)
    )

def validate_pqr(atoms, expected, roles):
    def heavy(items):
        return Counter((a.chain, a.number, a.element, *(round(float(v), 3) for v in a.xyz))
                       for a in items if a.element != "H")
    if heavy(atoms) != heavy(expected):
        raise ValueError("PDB2PQR moved, dropped or added unexpected heavy atoms")
    groups = {}
    for atom in atoms:
        key = (atom.chain, atom.number)
        if key not in roles:
            raise ValueError("PDB2PQR introduced an unexpected residue")
        atom.cap = roles[key] == "cap"
        groups.setdefault(key, {})[atom.name] = atom
    if groups.keys() != roles.keys():
        raise ValueError("PDB2PQR omitted a residue/cap")
    for key, actual in groups.items():
        residue = next(iter(actual.values())).residue
        prefix = "" if roles[key] == "cap" else roles[key]
        
        template_name, template = resolve_amber_template(
            prefix, residue, actual
        )
        if roles[key] == "cap" and abs(sum(a.charge for a in actual.values())) > 0.002:
            raise ValueError(f"Cap {key} is not electrically neutral")


def prepare_receptor(protein, metadata, model_ranges, scratch):
    source = OPTIONS.structure_dir / f"{protein}.pdb"
    chains = read_structure(source)
    coverage = classify_coverage(protein, source, chains, metadata, model_ranges)
    prepared, roles, boundary = prepared_structure(chains, coverage, scratch / "input.pdb")
    run_program([OPTIONS.pdb2pqr, "--ff=AMBER", "--ffout=AMBER", "--keep-chain", "--whitespace",
                 "--nodebump", "--noopt", "--titration-state-method=propka",
                 f"--with-ph={OPTIONS.ph}", "input.pdb", "receptor.pqr"],
                scratch, OPTIONS.timeout, "pdb2pqr")
    atoms = read_pqr(scratch / "receptor.pqr")
    validate_pqr(atoms, prepared, roles)
    is_fragment = any(not c["biological_n"] or not c["biological_c"] for c in coverage.values())
    provenance = {"structure_sha256": fingerprint(source),
                  "pqr_sha256": fingerprint(scratch / "receptor.pqr"),
                  "coverage": coverage, "fragment_only": is_fragment}
    return atoms, boundary, provenance


def grid_geometry(xyz, spacing, padding):
    lower, upper = xyz.min(axis=0) - padding, xyz.max(axis=0) + padding
    counts = np.maximum(33, np.ceil((upper - lower) / (32 * spacing)).astype(int) * 32 + 1)
    if int(counts.max()) > OPTIONS.max_grid_points:
        raise DiagnosticError("grid_size_error", f"Grid requires {counts.tolist()} points; "
                              f"limit is {OPTIONS.max_grid_points} per axis")
    lengths = (counts - 1) * spacing
    center = (lower + upper) / 2
    return counts, lengths, center


def solve_apbs(scratch, geometry, stem, coarse=False):
    counts, lengths, center = geometry
    vector = lambda values: " ".join(f"{float(v):.8g}" for v in values)
    read_map = "" if coarse else "    pot dx coarse.dx\n"
    boundary = "    bcfl sdh\n" if coarse else "    bcfl map\n    usemap pot 1\n"
    text = ("read\n    mol pqr receptor.pqr\n" + read_map + "end\n"
            "elec\n    mg-manual\n"
            f"    dime {' '.join(str(int(v)) for v in counts)}\n"
            f"    glen {vector(lengths)}\n    gcent {vector(center)}\n"
            "    mol 1\n    npbe\n" + boundary +
            f"    pdie {OPTIONS.protein_dielectric:g}\n"
            f"    sdie {OPTIONS.solvent_dielectric:g}\n"
            "    srfm smol\n    chgm spl2\n    sdens 10\n    srad 1.4\n    swin 0.3\n"
            f"    temp {OPTIONS.temperature:g}\n"
            f"    ion charge 1 conc {OPTIONS.salt:g} radius 2.0\n"
            f"    ion charge -1 conc {OPTIONS.salt:g} radius 2.0\n"
            "    calcenergy no\n    calcforce no\n"
            f"    write pot dx {stem}\nend\nquit\n")
    input_path = scratch / f"{stem}.in"
    input_path.write_text(text, encoding="ascii")
    log = run_program([OPTIONS.apbs, input_path.name], scratch, OPTIONS.timeout, "apbs")
    if re.search(r"(?:failed to converge|did not converge|not converged|FATAL|ERR\s*\(|"
                 r"exceeded.*iterations)", log, re.IGNORECASE):
        raise DiagnosticError("apbs_error", log[-1800:])
    outputs = [path for path in scratch.glob(f"{stem}*.dx") if path.is_file()]
    if len(outputs) != 1:
        raise DiagnosticError("apbs_error", "Expected one serial APBS potential DX output")
    target = scratch / f"{stem}.dx"
    if outputs[0] != target:
        outputs[0].replace(target)
    return target


def read_dx(path, geometry):
    """APBS OpenDX scalar order: z varies fastest, then y, then x."""
    counts, origin, deltas, items, initial = None, None, [], None, ""
    with path.open(encoding="ascii") as handle:
        for line in handle:
            fields = line.split()
            if "gridpositions" in fields and "counts" in fields:
                counts = tuple(int(v) for v in fields[fields.index("counts") + 1:])
            elif fields[:1] == ["origin"]:
                origin = np.array([float(v) for v in fields[1:]])
            elif fields[:1] == ["delta"]:
                deltas.append([float(v) for v in fields[1:]])
            elif "data follows" in line:
                header, initial = line.split("data follows", 1)
                match = re.search(r"rank\s+0\s+items\s+(\d+)", header)
                if match is None:
                    raise ValueError("Expected a scalar OpenDX array")
                items = int(match[1])
                break
        if (counts is None or len(counts) != 3 or min(counts) < 2 or origin is None
                or origin.shape != (3,) or np.shape(deltas) != (3, 3)
                or items != math.prod(counts)):
            raise ValueError("Incomplete OpenDX grid header")
        expected_counts, lengths, center = geometry
        delta = np.asarray(deltas)
        if (tuple(expected_counts) != counts
                or not np.allclose(delta, np.diag(lengths / (expected_counts - 1)), rtol=0, atol=1e-5)
                or not np.allclose(origin, center - lengths / 2, rtol=0, atol=1e-4)):
            raise ValueError("APBS output grid differs from requested geometry")
        values = np.empty(items, dtype=np.float64)
        offset, chunk, n_tokens = 0, [], 0
        if initial.strip():
            chunk.append(initial)
            n_tokens += len(initial.split())
        # Batches avoid a Python float object for every grid value.
        while offset < items:
            if n_tokens and (n_tokens >= 12288 or offset + n_tokens >= items):
                block = np.fromstring(" ".join(chunk), sep=" ")
                if len(block) != n_tokens or offset + n_tokens > items:
                    raise ValueError("Invalid OpenDX potential values")
                values[offset:offset + n_tokens] = block
                offset += n_tokens
                chunk, n_tokens = [], 0
                continue
            line = handle.readline()
            if not line:
                raise ValueError("Truncated OpenDX potential array")
            if line.lstrip().startswith("#") or not line.strip():
                continue
            if re.search(r"[A-Za-z]", line.replace("e", "").replace("E", "")):
                raise ValueError("Nonfinite or incomplete OpenDX potential array")
            chunk.append(line)
            n_tokens += len(line.split())
    if not np.isfinite(values).all():
        raise ValueError("Nonfinite APBS potential")
    return values.reshape(counts), origin, np.diag(delta)


def sample_grid(grid, xyz):
    values, origin, step = grid
    index = (xyz - origin) / step
    maximum = np.array(values.shape) - 1
    if not np.isfinite(index).all() or ((index < -1e-7) | (index > maximum + 1e-7)).any():
        raise DiagnosticError("outside_grid", "A requested sample lies outside the APBS grid")
    index = np.minimum(np.maximum(index, 0), maximum)
    lower = np.minimum(np.floor(index).astype(int), maximum - 1)
    fraction = index - lower
    result = np.zeros(len(xyz), dtype=float)
    for dx in (0, 1):
        for dy in (0, 1):
            for dz in (0, 1):
                shift = np.array([dx, dy, dz])
                weight = np.prod(np.where(shift, fraction, 1 - fraction), axis=1)
                i, j, k = (lower + shift).T
                result += weight * values[i, j, k]
    if not np.isfinite(result).all():
        raise ValueError("Nonfinite interpolated potential")
    return result


def vina_path(protein, pocket, aa):
    root = OPTIONS.vina_d_dir if aa in D_AAS else OPTIONS.vina_dir
    name = f"{protein}_{pocket}_{aa}.pdbqt"
    candidates = [root / folder / name for folder in (aa.lower(), "docking_" + aa.lower())]
    found = [path for path in candidates if path.is_file()]
    if len(found) > 1:
        raise ValueError(f"Ambiguous Vina output: {found}")
    return found[0] if found else candidates[0]


def too_close(xyz, boundary):
    return bool(len(boundary) and (cKDTree(boundary).query(xyz)[0] <= OPTIONS.boundary_cutoff).any())


def has_clash(ligand, receptor):
    """Heavy-atom overlap; cap proximity is handled independently."""
    tree, radii = receptor
    table = Chem.GetPeriodicTable()
    conf = ligand.GetConformer()
    for atom in ligand.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        xyz = np.array(conf.GetAtomPosition(atom.GetIdx()))
        radius = table.GetRvdw(atom.GetAtomicNum())
        for index in tree.query_ball_point(xyz, radius + float(radii.max())):
            if radius + radii[index] - np.linalg.norm(xyz - tree.data[index]) > OPTIONS.clash_overlap:
                return True
    return False


def failed_pose(pid, status, message):
    return (pid, "", "", "", "", "", "", status, message)


def summary_row(protein, rank, xyz, potential, status, message=""):
    base = (protein, rank, len(xyz), point_hash(xyz) if len(xyz) else "")
    if potential is None:
        return (*base, "", "", "", "", "", status, message)
    return (*base, f"{potential.mean():.4f}", f"{potential.min():.4f}", f"{potential.max():.4f}",
            f"{np.mean(potential > 0):.6f}", f"{np.mean(potential < 0):.6f}", status, message)


def init_worker(args, amber, ligand_records):
    global OPTIONS, AMBER, LIGANDS
    OPTIONS, AMBER = args, amber
    dependencies(not args.pockets_only)
    LIGANDS = {}
    for aa, record in ligand_records.items():
        mol = Chem.MolFromMolBlock(record["molblock"], removeHs=False, sanitize=True)
        if mol is None:
            raise ValueError(f"Cannot restore prepared {aa} chemistry")
        LIGANDS[aa] = (mol, np.array(record["charges"]))


def process_protein(task):
    protein, groups, points_path, metadata, model_ranges, stage = task
    filename = f"points_{protein}.tsv.gz"
    summaries, poses, provenance, grid_files = [], [], {}, []
    empty = np.empty((0, 3))
    points, point_error = {}, None
    try:
        points = read_points(points_path, groups)
        provenance["points_source_sha256"] = fingerprint(points_path)
    except (ValueError, OSError) as exc:
        point_error = error_info(exc, "points_error")
    with gzip_writer(stage / filename, POINT_FIELDS) as point_writer:
        def fail_pocket(rank, error):
            xyz = points.get(rank, empty)
            summaries.append(summary_row(protein, rank, xyz, None, *error))
            point_writer.writerows((rank, index, "") for index in range(len(xyz)))
            poses.extend((bi, aa, failed_pose(pid, *error)) for bi, aa, pid in groups[rank]["poses"])

        if point_error:
            for rank in sorted(groups):
                fail_pocket(rank, point_error)
            return filename, summaries, poses, provenance, grid_files

        with tempfile.TemporaryDirectory(prefix="aa_electro_", dir=OPTIONS.work_dir) as work:
            scratch = Path(work)
            try:
                atoms, boundary, receptor_info = prepare_receptor(protein, metadata, model_ranges, scratch)
                provenance.update(receptor_info)
            except (ValueError, OSError) as exc:
                for rank in sorted(groups):
                    fail_pocket(rank, error_info(exc, "receptor_error"))
                return filename, summaries, poses, provenance, grid_files

            status = "fragment_only" if provenance["fragment_only"] else "success"
            note = "Field excludes unmodeled sequence; fragment ends were capped" if status == "fragment_only" else ""
            receptor = None
            if not OPTIONS.pockets_only:
                heavy = [atom for atom in atoms if atom.element != "H" and not atom.cap]
                periodic = Chem.GetPeriodicTable()
                receptor = (cKDTree([a.xyz for a in heavy]),
                            np.array([periodic.GetRvdw(periodic.GetAtomicNumber(a.element)) for a in heavy]))
            prepared, geometries = {}, {}
            for rank in sorted(groups):
                xyz = points.get(rank, empty)
                if not len(xyz):
                    fail_pocket(rank, ("missing_points", f"No P2Rank points for rank {rank}"))
                    continue
                if too_close(xyz, boundary):
                    fail_pocket(rank, ("unreliable_fragment_boundary", "Pocket points are near an artificial chain end"))
                    continue
                entries, coordinates = [], [xyz]
                for bi, aa, pid in groups[rank]["poses"]:
                    try:
                        original, charges = LIGANDS[aa]
                        ligand = reconstruct_ligand(vina_path(protein, groups[rank]["pocket"], aa), original)
                        pose_xyz = ligand.GetConformer().GetPositions()
                        if too_close(pose_xyz, boundary):
                            raise DiagnosticError("unreliable_fragment_boundary", "Ligand is near an artificial chain end")
                        if has_clash(ligand, receptor):
                            raise DiagnosticError("unreliable_clash", "Heavy-atom overlap exceeds configured cutoff")
                        entries.append((bi, aa, pid, pose_xyz, charges))
                        coordinates.append(pose_xyz)
                    except (ValueError, OSError, RuntimeError, KeyError) as exc:
                        poses.append((bi, aa, failed_pose(pid, *error_info(exc, "ligand_error"))))
                try:
                    geometries[rank] = grid_geometry(np.concatenate(coordinates), OPTIONS.grid_spacing, OPTIONS.padding)
                    prepared[rank] = entries
                except ValueError as exc:
                    error = error_info(exc, "grid_error")
                    summaries.append(summary_row(protein, rank, xyz, None, *error))
                    point_writer.writerows((rank, index, "") for index in range(len(xyz)))
                    poses.extend((bi, aa, failed_pose(pid, *error)) for bi, aa, pid, *_ in entries)

            if geometries:
                try:
                    bounds = [np.array([a.xyz for a in atoms])]
                    for counts, lengths, center in geometries.values():
                        bounds.append(np.array([center - lengths / 2, center + lengths / 2]))
                    coarse_geometry = grid_geometry(np.concatenate(bounds), OPTIONS.max_coarse_spacing,
                                                    OPTIONS.coarse_padding)
                    coarse_file = solve_apbs(scratch, coarse_geometry, "coarse", coarse=True)
                    # Validate before APBS consumes the external boundary map.
                    coarse_grid = read_dx(coarse_file, coarse_geometry)
                    del coarse_grid
                    provenance["coarse_grid_counts"] = [int(n) for n in coarse_geometry[0]]
                except (ValueError, OSError) as exc:
                    error = error_info(exc, "apbs_error")
                    for rank, entries in prepared.items():
                        xyz = points[rank]
                        summaries.append(summary_row(protein, rank, xyz, None, *error))
                        point_writer.writerows((rank, index, "") for index in range(len(xyz)))
                        poses.extend((bi, aa, failed_pose(pid, *error)) for bi, aa, pid, *_ in entries)
                    return filename, summaries, poses, provenance, grid_files

            for rank, entries in prepared.items():
                xyz, grid_file = points[rank], None
                try:
                    grid_file = solve_apbs(scratch, geometries[rank], f"pocket_{rank}")
                    grid = read_dx(grid_file, geometries[rank])
                    potential = sample_grid(grid, xyz)
                    # Stage every pose first so a failure cannot leave duplicate rows.
                    sampled = []
                    for bi, aa, pid, pose_xyz, charges in entries:
                        try:
                            phi = sample_grid(grid, pose_xyz)
                            score = float(np.dot(charges, phi))
                            if not math.isfinite(score):
                                raise ValueError("Nonfinite qphi")
                            row = (pid, 1, len(phi), f"{phi.mean():.4f}", f"{phi.min():.4f}",
                                   f"{phi.max():.4f}", f"{score:.4f}", status, note)
                        except ValueError as exc:
                            row = failed_pose(pid, *error_info(exc, "sampling_error"))
                        sampled.append((bi, aa, row))
                    del grid
                    if OPTIONS.keep_grids:
                        grid_name = f"grid_{protein}_{rank}.dx.gz"
                        with grid_file.open("rb") as source, gzip.open(stage / grid_name, "wb") as target:
                            shutil.copyfileobj(source, target)
                        grid_files.append(grid_name)
                    summaries.append(summary_row(protein, rank, xyz, potential, status, note))
                    point_writer.writerows((rank, index, f"{value:.2f}") for index, value in enumerate(potential))
                    poses.extend(sampled)
                except (ValueError, OSError) as exc:
                    error = error_info(exc, "apbs_error")
                    summaries.append(summary_row(protein, rank, xyz, None, *error))
                    point_writer.writerows((rank, index, "") for index in range(len(xyz)))
                    poses.extend((bi, aa, failed_pose(pid, *error)) for bi, aa, pid, *_ in entries)
                finally:
                    if grid_file is not None:
                        grid_file.unlink(missing_ok=True)
    return filename, summaries, poses, provenance, grid_files


def completed_proteins(args, tasks, amber, ligands):
    if args.jobs == 1:
        init_worker(args, amber, ligands)
        for task in tasks:
            yield task[0], process_protein(task)
        return
    # Keep task submissions bounded rather than queueing a whole proteome.
    with ProcessPoolExecutor(max_workers=args.jobs, initializer=init_worker,
                             initargs=(args, amber, ligands)) as pool:
        pending = deque()
        for task in tasks:
            pending.append((task[0], pool.submit(process_protein, task)))
            if len(pending) >= args.jobs * 2:
                protein, future = pending.popleft()
                yield protein, future.result()
        while pending:
            protein, future = pending.popleft()
            yield protein, future.result()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bundle", type=Path, nargs="+", required=True,
                        help="Compact atlas bundle(s); select L and D together to share receptor solves")
    parser.add_argument("--points-dir", type=Path, required=True, help="Original P2Rank *_points.pdb[.gz] directory")
    parser.add_argument("--retained", type=Path, default=Path("results/pockets_for_docking.tsv"))
    parser.add_argument("--structure-dir", type=Path, default=Path("structures_clean"))
    parser.add_argument("--uniprot-tsv", type=Path, help="UniProt Entry/Length TSV(.gz); auto-select only if unambiguous")
    parser.add_argument("--model-ranges", type=Path, help="Explicit modeled ranges if original DBREF coverage is unavailable")
    parser.add_argument("--vina-dir", type=Path, default=Path("docking_vina"))
    parser.add_argument("--vina-d-dir", type=Path, default=Path("docking_vina_d"))
    parser.add_argument("--ligand-dir", type=Path, default=Path("ligands/amino_acids_sdf"))
    parser.add_argument("--ligand-d-dir", type=Path, default=Path("ligands/d_amino_acids_sdf"))
    parser.add_argument("--output-dir", type=Path, help="Shared output; default <bundle-parent>/electrostatics")
    parser.add_argument("--pockets-only", action="store_true", help="Skip ligand charges, poses and pose tables")
    parser.add_argument("--protein", action="append", help="Process just this structure stem; repeat for multiple candidates")
    parser.add_argument("--keep-grids", action="store_true", help="Retain focused DX.gz files; requires --protein")
    parser.add_argument("--overwrite", action="store_true", help="Replace this exporter's existing output files")
    parser.add_argument("--work-dir", type=Path, help="Existing scratch parent; default system temporary directory")
    parser.add_argument("--jobs", type=int, default=1, help="Concurrent proteins; APBS can require substantial RAM per worker")
    parser.add_argument("--timeout", type=float, default=1800, help="Timeout in seconds for each external-tool invocation")
    parser.add_argument("--pdb2pqr", default="pdb2pqr", help="PDB2PQR executable (e.g. pdb2pqr30)")
    parser.add_argument("--apbs", default="apbs")
    parser.add_argument("--antechamber", default="antechamber")
    parser.add_argument("--ph", type=float, default=7.0)
    parser.add_argument("--temperature", type=float, default=298.15, help="Kelvin")
    parser.add_argument("--salt", type=float, default=0.15, help="Molar concentration of each monovalent ion species")
    parser.add_argument("--protein-dielectric", type=float, default=4.0)
    parser.add_argument("--solvent-dielectric", type=float, default=78.5)
    parser.add_argument("--grid-spacing", type=float, default=0.5, help="Focused APBS grid spacing, Angstrom")
    parser.add_argument("--padding", type=float, default=8.0, help="Fine-grid margin around all pocket/pose samples, Angstrom")
    parser.add_argument("--coarse-padding", type=float, default=20.0, help="Coarse-grid margin around receptor and fine grids, Angstrom")
    parser.add_argument("--max-coarse-spacing", type=float, default=2.0, help="Coarse grid spacing, Angstrom; never exceeded")
    parser.add_argument("--max-grid-points", type=int, default=257, help="Maximum points per axis (32*n+1)")
    parser.add_argument("--boundary-cutoff", type=float, default=6.0, help="Unreliable proximity to artificial ends/caps, Angstrom")
    parser.add_argument("--clash-overlap", type=float, default=0.4, help="Maximum allowed heavy-atom vdW overlap, Angstrom")
    args = parser.parse_args()
    positive = ("temperature", "protein_dielectric", "solvent_dielectric", "grid_spacing", "padding",
                "coarse_padding", "max_coarse_spacing", "boundary_cutoff", "clash_overlap", "timeout")
    for key in positive:
        value = getattr(args, key)
        if not math.isfinite(value) or value <= 0:
            parser.error(f"--{key.replace('_', '-')} must be positive and finite")
    if not math.isfinite(args.ph) or not 0 <= args.ph <= 14:
        parser.error("--ph must be between 0 and 14")
    if not math.isfinite(args.salt) or args.salt < 0:
        parser.error("--salt must be nonnegative and finite")
    if args.jobs < 1 or args.max_grid_points < 33 or (args.max_grid_points - 1) % 32:
        parser.error("--jobs must be positive; --max-grid-points must be 32*n+1, n>=1")
    if args.max_coarse_spacing < args.grid_spacing:
        parser.error("--max-coarse-spacing must be at least --grid-spacing")
    if args.keep_grids and not args.protein:
        parser.error("--keep-grids requires an explicit --protein selection")
    for key, value in vars(args).items():
        if isinstance(value, Path):
            setattr(args, key, value.resolve())
    args.bundle = [path.resolve() for path in args.bundle]
    if args.work_dir is not None and not args.work_dir.is_dir():
        parser.error("--work-dir must be an existing directory")
    return args


def check_outputs(args, bundles, proteins, shared):
    if shared in args.bundle or shared == Path(args.retained).parent:
        raise ValueError("Use a dedicated shared electrostatics output directory")
    pose_targets = {}
    if not args.pockets_only:
        for bi, info in enumerate(bundles):
            for aa in info["codes"]:
                pose_targets[(bi, aa)] = Path(info["path"]) / "pose_electrostatics" / f"electrostatics_{aa.lower()}.tsv.gz"
    targets = [shared / "manifest.json", shared / "pocket_summary.tsv.gz", *pose_targets.values()]
    targets.extend(shared / f"points_{protein}.tsv.gz" for protein in proteins)
    if args.keep_grids:
        targets.extend(shared / f"grid_{protein}_{rank}.dx.gz" for protein, groups in proteins.items() for rank in groups)
    for target in targets:
        if target.exists() and (not args.overwrite or not target.is_file()):
            raise ValueError(f"Output exists (use --overwrite for files): {target}")
    previous = shared / "manifest.json"
    if previous.is_file():
        old = json.loads(previous.read_text(encoding="utf-8"))
        if old.get("format") != FORMAT:
            raise ValueError(f"Refusing to replace another format's manifest: {previous}")
        previous_bundles = {item["path"] for item in old["bundles"]}
        if not previous_bundles.issubset(info["path"] for info in bundles):
            raise ValueError("Include all previously exported L/D bundles when replacing shared electrostatics")
        if not set(old["proteins"]).issubset(proteins):
            raise ValueError("This selection would discard previous proteins; use another --output-dir")
        if not old["pockets_only"] and args.pockets_only:
            raise ValueError("Existing shared data include poses; include the same bundles/poses or use another --output-dir")
    # A candidate-only run must never replace full per-AA pose tables.
    if args.protein and any(target.exists() for target in pose_targets.values()):
        if not previous.is_file() or set(json.loads(previous.read_text(encoding="utf-8"))["proteins"]) != set(proteins):
            raise ValueError("Candidate selection cannot overwrite broader pose tables; use --pockets-only for candidate grids")
    return pose_targets


def software_info(args):
    packages = {}
    names = ["numpy", "scipy", "pdb2pqr", "propka"]
    if not args.pockets_only:
        names += ["rdkit", "meeko"]
    for name in names:
        try:
            packages[name] = version(name)
        except PackageNotFoundError:
            packages[name] = "unknown"
    executables = {}
    for name in (["pdb2pqr", "apbs"] if args.pockets_only else ["pdb2pqr", "apbs", "antechamber"]):
        requested = getattr(args, name)
        resolved = shutil.which(requested)
        if resolved is None:
            raise ValueError(f"Missing {name} executable: {requested}")
        executable = Path(resolved).resolve()
        setattr(args, name, str(executable))
        executables[name] = {"path": str(executable), "sha256": fingerprint(executable)}
    return {"packages": packages, "executables": executables, "python": sys.version.split()[0]}


def make_manifest(args, bundles, software, amber_hash, ligands, proteins, counts, grids):
    settings = {key: getattr(args, key) for key in (
        "ph", "temperature", "salt", "protein_dielectric", "solvent_dielectric", "grid_spacing",
        "padding", "coarse_padding", "max_coarse_spacing", "max_grid_points", "boundary_cutoff", "clash_overlap")}
    return {
        "format": FORMAT, "version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
        "pockets_only": args.pockets_only, "bundles": bundles, "proteins": proteins,
        "settings": settings, "software": software, "amber_parameters_sha256": amber_hash,
        "ligand_charge_models": ligands, "status_counts": dict(counts),
        "retained_pockets_sha256": fingerprint(args.retained),
        "model_ranges_sha256": fingerprint(args.model_ranges) if args.model_ranges else None,
        "uniprot_metadata_sha256": fingerprint(args.uniprot_tsv),
        "summary": {"file": "pocket_summary.tsv.gz", "columns": SUMMARY_FIELDS},
        "points": {"columns": POINT_FIELDS, "file_template": "points_{protein}.tsv.gz",
                   "join": ["protein", "pocket_rank"], "point_index_base": 0,
                   "order": "ATOM/HETATM encounter order within the original P2Rank pocket rank",
                   "coordinate_fingerprint": "SHA256 of ASCII x,y,z;x,y,z;...; fixed 3 decimals, -0.000 normalized to 0.000; no newline",
                   "potential_units": "kT/e", "decimals": 2,
                   "color_scale": {"minimum": -5, "midpoint": 0, "maximum": 5,
                                   "colors": ["#d73027", "#ffffff", "#4575b4"]}},
        "poses": {"columns": POSE_FIELDS, "file_template": "<bundle>/pose_electrostatics/electrostatics_<aa>.tsv.gz",
                  "pose": "Vina MODEL 1", "atom_selection": "all ligand atoms including H",
                  "charge_model": "AM1-BCC at fixed original AA SDF geometry, maxcyc=0, full chemical-graph mapping",
                  "qphi_definition": "sum_i (q_i/e) * (phi_receptor(r_i)/(kT/e))",
                  "qphi_units": "kT", "interpretation": "Negative is favorable in the fixed receptor field; not binding free energy or a desolvation-corrected score"},
        "solver": {"equation": "nonlinear Poisson-Boltzmann", "coarse_bc": "sdh", "focused_bc": "map",
                   "charge_discretization": "spl2", "surface": "smol", "probe_radius_A": 1.4,
                   "surface_window_A": 0.3, "surface_density": 10, "ion_radius_A": 2.0,
                   "sampling": "trilinear; no extrapolation; all required samples must be finite"},
        "chemistry": {"receptor": "PDB2PQR AMBER + PROPKA; no debumping/no H-network optimization",
                      "termini": "Biological ends handled by PDB2PQR at requested pH; artificial ends ACE/NME capped before PDB2PQR",
                      "caps_contribute_to_field": True, "fragment_status": "fragment_only even beyond the boundary cutoff",
                      "steric_radii": "RDKit periodic-table vdW radii; heavy atoms only"},
        "retained_grids": grids,
        "lazy_loading": "Load the summary once, the selected protein's point file on structure inspection, and only the selected bundle/AA pose table. Verify point count and coordinate fingerprint before applying colors. Blank values are unavailable, never zero.",
    }


def main():
    args = parse_args()
    dependencies(not args.pockets_only)
    software = software_info(args)
    amber, amber_hash = load_amber_parameters()
    metadata = load_uniprot_metadata(args)
    ranges = load_model_ranges(args.model_ranges)
    bundles, proteins, parent = load_selection(args)
    shared = args.output_dir or parent / "electrostatics"
    pose_targets = check_outputs(args, bundles, proteins, shared)
    points = point_files(args.points_dir, proteins)
    total_pockets = sum(len(groups) for groups in proteins.values())
    total_poses = sum(sum(info["rows"].values()) for info in bundles)
    print(f"Electrostatics: {len(proteins):,} proteins; {total_pockets:,} pockets; {total_poses:,} poses; {args.jobs} worker(s)", flush=True)
    print(f"Shared output: {shared}", flush=True)
    for info in bundles:
        if not args.pockets_only:
            print(f"Pose output: {Path(info['path']) / 'pose_electrostatics'}", flush=True)
    codes = sorted({aa for _, aa in pose_targets})
    with tempfile.TemporaryDirectory(prefix="aa_charges_", dir=args.work_dir) as work:
        ligands = prepare_ligands(args, codes, Path(work)) if codes else {}
    # Stage under the compact parent so publication normally stays on one filesystem.
    with tempfile.TemporaryDirectory(prefix=".electrostatics_", dir=parent) as work:
        stage = Path(work)
        counts, seen, provenance, grid_files = Counter(), Counter(), {}, []
        summary_path = stage / "pocket_summary.tsv.gz"
        staged_poses = {key: stage / f"pose_{key[0]}_{key[1].lower()}.tsv.gz" for key in pose_targets}
        done_pockets, done_poses, last_percent = 0, 0, -1
        last_update = time.monotonic()
        with ExitStack() as stack:
            summary_writer = stack.enter_context(gzip_writer(summary_path, SUMMARY_FIELDS))
            writers = {key: stack.enter_context(gzip_writer(path, POSE_FIELDS)) for key, path in staged_poses.items()}
            tasks = ((protein, groups, points.get(protein), metadata.get(alphafold_accession(protein)),
                      ranges.get(protein, {}), stage) for protein, groups in sorted(proteins.items()))
            print(f"Electrostatics: 0% (0/{len(proteins):,} proteins)", flush=True)
            for done, (protein, result) in enumerate(completed_proteins(args, tasks, amber, ligands), 1):
                filename, summaries, poses, info, grids = result
                if len(summaries) != len(proteins[protein]):
                    raise ValueError(f"Incomplete pocket export for {protein}")
                expected = {(bi, aa, pid) for group in proteins[protein].values() for bi, aa, pid in group["poses"]}
                actual = [(bi, aa, row[0]) for bi, aa, row in poses]
                if len(actual) != len(expected) or set(actual) != expected:
                    raise ValueError(f"Duplicate or incomplete pose export for {protein}")
                summary_writer.writerows(sorted(summaries, key=lambda row: row[1]))
                for bi, aa, row in sorted(poses, key=lambda item: (item[0], item[1], item[2][0])):
                    writers[(bi, aa)].writerow(row)
                    seen[(bi, aa)] += 1
                for row in [*summaries, *(row for _, _, row in poses)]:
                    status, error = row[-2:]
                    counts[status] += 1
                    if status not in GOOD and counts[status] <= 3:
                        print(f"{status}: {protein}: {error}", file=sys.stderr, flush=True)
                provenance[protein] = {**info, "points_file": filename,
                                       "points_file_sha256": fingerprint(stage / filename)}
                grid_files.extend(grids)
                done_pockets += len(summaries)
                done_poses += len(poses)
                percent = done * 100 // len(proteins)
                if percent != last_percent or time.monotonic() - last_update >= 30:
                    print(f"Electrostatics: {percent}% ({done:,}/{len(proteins):,} proteins; "
                          f"{done_pockets:,} pockets; {done_poses:,} poses)", flush=True)
                    last_percent = percent
                    last_update = time.monotonic()
        for (bi, aa), path in staged_poses.items():
            if seen[(bi, aa)] != bundles[bi]["rows"][aa]:
                raise ValueError(f"Incomplete pose table for bundle {bi}, {aa}")
            bundles[bi].setdefault("pose_files", {})[aa] = {"file": f"pose_electrostatics/electrostatics_{aa.lower()}.tsv.gz",
                                                          "rows": seen[(bi, aa)], "sha256": fingerprint(path)}
        manifest = make_manifest(args, bundles, software, amber_hash, ligands, provenance, counts, grid_files)
        manifest["summary"]["sha256"] = fingerprint(summary_path)
        write_json(stage / "manifest.json", manifest)
        shared.mkdir(parents=True, exist_ok=True)
        for key, source in staged_poses.items():
            target = pose_targets[key]
            target.parent.mkdir(parents=True, exist_ok=True)
            os.replace(source, target)
        # Copy to a temporary sibling for --output-dir on a different filesystem.
        publish = ["pocket_summary.tsv.gz", *(v["points_file"] for v in provenance.values()), *grid_files, "manifest.json"]
        for name in publish:
            with tempfile.NamedTemporaryFile(prefix=".electrostatics_", dir=shared, delete=False) as handle:
                temporary = Path(handle.name)
            try:
                shutil.copyfile(stage / name, temporary)
                os.replace(temporary, shared / name)
            finally:
                temporary.unlink(missing_ok=True)
    print("Status counts (pocket + pose rows): " + ", ".join(f"{key}={value:,}" for key, value in sorted(counts.items())), flush=True)
    print(f"Saved electrostatics manifest: {shared / 'manifest.json'}", flush=True)
    return 0 if set(counts) <= {"success"} else 2


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Export interrupted; check the last published manifest before using output.", file=sys.stderr)
        sys.exit(130)
    except (ValueError, OSError, RuntimeError, KeyError, ImportError) as exc:
        print(f"Export failed: {exc}", file=sys.stderr)
        sys.exit(1)
