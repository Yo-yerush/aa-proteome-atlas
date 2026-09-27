#!/usr/bin/env python3
"""Export Vina MODEL 1 and saved SFCT-selected poses for a compact atlas bundle.

Writes aa_positions/positions_<aa>.tsv.gz with pocket_id, pose_id and positions.
pose_id is the Vina MODEL number; the SFCT selection is sfct_best_pose + 1.
Each pocket/pose pair is written once, including when SFCT selects MODEL 1.
Reads the existing score tables and manifest without modifying them.
Missing or invalid models are reported and leave positions empty; no substitute
pose is used. Missing or invalid SFCT indices are reported and omitted.
Supports both <root>/<aa>/ and <root>/docking_<aa>/ docking folders.
Run separately with --bundle results/compact/L or --bundle results/compact/D.
"""

import argparse
import csv
from decimal import Decimal, InvalidOperation
import gzip
import json
import math
import os
from pathlib import Path
import sys
import tempfile

from prepare_atlas_data import D_LIGANDS, FORMAT, LIGANDS, write_tsv


def with_progress(items, label, total):
    print(f"{label}: 0% (0/{total:,})", flush=True)
    last_percent = 0
    for count, item in enumerate(items, 1):
        yield item
        percent = count * 100 // total if total else 100
        if percent > last_percent:
            print(f"{label}: {percent}% ({count:,}/{total:,})", flush=True)
            last_percent = percent


def read_coordinates(pdbqt, pose_ids):
    try:
        handle = pdbqt.open()
    except FileNotFoundError:
        return None

    coordinates = {pose_id: "" for pose_id in pose_ids}
    errors = {pose_id: f"MODEL {pose_id} is missing" for pose_id in pose_ids}
    seen = set()
    current = None
    with handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith("MODEL"):
                try:
                    model = int(line.split()[1])
                    if model < 1:
                        raise ValueError
                except (IndexError, ValueError) as exc:
                    raise ValueError(f"line {line_number}: invalid MODEL number") from exc
                current = model if model in coordinates else None
                atoms = []
                invalid = ""
                if current is not None:
                    if current in seen:
                        coordinates[current] = ""
                        errors[current] = f"MODEL {current} occurs more than once"
                        current = None
                    else:
                        seen.add(current)
                        errors[current] = f"MODEL {current} is incomplete (missing ENDMDL)"
                continue
            if current is None:
                continue
            if line.startswith("ENDMDL"):
                if invalid:
                    errors[current] = invalid
                elif not atoms:
                    errors[current] = f"MODEL {current} contains no atoms"
                else:
                    coordinates[current] = ";".join(atoms)
                    errors.pop(current)
                current = None
            elif line.startswith(("ATOM", "HETATM")):
                try:
                    atom_name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    if not atom_name or not all(math.isfinite(value) for value in (x, y, z)):
                        raise ValueError
                except ValueError:
                    invalid = invalid or f"MODEL {current}, line {line_number}: invalid atom coordinates"
                    continue
                atoms.append(f"{atom_name}:{x:.3f},{y:.3f},{z:.3f}")
    return coordinates, errors


def sfct_model_id(value):
    if not value.strip():
        raise ValueError("missing sfct_best_pose")
    try:
        index = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError(f"invalid sfct_best_pose: {value!r}") from exc
    if not index.is_finite() or index < 0 or index != index.to_integral_value():
        raise ValueError(f"invalid sfct_best_pose: {value!r}; expected a nonnegative integer")
    return int(index) + 1


def bundle_file(bundle, filename):
    path = bundle / filename
    if path.resolve().parent != bundle:
        raise ValueError(f"Expected a file directly inside the bundle: {filename}")
    return path


def read_pockets(path):
    pockets = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        if not {"pocket_id", "protein", "pocket"}.issubset(reader.fieldnames or []):
            raise ValueError(f"{path}: missing pocket_id, protein or pocket column")
        for row in reader:
            key = row["pocket_id"]
            if not key or not row["protein"] or not row["pocket"] or key in pockets:
                raise ValueError(f"{path}, line {reader.line_num}: invalid or duplicate pocket ID")
            pockets[key] = (row["protein"], row["pocket"])
    return pockets


def export_positions(source, target, aa, docking_dir, pockets, total):
    # Older runs use docking_ala/; newer workflows use ala/.
    ligand_dirs = [path for path in (
        docking_dir / aa.lower(),
        docking_dir / f"docking_{aa.lower()}",
    ) if path.is_dir()]
    if not ligand_dirs:
        raise ValueError(
            f"{aa}: no {aa.lower()}/ or docking_{aa.lower()}/ folder in {docking_dir.resolve()}"
        )
    found = 0
    processed = 0
    successful = 0
    missing_files = 0
    missing_poses = 0
    invalid_selections = 0
    seen = set()
    with gzip.open(source, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        headers = reader.fieldnames or []
        if not {"pocket_id", "sfct_best_pose"}.issubset(headers) or len(headers) != len(set(headers)):
            raise ValueError(f"{source}: missing pocket_id/sfct_best_pose or duplicate columns")

        def rows():
            nonlocal found, processed, successful, missing_files, missing_poses, invalid_selections
            for row in with_progress(reader, aa, total):
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"{source}, line {reader.line_num}: invalid TSV row")
                pocket_id = row["pocket_id"]
                if pocket_id not in pockets:
                    raise ValueError(f"{source}: unknown pocket_id {pocket_id}")
                if pocket_id in seen:
                    raise ValueError(f"{source}, line {reader.line_num}: duplicate pocket_id {pocket_id}")
                seen.add(pocket_id)
                processed += 1
                pose_ids = [1]
                try:
                    selected = sfct_model_id(row["sfct_best_pose"])
                except ValueError as exc:
                    invalid_selections += 1
                    print(f"{aa}: pocket_id={pocket_id}: {exc}; SFCT-selected pose omitted",
                          file=sys.stderr, flush=True)
                else:
                    if selected != 1:
                        pose_ids.append(selected)
                protein, pocket = pockets[pocket_id]
                for ligand_dir in ligand_dirs:
                    pdbqt = ligand_dir / f"{protein}_{pocket}_{aa}.pdbqt"
                    try:
                        pose_data = read_coordinates(pdbqt, pose_ids)
                    except ValueError as exc:
                        pose_data = ({}, {pose_id: str(exc) for pose_id in pose_ids})
                    if pose_data is not None:
                        break
                if pose_data is None:
                    missing_files += 1
                    coordinates = {}
                    errors = {pose_id: "missing file" for pose_id in pose_ids}
                else:
                    coordinates, errors = pose_data
                successful += row.get("vina_status", "").strip().lower() == "success"
                for pose_id in pose_ids:
                    positions = coordinates.get(pose_id, "")
                    if pose_id in errors:
                        missing_poses += 1
                        print(f"{aa}: pocket_id={pocket_id}, pose_id={pose_id}: "
                              f"{errors[pose_id]} ({pdbqt.resolve()})",
                              file=sys.stderr, flush=True)
                    found += bool(positions)
                    yield (pocket_id, pose_id, positions)

        result = write_tsv(target, ("pocket_id", "pose_id", "positions"), rows())
    if processed != total:
        raise ValueError(f"{source}: row count differs from manifest.json")
    print(f"{aa}: {found:,} poses exported; {missing_files:,} missing files; "
          f"{missing_poses:,} missing or invalid poses; "
          f"{invalid_selections:,} missing or invalid SFCT selections", flush=True)
    if successful and not found:
        raise ValueError(
            f"{aa}: no coordinates recovered despite {successful:,} successful Vina rows. "
            "Check the reported paths and pose format above. No position tables were replaced."
        )
    return result


def export(bundle, vina_dir, vina_d_dir):
    bundle = bundle.resolve()
    positions_dir = bundle / "aa_positions"
    manifest_path = bundle / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("format") != FORMAT:
        raise ValueError("The selected folder is not a compact atlas bundle.")

    jobs = []
    seen = set()
    for entry in manifest["ligands"]:
        aa = entry["code"]
        if aa not in LIGANDS or aa in seen:
            raise ValueError(f"Unsupported or duplicate amino-acid code: {aa}")
        seen.add(aa)
        source = bundle_file(bundle, entry["file"])
        if source.name != f"scores_{aa.lower()}.tsv.gz":
            raise ValueError(f"Unexpected score filename for {aa}: {source.name}")
        docking_dir = vina_d_dir if aa in D_LIGANDS else vina_dir
        if not docking_dir.is_dir():
            raise ValueError(f"Docking directory does not exist: {docking_dir}")
        target = bundle_file(positions_dir, f"positions_{aa.lower()}.tsv.gz")
        jobs.append((entry, source, docking_dir, target))

    if not jobs:
        raise ValueError("The selected bundle has no amino-acid score tables.")
    print(f"Bundle: {bundle}", flush=True)
    pockets = read_pockets(bundle_file(bundle, manifest["pockets"]["file"]))

    # Finish all position tables before publishing them.
    with tempfile.TemporaryDirectory(prefix=".vina-positions-", dir=bundle) as temporary:
        stage = Path(temporary)
        for entry, source, docking_dir, target in jobs:
            export_positions(
                source, stage / target.name, entry["code"],
                docking_dir, pockets, entry["rows"],
            )
        positions_dir.mkdir(exist_ok=True)
        for _, _, _, target in jobs:
            os.replace(stage / target.name, target)

    print(f"Wrote {len(jobs)} position tables to {positions_dir}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=Path("results/compact/L"),
                        help="Existing compact bundle; position tables are saved in its aa_positions subfolder "
                             "(default: results/compact/L)")
    parser.add_argument("--vina-dir", type=Path, default=Path("docking_vina"),
                        help="L-AA Vina docking root")
    parser.add_argument("--vina-d-dir", type=Path, default=Path("docking_vina_d"),
                        help="D-AA Vina docking root")
    args = parser.parse_args()
    try:
        export(args.bundle, args.vina_dir, args.vina_d_dir)
    except (OSError, ValueError, KeyError, csv.Error, EOFError) as exc:
        print(f"Export failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
