#!/usr/bin/env python3
"""Export retained P2Rank pocket points as pocket_id/points in one TSV.

Uses the selected compact bundle's existing pocket IDs and keeps only pockets
listed in the retained-pocket table. Reads P2Rank *_points.pdb(.gz) files
recursively under --points-dir. These files must match the original predictions.
Writes pocket_points.tsv in the bundle's parent directory (default: results/compact).
"""

import argparse
import csv
import gzip
import json
import math
import os
from pathlib import Path
import sys
import tempfile


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
            if any(not row[field].strip() for field in required):
                raise ValueError(f"{path}, line {reader.line_num}: empty required field")
            yield row


def retained_pockets(bundle, retained):
    manifest = json.loads((bundle / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("format") != "aa-proteome-atlas-compact":
        raise ValueError("The selected folder is not a compact atlas bundle.")
    metadata = (bundle / manifest["pockets"]["file"]).resolve()
    if metadata.parent != bundle:
        raise ValueError("The pocket metadata file must be inside the bundle.")

    keep = {(row["protein"], row["pocket"])
            for row in table_rows(retained, {"protein", "pocket"})}
    groups = {}
    seen = set()
    for row in table_rows(metadata, {"pocket_id", "protein", "pocket", "rank"}):
        pocket_id = row["pocket_id"]
        if pocket_id in seen:
            raise ValueError(f"{metadata}: duplicate pocket_id {pocket_id}")
        seen.add(pocket_id)
        if (row["protein"], row["pocket"]) not in keep:
            continue
        rank = float(row["rank"])
        if not math.isfinite(rank) or rank < 1 or not rank.is_integer():
            raise ValueError(f"{metadata}: invalid pocket rank for pocket_id {pocket_id}")
        groups.setdefault(row["protein"], []).append(
            (pocket_id, row["pocket"], int(rank))
        )
    if not groups:
        raise ValueError("No retained pockets match the selected compact bundle.")
    return groups


def find_point_files(points_dir, proteins):
    if not points_dir.is_dir():
        raise ValueError(f"P2Rank point directory does not exist: {points_dir}")
    files = {}
    for path in points_dir.rglob("*_points.pdb*"):
        suffix = next((suffix for suffix in ("_points.pdb.gz", "_points.pdb")
                       if path.name.endswith(suffix)), None)
        if suffix is None or not path.is_file():
            continue
        structure = path.name[:-len(suffix)]
        # P2Rank includes the input structure extension, e.g. protein.pdb_points.pdb.gz.
        protein = structure if structure in proteins else Path(structure).stem
        if protein not in proteins:
            continue
        if protein in files:
            raise ValueError(
                f"Multiple point files for {protein}: {files[protein]} and {path}. "
                "Select a points directory containing only the intended prediction run."
            )
        files[protein] = path

    missing = [protein for protein in proteins if protein not in files]
    if missing:
        example = missing[0]
        raise ValueError(
            f"Missing P2Rank point files for {len(missing):,} retained proteins under "
            f"{points_dir}. Example: {example}.pdb_points.pdb.gz. "
            "Point files are normally in visualizations/data/ and are not generated "
            "with -visualizations 0. No output file was replaced."
        )
    return files


def read_points(path, ranks):
    points = {rank: [] for rank in ranks}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.startswith(("ATOM", "HETATM")):
                continue
            # P2Rank stores the pocket rank in PDB residue-number columns 23-26.
            try:
                rank = int(line[22:26])
            except ValueError as exc:
                raise ValueError(f"{path}, line {line_number}: invalid point pocket rank") from exc
            if rank not in points:  # Excludes rank 0 and all unretained pockets.
                continue
            values = [line[30:38].strip(), line[38:46].strip(), line[46:54].strip()]
            try:
                valid = all(math.isfinite(float(value)) for value in values)
            except ValueError:
                valid = False
            if not valid:
                raise ValueError(f"{path}, line {line_number}: invalid point coordinates")
            # Preserve the coordinate precision; discard only PDB padding.
            points[rank].append(",".join(values))
    return {rank: ";".join(values) for rank, values in points.items()}


def export(bundle, points_dir, retained):
    bundle = bundle.resolve()
    points_dir = points_dir.resolve()
    groups = retained_pockets(bundle, retained)
    total = sum(len(pockets) for pockets in groups.values())
    print(f"Retained pockets: {total:,} across {len(groups):,} proteins", flush=True)
    print(f"Finding P2Rank point files under {points_dir}", flush=True)
    files = find_point_files(points_dir, groups)
    output = bundle.parent / "pocket_points.tsv"

    completed = 0
    point_count = 0
    last_percent = 0
    print(f"Exporting points: 0% (0/{total:,})", flush=True)
    # Publish only after every retained pocket has been matched successfully.
    with tempfile.TemporaryDirectory(prefix=".pocket-points-", dir=bundle) as temporary:
        stage = Path(temporary) / output.name
        with stage.open("wb") as out:
            out.write(b"pocket_id\tpoints\n")
            for protein, pockets in groups.items():
                points = read_points(files[protein], {rank for _, _, rank in pockets})
                for pocket_id, pocket, rank in pockets:
                    coordinates = points[rank]
                    if not coordinates:
                        raise ValueError(
                            f"No points for pocket_id {pocket_id} ({protein}, {pocket}, "
                            f"rank {rank}) in {files[protein]}. "
                            "Check that the point files match the retained predictions. "
                            "No output file was replaced."
                        )
                    out.write(f"{pocket_id}\t{coordinates}\n".encode("utf-8"))
                    point_count += coordinates.count(";") + 1
                    completed += 1
                    percent = completed * 100 // total
                    if percent > last_percent:
                        print(f"Exporting points: {percent}% ({completed:,}/{total:,})",
                              flush=True)
                        last_percent = percent
        os.replace(stage, output)
    print(f"Wrote {completed:,} pocket rows and {point_count:,} points to {output}", flush=True)
    print(f"File size: {output.stat().st_size / 1_000_000:.2f} MB", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=Path("results/compact/L"),
                        help="Compact bundle supplying pocket IDs (default: results/compact/L)")
    parser.add_argument("--retained", type=Path, default=Path("results/pockets_for_docking.tsv"),
                        help="Retained-pocket TSV(.gz) (default: results/pockets_for_docking.tsv)")
    parser.add_argument("--points-dir", type=Path, default=Path("pockets_p2rank"),
                        help="Directory searched recursively for *_points.pdb(.gz) files")
    args = parser.parse_args()
    try:
        export(args.bundle, args.points_dir, args.retained)
    except (OSError, ValueError, KeyError, csv.Error, EOFError) as exc:
        print(f"Export failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
