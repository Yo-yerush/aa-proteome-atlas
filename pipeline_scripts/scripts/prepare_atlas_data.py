#!/usr/bin/env python3
"""Export compact AA atlas data using only the Python standard library.

Keep full pipeline results as the archive; publish the generated compact folder.
No scores are rounded, combined, ranked, or quality-filtered by this exporter.
Read collected *_all_pockets.tsv(.gz) tables from results/ for L-AAs
and results/d_amino_acids/ for D-AAs. Legacy DMET tables in results/ are ignored.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile


FORMAT = "aa-proteome-atlas-compact"
VERSION = 1
L_LIGANDS = "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split()
D_LIGANDS = [f"D{aa}" for aa in L_LIGANDS if aa != "GLY"]
LIGANDS = L_LIGANDS + D_LIGANDS
METADATA = (
    "uniprot_id", "tair_id", "protein", "pocket", "rank", "score", "probability",
    "center_x", "center_y", "center_z", "sas_points", "surf_atoms",
    "n_pocket_residues", "n_plddt_matched", "mean_pocket_plddt", "min_pocket_plddt",
    "fraction_plddt_ge70", "fraction_plddt_ge90", "residue_ids",
    "gene_id",
)
SCORES = (
    "vina_affinity", "vina_status", "sfct_vina_score", "sfct_score",
    "vina_sfct_combined", "sfct_best_pose", "sfct_n_poses", "sfct_status",
)
NUMERIC_MATCH = ("rank", "probability", "mean_pocket_plddt", "center_x", "center_y", "center_z")
SFCT_FIELDS = ("sfct_vina_score", "sfct_score", "vina_sfct_combined", "sfct_best_pose", "sfct_n_poses")


class ExportError(ValueError):
    pass


def find_input(directory: Path, source: str, aa: str) -> Path | None:
    stem = f"{source}_{aa.lower()}_all_pockets.tsv"
    matches = [path for path in (directory / stem, directory / (stem + ".gz")) if path.is_file()]
    if len(matches) > 1:
        raise ExportError(f"Both compressed and plain inputs exist for {stem}; keep only the intended input in this directory.")
    return matches[0] if matches else None


def source_rows(path: Path, source: str, aa: str):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t", strict=True)
        headers = reader.fieldnames or []
        if len(set(headers)) != len(headers):
            raise ExportError(f"{path}: duplicate column names")
        required = {"protein", "pocket", "probability", "mean_pocket_plddt"}
        required.update({f"vina_{aa.lower()}_affinity"} if source == "vina" else SFCT_FIELDS[:3])
        missing = required.difference(headers)
        if "status" not in headers and f"{source}_status" not in headers:
            missing.add(f"status or {source}_status")
        if missing:
            raise ExportError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ExportError(f"{path}, line {reader.line_num}: unexpected number of TSV fields")
            protein = row["protein"]
            if not row.get("uniprot_id", "").strip():
                # Supports the original pipeline collectors before adding annotation columns.
                match = re.fullmatch(r"AF-(.+)-F\d+-model_\S+", protein)
                if not match:
                    raise ExportError(f"{path}, line {reader.line_num}: no UniProt ID and unrecognized AlphaFold model name")
                row["uniprot_id"] = match[1]
            if not all(row.get(key, "").strip() for key in ("uniprot_id", "protein", "pocket")):
                raise ExportError(f"{path}, line {reader.line_num}: empty protein or pocket identifier")
            if not row.get("gene_id", "").strip():
                row["gene_id"] = row.get("tair_id", "")  # Support older annotated tables.
            yield row


def number(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return math.nan


def same_number(left: str, right: str) -> bool:
    a, b = number(left), number(right)
    return a == b or (math.isnan(a) and math.isnan(b)) or abs(a - b) <= 1e-9


def check_same_pocket(previous: dict, row: dict, aa: str, path: Path) -> None:
    consistent = previous["protein"] == row["protein"]
    consistent &= all(same_number(previous.get(key, ""), row.get(key, "")) for key in NUMERIC_MATCH)
    consistent &= sorted(previous.get("residue_ids", "").split()) == sorted(row.get("residue_ids", "").split())
    if not consistent:
        raise ExportError(f"{aa} / {row['uniprot_id']} / {row['pocket']}: Vina/SFCT model or pocket metadata mismatch ({path.name})")


def merge_ligand(aa: str, files: dict[str, Path | None]) -> dict:
    # Preserve the legacy loader's order: SFCT rows first, then Vina-only pockets.
    merged = {}
    for source in ("sfct", "vina"):
        path = files[source]
        if path is None:
            continue
        seen = set()
        for row in source_rows(path, source, aa):
            key = (row["uniprot_id"], row["pocket"])
            if key in seen:
                raise ExportError(f"{path}: duplicate protein/pocket: {' / '.join(key)}")
            seen.add(key)
            if key in merged:
                check_same_pocket(merged[key]["metadata"], row, aa, path)
            else:
                merged[key] = {
                    "metadata": {field: row.get(field, "") for field in METADATA},
                    "scores": {field: "missing" if field.endswith("_status") else "" for field in SCORES},
                }
            scores = merged[key]["scores"]
            scores[f"{source}_status"] = row.get("status", row.get(f"{source}_status", "")).strip().lower()
            if source == "vina":
                scores["vina_affinity"] = row[f"vina_{aa.lower()}_affinity"]
            else:
                for field in SFCT_FIELDS:
                    scores[field] = row.get(field, "")
    return merged


def write_tsv(path: Path, headers, rows) -> dict:
    digest = hashlib.sha256()
    count = 0
    # A fixed timestamp/name makes gzip output reproducible across identical runs.
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", fileobj=raw, mode="wb", compresslevel=9, mtime=0) as compressed:
            def emit(values):
                cells = [str(value) for value in values]
                if any(any(char in cell for char in "\t\r\n\x00") for cell in cells):
                    raise ExportError(f"{path.name}: embedded tab/newline/NUL is not supported in atlas fields")
                data = ("\t".join(cells) + "\n").encode("utf-8")
                digest.update(data)
                compressed.write(data)

            emit(headers)
            for values in rows:
                emit(values)
                count += 1
    return {"file": path.name, "rows": count, "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def validate_output(output: Path, inputs: list[Path], overwrite: bool) -> None:
    # Allow a dedicated output subfolder, but never replace an input directory.
    for directory in inputs:
        if output == directory or output in directory.parents:
            raise ExportError("Output cannot be an input directory or contain an input directory.")
    if output.exists():
        if not output.is_dir():
            raise ExportError(f"Output is not a directory: {output}")
        if not any(output.iterdir()):
            return
        if not overwrite:
            raise ExportError(f"Output already contains files: {output}. Use --overwrite to update a previously generated bundle.")
        try:
            previous = json.loads((output / "manifest.json").read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            raise ExportError("Refusing to overwrite a nonempty folder without a valid atlas manifest.") from exc
        if previous.get("format") != FORMAT or previous.get("generator") != Path(__file__).name:
            raise ExportError("Refusing to overwrite a folder not created by this exporter.")


def export(
    vina_dir: Path,
    sfct_dir: Path,
    output: Path,
    overwrite: bool = False,
    *,
    vina_d_dir: Path = Path("results/d_amino_acids"),
    sfct_d_dir: Path = Path("results/d_amino_acids"),
    configuration: str = "both",
) -> dict:
    vina_dir, sfct_dir, output = vina_dir.resolve(), sfct_dir.resolve(), output.resolve()
    vina_d_dir, sfct_d_dir = vina_d_dir.resolve(), sfct_d_dir.resolve()
    configurations = {
        "L": (L_LIGANDS, [vina_dir, sfct_dir]),
        "D": (D_LIGANDS, [vina_d_dir, sfct_d_dir]),
        "both": (LIGANDS, [vina_dir, sfct_dir, vina_d_dir, sfct_d_dir]),
    }
    if configuration not in configurations:
        raise ExportError("Configuration must be L, D, or both.")
    selected_ligands, input_dirs = configurations[configuration]
    for directory in input_dirs:
        if not directory.is_dir():
            raise ExportError(f"Input directory does not exist: {directory}")
    validate_output(output, input_dirs, overwrite)
    files = {}
    for aa in selected_ligands:
        # D-AAs, including DMET, are read only from the D-AA roots.
        aa_vina_dir = vina_d_dir if aa in D_LIGANDS else vina_dir
        aa_sfct_dir = sfct_d_dir if aa in D_LIGANDS else sfct_dir
        files[aa] = {
            "sfct": find_input(aa_sfct_dir, "sfct", aa),
            "vina": find_input(aa_vina_dir, "vina", aa),
        }
    files = {aa: pair for aa, pair in files.items() if any(pair.values())}
    if not files:
        raise ExportError("No supported *_all_pockets.tsv or .tsv.gz inputs found. Do not use best-per-protein or filtered top-hit tables.")
    source_bytes = sum(path.stat().st_size for pair in files.values() for path in pair.values() if path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".atlas-export-", dir=output.parent) as temporary:
        stage = Path(temporary)
        pocket_ids, pocket_metadata, ligands = {}, [], []
        total_rows = 0
        for aa, pair in files.items():
            for source, path in pair.items():
                if path is None:
                    print(f"Warning: {aa} has no {source} file; those scores will be missing.", file=sys.stderr)
            merged = merge_ligand(aa, pair)

            def score_rows():
                for record in merged.values():
                    metadata = tuple(record["metadata"][field] for field in METADATA)
                    # Deduplicate exact metadata snapshots, not rounded values. Any per-AA
                    # representation differences get their own record and remain lossless.
                    if metadata not in pocket_ids:
                        pocket_ids[metadata] = len(pocket_metadata) + 1
                        pocket_metadata.append(metadata)
                    yield (pocket_ids[metadata], *(record["scores"][field] for field in SCORES))

            entry = write_tsv(stage / f"scores_{aa.lower()}.tsv.gz", ("pocket_id", *SCORES), score_rows())
            entry["code"] = aa
            entry["sources"] = {source: path.name if path else None for source, path in pair.items()}
            ligands.append(entry)
            total_rows += entry["rows"]
            print(f"{aa:4s} {entry['rows']:>8,d} pocket rows -> {entry['bytes'] / 1_000_000:.3f} MB", flush=True)
        if not pocket_metadata:
            raise ExportError("Input tables contain no pocket rows.")
        pockets = write_tsv(stage / "pockets.tsv.gz", ("pocket_id", *METADATA),
                            ((index, *values) for index, values in enumerate(pocket_metadata, 1)))
        manifest = {
            "format": FORMAT, "version": VERSION, "generator": Path(__file__).name,
            "pockets": pockets, "ligands": ligands,
            "input_bytes": source_bytes, "protein_aa_pocket_rows": total_rows,
            "score_precision": "Source score strings retained without rounding or recalculation",
            "quality_filters_applied": False, "failed_rows_retained": True,
            "sha256_basis": "Uncompressed UTF-8 TSV bytes",
        }
        (stage / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        output.mkdir(exist_ok=True)
        # Validation and generation finish before replacing any published data. Publish
        # the manifest last; its hashes let clients reject a partially updated bundle.
        for path in stage.iterdir():
            if path.name != "manifest.json":
                os.replace(path, output / path.name)
        os.replace(stage / "manifest.json", output / "manifest.json")
    generated_bytes = pockets["bytes"] + sum(entry["bytes"] for entry in ligands) + (output / "manifest.json").stat().st_size
    reduction = 100 * (1 - generated_bytes / source_bytes) if source_bytes else 0
    print(f"\nOriginal inputs: {source_bytes / 1_000_000:.2f} MB")
    print(f"Compact bundle:  {generated_bytes / 1_000_000:.2f} MB ({reduction:.1f}% smaller)")
    print(f"Shared metadata: {len(pocket_metadata):,} records; {total_rows:,} AA/pocket rows preserved")
    print(f"Output: {output}")
    print("Original input files were not modified. Publish only the compact bundle, not the full input folders.")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--vina-dir", type=Path, default=Path("results"), help="L-AA directory containing Vina all-pocket TSV(.gz) files")
    parser.add_argument("--sfct-dir", type=Path, default=Path("results"), help="L-AA directory containing SFCT all-pocket TSV(.gz) files")
    parser.add_argument("--vina-d-dir", type=Path, default=Path("results/d_amino_acids"), help="D-AA directory containing Vina all-pocket TSV(.gz) files")
    parser.add_argument("--sfct-d-dir", type=Path, default=Path("results/d_amino_acids"), help="D-AA directory containing SFCT all-pocket TSV(.gz) files")
    parser.add_argument("--configuration", choices=("L", "D", "both"), default="both", help="Which AA configuration to export (default: both)")
    parser.add_argument("--output", type=Path, default=Path("results/compact"), help="Separate generated bundle directory")
    parser.add_argument("--overwrite", action="store_true", help="Replace previously generated bundle files; never delete original inputs")
    args = parser.parse_args()
    try:
        export(
            args.vina_dir, args.sfct_dir, args.output, args.overwrite,
            vina_d_dir=args.vina_d_dir, sfct_d_dir=args.sfct_d_dir,
            configuration=args.configuration,
        )
    except (ExportError, OSError, csv.Error, EOFError) as exc:
        print(f"Export failed: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
