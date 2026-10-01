"""Fragment identity regressions; only temporary synthetic pipeline data is written."""
import contextlib
import csv
import gzip
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "pipeline_scripts/scripts/prepare_atlas_data.py"
SPEC = importlib.util.spec_from_file_location("prepare_atlas_data", SCRIPT)
atlas = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(atlas)


class FragmentExportTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="atlas-fragment-test-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.source = self.root / "results"
        self.source.mkdir()
        self.output = self.root / "compact/L"
        self.rows = [dict(
            uniprot_id="A2VEC9", protein=f"AF-A2VEC9-F{fragment}-model_v6", pocket="pocket3",
            rank="3", probability="0.9", mean_pocket_plddt="95", center_x="1", center_y="2",
            center_z="3", residue_ids="A_1 A_2", vina_ala_affinity=f"-{fragment}.123456789012345",
            vina_status="success", sfct_status="missing" if fragment == 17 else "success",
            sfct_vina_score="" if fragment == 17 else f"-{fragment}.2",
            sfct_score="" if fragment == 17 else f"-{fragment}.4",
            vina_sfct_combined="" if fragment == 17 else f"-{fragment}.36",
            sfct_best_pose="" if fragment == 17 else "0", sfct_n_poses="" if fragment == 17 else "9",
        ) for fragment in [1, 17, 20, 5, 6]]

    def write(self, source, rows):
        path = self.source / f"{source}_ala_all_pockets.tsv"
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=self.rows[0], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        return path

    def export(self):
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            return atlas.export(self.source, self.source, self.output, configuration="L")

    def read_table(self, entry):
        data = gzip.decompress((self.output / entry["file"]).read_bytes())
        self.assertEqual(hashlib.sha256(data).hexdigest(), entry["sha256"])
        rows = list(csv.DictReader(io.StringIO(data.decode()), delimiter="\t"))
        self.assertEqual(len(rows), entry["rows"])
        return rows

    def test_l_only_export_preserves_fragments_missing_status_and_precision(self):
        paths = [self.write("sfct", self.rows), self.write("vina", list(reversed(self.rows)))]
        before = [path.read_bytes() for path in paths]
        manifest = self.export()  # There is deliberately no D directory.
        self.assertEqual(manifest["version"], 1)
        self.assertEqual([entry["code"] for entry in manifest["ligands"]], ["ALA"])
        pockets = {row["pocket_id"]: row for row in self.read_table(manifest["pockets"])}
        scores = self.read_table(manifest["ligands"][0])
        self.assertEqual(len(pockets), 5)
        self.assertEqual(len(scores), 5)
        by_model = {row["protein"]: row for row in self.rows}
        for score in scores:
            pocket = pockets[score["pocket_id"]]
            original = by_model[pocket["protein"]]
            self.assertEqual(pocket["uniprot_id"], "A2VEC9")
            self.assertEqual(pocket["pocket"], "pocket3")
            self.assertEqual(score["vina_affinity"], original["vina_ala_affinity"])
            for field in atlas.SFCT_FIELDS + ("sfct_status", "vina_status"):
                self.assertEqual(score[field], original[field])
        self.assertEqual(before, [path.read_bytes() for path in paths])
        self.assertEqual(json.loads((self.output / "manifest.json").read_text()), manifest)

    def test_true_duplicate_still_fails_in_either_source(self):
        for source in ["vina", "sfct"]:
            with self.subTest(source=source):
                self.write(source, self.rows + [self.rows[0]])
                with self.assertRaisesRegex(atlas.ExportError, "duplicate protein/pocket.*AF-A2VEC9-F1-model_v6"):
                    self.export()
                (self.source / f"{source}_ala_all_pockets.tsv").unlink()
                self.assertFalse(self.output.exists())

    def test_same_model_geometry_mismatch_still_fails(self):
        self.write("sfct", self.rows)
        self.write("vina", [dict(row, center_x="99") for row in self.rows])
        with self.assertRaisesRegex(atlas.ExportError, "metadata mismatch"):
            self.export()

    def test_different_models_never_acquire_each_others_scores(self):
        self.write("sfct", [self.rows[0]])
        self.write("vina", [self.rows[2]])
        manifest = self.export()
        rows = self.read_table(manifest["ligands"][0])
        self.assertEqual([(row["sfct_status"], row["vina_status"]) for row in rows],
                         [("success", "missing"), ("missing", "success")])

    def test_versions_stay_distinct_and_accession_can_be_inferred(self):
        rows = [dict(self.rows[0], uniprot_id="", protein=f"AF-A2VEC9-F1-model_v{v}") for v in [4, 6]]
        self.write("sfct", rows)
        manifest = self.export()
        pockets = self.read_table(manifest["pockets"])
        self.assertEqual(len(pockets), 2)
        self.assertEqual({row["uniprot_id"] for row in pockets}, {"A2VEC9"})


if __name__ == "__main__":
    unittest.main()
