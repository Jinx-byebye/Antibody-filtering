import argparse
import json
import tempfile
import unittest
from pathlib import Path

from scripts.prepare import build, optional_number


class PreparationTests(unittest.TestCase):
    def test_path_map_validation(self):
        from scripts.prepare import chain_sequences
        with self.assertRaises(FileNotFoundError):
            chain_sequences({"protenix_json_path": "/missing/input.json"}, [])

    def test_missing_and_nonfinite_measurements(self):
        self.assertIsNone(optional_number(""))
        self.assertIsNone(optional_number("nan"))
        self.assertEqual(optional_number("2.5"), 2.5)

    def test_candidate_source_join(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prediction = root / "prediction.json"
            prediction.write_text(json.dumps([{"sequences": [
                {"id": "H", "proteinChain": {"sequence": "AA"}},
                {"id": "A", "proteinChain": {"sequence": "A"}},
            ]}]))
            paths = {name: root / name for name in (
                "records.tsv", "panel.jsonl", "yeast.jsonl", "igdesign.jsonl",
                "source.jsonl", "affinity.csv", "cdr.csv",
            )}
            paths["records.tsv"].write_text(
                "benchmark\tpair_id\ttarget\tlabel\tantigen\theavy\tlight\n"
                "Yeast-display enrichment\tpair\ttarget\tpositive\tA\tAA\t\n")
            paths["panel.jsonl"].write_text(json.dumps({"pair_id": "pair", "metadata": {
                "label": "binder", "grade_semantics": "ordinal_enrichment", "evidence_score": 3}}) + "\n")
            paths["yeast.jsonl"].write_text(json.dumps({
                "pair_id": "pair", "protenix_json_path": str(prediction),
                "antibody_chain_ids": ["H"], "antigen_chain_ids": ["A"],
            }) + "\n")
            paths["igdesign.jsonl"].write_text("")
            paths["source.jsonl"].write_text("")
            paths["affinity.csv"].write_text("pair_id,label\n")
            paths["cdr.csv"].write_text("Id,seqstart_index,27,39\npair__chain_H,0,A,A\n")
            args = argparse.Namespace(records=paths["records.tsv"], yeast_panel=paths["panel.jsonl"],
                yeast_predictions=paths["yeast.jsonl"], igdesign_predictions=paths["igdesign.jsonl"],
                source_predictions=paths["source.jsonl"], source_affinity=paths["affinity.csv"],
                cdr_yeast=[paths["cdr.csv"]], cdr_igdesign=[paths["cdr.csv"]],
                cdr_source=[paths["cdr.csv"]], expected=1)
            candidate, = build(args)
            self.assertEqual(candidate.antibody_sequences, {"H": "AA"})
            self.assertEqual(candidate.cdr_mask, (True, False))
            self.assertEqual(candidate.evidence["scale"], "ordinal_enrichment")

    def test_igdesign_negative_with_moved_input(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            moved = root / "moved"
            moved.mkdir()
            (moved / "input.json").write_text(json.dumps([{"sequences": [
                {"id": "H", "proteinChain": {"sequence": "AA"}},
                {"id": "A", "proteinChain": {"sequence": "G"}},
            ]}]))
            (root / "records.tsv").write_text(
                "benchmark\tpair_id\ttarget\tlabel\tantigen\theavy\tlight\n"
                "IgDesign SPR\tpair\ttarget\tnegative\tG\tAA\t\n")
            (root / "igdesign.jsonl").write_text(json.dumps({
                "job_id": "pair", "record_id": "record", "label": "negative",
                "replicates": 3, "binding_replicates": 0,
                "protenix_json_path": "/previous/input.json",
                "antibody_chain_ids": ["H"], "antigen_chain_ids": ["A"],
            }) + "\n")
            (root / "cdr.csv").write_text("Id,seqstart_index,27,39\nrecord__chain_H,0,A,A\n")
            for name in ("panel.jsonl", "yeast.jsonl", "source.jsonl"):
                (root / name).write_text("")
            (root / "affinity.csv").write_text("pair_id,label\n")
            args = argparse.Namespace(
                records=root / "records.tsv", yeast_panel=root / "panel.jsonl",
                yeast_predictions=root / "yeast.jsonl", igdesign_predictions=root / "igdesign.jsonl",
                source_predictions=root / "source.jsonl", source_affinity=root / "affinity.csv",
                cdr_yeast=[root / "cdr.csv"], cdr_igdesign=[root / "cdr.csv"],
                cdr_source=[root / "cdr.csv"], path_map=[f"/previous={moved}"], expected=1,
            )
            candidate, = build(args)
            self.assertEqual(candidate.label, 0)
            self.assertEqual(candidate.evidence["binding_replicates"], 0)

    def test_source_heavy_and_light_annotations(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "input.json").write_text(json.dumps([{"sequences": [
                {"id": "H", "proteinChain": {"sequence": "AA"}},
                {"id": "L", "proteinChain": {"sequence": "GG"}},
                {"id": "A", "proteinChain": {"sequence": "T"}},
            ]}]))
            (root / "records.tsv").write_text(
                "benchmark\tpair_id\ttarget\tlabel\tantigen\theavy\tlight\n"
                "Literature/patent affinity\tpair\ttarget\tpositive\tT\tAA\tGG\n")
            (root / "source.jsonl").write_text(json.dumps({
                "pair_id": "pair", "protenix_json_path": str(root / "input.json"),
                "antibody_chain_ids": ["H", "L"], "antigen_chain_ids": ["A"],
            }) + "\n")
            (root / "affinity.csv").write_text(
                "pair_id,label,kd_median_nm,kd_measurement_count\npair,binder,4.5,2\n")
            (root / "heavy.csv").write_text("Id,seqstart_index,27,39\npair__chain_H,0,A,A\n")
            (root / "light.csv").write_text("Id,seqstart_index,27,39\npair__chain_L,0,G,G\n")
            for name in ("panel.jsonl", "yeast.jsonl", "igdesign.jsonl"):
                (root / name).write_text("")
            args = argparse.Namespace(
                records=root / "records.tsv", yeast_panel=root / "panel.jsonl",
                yeast_predictions=root / "yeast.jsonl", igdesign_predictions=root / "igdesign.jsonl",
                source_predictions=root / "source.jsonl", source_affinity=root / "affinity.csv",
                cdr_yeast=[root / "heavy.csv"], cdr_igdesign=[root / "heavy.csv"],
                cdr_source=[root / "heavy.csv", root / "light.csv"], expected=1,
            )
            candidate, = build(args)
            self.assertEqual(candidate.cdr_mask, (True, False, True, False))
            self.assertEqual(candidate.evidence["kd_median"], 4.5)


if __name__ == "__main__":
    unittest.main()
