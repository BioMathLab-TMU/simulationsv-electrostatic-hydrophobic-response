"""Small, self-contained integrity and statistical-contract tests."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "analysis"))
from force_statistics import paired_statistics, exact_sign_flip
from quality_control import evaluate_runtime_diagnostics


class ReleaseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest = json.loads((ROOT / "configs/manifest.json").read_text())
        cls.available = json.loads((ROOT / "data/availability.json").read_text())

    def test_registered_design_and_inputs(self):
        rows = self.manifest["runs"]
        self.assertEqual(len(rows), 12)
        self.assertEqual({r["seed"] for r in rows}, {9310002, 9310003, 9310004})
        self.assertEqual({(r["electrostatic_enabled"], r["hydrophobic_enabled"]) for r in rows},
                         {(False, False), (False, True), (True, False), (True, True)})
        for row in rows:
            config = ROOT / row["config"]
            self.assertEqual(hashlib.sha256(config.read_bytes()).hexdigest(), row["config_sha256"])
            payload = json.loads(config.read_text())
            forces = payload["forces_parms"][0]
            self.assertEqual(forces["enable_electrostatic"], row["electrostatic_enabled"])
            self.assertEqual(forces["enable_hydrophobic"], row["hydrophobic_enabled"])
            self.assertTrue(forces["enable_van_der_waals"])
            self.assertEqual(payload["times"][0]["time_step"], 0.0007)

    def test_released_data_integrity_and_schedule(self):
        self.assertEqual(self.available["verified_released_runs"], 12)
        for row in self.available["runs"]:
            folder = ROOT / "data/raw" / row["run_id"]
            h5path = folder / "trajectory.h5"
            self.assertEqual(hashlib.sha256(h5path.read_bytes()).hexdigest(), row["publication_h5_sha256"])
            self.assertEqual(evaluate_runtime_diagnostics(folder, 40000)["status"], "passed")
            with h5py.File(h5path, "r") as h5:
                frames = {name for name in h5 if name.startswith("step_")}
                self.assertEqual(frames, {f"step_{step}" for step in range(200, 40001, 200)})
                self.assertEqual(h5["step_200/positions"].shape, (2007, 3))
                self.assertEqual(h5["run_completion/run_completed"][()].item(), 1)
                self.assertNotIn("provenance", h5)
                self.assertNotIn("rng_state_mt19937_ascii", h5["step_200"])

    def test_contrasts_use_seed_not_frames(self):
        records = []
        for seed in (2, 3, 4):
            for e in (0, 1):
                for h in (0, 1):
                    records.append({"seed": seed, "electrostatic_enabled": e,
                                    "hydrophobic_enabled": h,
                                    "late_mean_largest_component_fraction": seed / 100 + 0.1*e + 0.2*h + 0.3*e*h,
                                    "late_mean_synapsin_recruited_fraction": 0.2 + 0.1*e + 0.1*h,
                                    "late_mean_largest_component_vesicle_fraction": 0.1 + 0.1*e + 0.2*h})
        _, contrasts = paired_statistics(pd.DataFrame(records))
        row = contrasts[(contrasts.metric == "late_mean_largest_component_fraction") &
                        (contrasts.contrast == "E_H_interaction")].iloc[0]
        self.assertEqual(row.n_seeds, 3)
        self.assertAlmostEqual(row["mean"], 0.3)
        self.assertEqual(exact_sign_flip([0.3, 0.3, 0.3]), 0.25)

    def test_processed_state_and_released_labels(self):
        processed = ROOT / "data/processed"
        status = json.loads((processed / "analysis_status.json").read_text())
        self.assertTrue(status["confirmatory_inference_ready"])
        frame = pd.read_csv(processed / "run_endpoints.csv")
        self.assertEqual(len(frame), 12)
        self.assertEqual(frame.groupby(["E", "H"]).size().tolist(), [3, 3, 3, 3])
        self.assertEqual(frame.frame_count.unique().tolist(), [200])


if __name__ == "__main__":
    unittest.main()
