"""Analytic/unit checks for the post hoc article diagnostics."""
import json
from pathlib import Path
import sys
import unittest
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"analysis"))
import interface_diagnostics
import mobility_diagnostic
from statistical_audit import holm


class DiagnosticTest(unittest.TestCase):
    def test_interface_exact_geometry(self):
        interface_diagnostics.self_test()

    def test_mobility_analytic_limits(self):
        mobility_diagnostic.self_test()

    def test_holm_preserves_original_order(self):
        np.testing.assert_allclose(holm([.04, .01, .03]), [.06, .03, .06])

    def test_reference_is_not_relabelled_as_model_validation(self):
        ref = json.loads((ROOT/"data/external/hoffmann_mobility_reference.json").read_text())
        self.assertEqual(ref["unit"], "um^2/s")
        self.assertEqual(ref["article_doi"], "10.1038/s41467-023-42372-6")
        self.assertGreater(ref["without_SV"]["mean"], ref["with_SV"]["mean"])
        self.assertGreaterEqual(len(ref["limitations"]), 4)

    def test_optional_generated_interface_reconciles_primary(self):
        path = ROOT/"results/interface_diagnostics/interface_run_metrics.csv"
        if not path.exists():
            self.skipTest("run interface_diagnostics.py to create numerical outputs")
        got = pd.read_csv(path).set_index("run_id").sort_index()
        old = pd.read_csv(ROOT/"data/processed/run_endpoints.csv").set_index("run_id").sort_index()
        self.assertEqual(list(got.index), list(old.index))
        for fresh, declared in [("largest_component_fraction", "late_mean_largest_component_fraction"),
              ("syn_nonisolated_fraction", "late_mean_synapsin_recruited_fraction"),
              ("largest_component_sv_fraction", "late_mean_largest_component_vesicle_fraction")]:
            np.testing.assert_allclose(got[fresh], old[declared], atol=1e-12, rtol=0)


if __name__ == "__main__":
    unittest.main()
