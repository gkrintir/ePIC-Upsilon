"""Proton selection is explicit, schema-checked, and numerator-only."""
import copy
import unittest

import numpy as np

from diagnostics.core import Accumulator
from diagnostics.final_maps import FinalMapAccumulator
from diagnostics.proton import policy, proton_mask
from run_efficiency_diagnostics import required_branches
from test_core import arrays, settings, sample


def tagged_arrays():
    a = arrays(5)
    a["reco_topology"] = np.array([0, 1, 2, 3, 4])
    a["reco_tagged_proton_p"] = np.array([-999., 250., 230., -999., -999.])
    a["reco_mass"][0] = -999.
    return a


def cfg(mode):
    result = settings()
    if mode is None:
        result.pop("proton_selection")
    else:
        result["proton_selection"] = mode
    return result


class ProtonSelectionTests(unittest.TestCase):
    def test_default_is_required_and_adds_two_branches(self):
        self.assertEqual(policy({}), "required")
        self.assertEqual(required_branches({})[-2:], ("reco_topology", "reco_tagged_proton_p"))
        self.assertEqual(len(required_branches({"proton_selection": "inclusive"})), 8)

    def test_missing_schema_fails_both_accumulators(self):
        for cls in (Accumulator, FinalMapAccumulator):
            with self.assertRaisesRegex(ValueError, "requires branches"):
                cls(sample(), cfg(None)).consume(arrays())

    def test_known_topologies_define_mask(self):
        np.testing.assert_array_equal(proton_mask(tagged_arrays(), cfg(None), 5), [False, True, True, False, False])

    def test_inconsistent_or_unknown_codes_fail_even_outside_truth(self):
        for topology, p in [(1, -999), (2, 0), (1, np.nan), (2, np.inf),
                            (0, 240), (3, 0), (4, np.nan), (5, -999),
                            (1.5, 240), (np.nan, -999), (np.inf, -999)]:
            a = arrays(1)
            a["true_mass"][:] = -999
            a["reco_topology"] = np.array([topology])
            a["reco_tagged_proton_p"] = np.array([p])
            for cls in (Accumulator, FinalMapAccumulator):
                with self.assertRaises(ValueError, msg=str((topology, p))):
                    cls(sample(), cfg(None)).consume(a)

    def test_inclusive_bypasses_schema_and_matches_old_numerator(self):
        a = tagged_arrays()
        del a["reco_topology"]
        del a["reco_tagged_proton_p"]
        for cls in (Accumulator, FinalMapAccumulator):
            acc = cls(sample(), cfg("inclusive"))
            acc.consume(a)
            result = acc.finish()
            audit = result["audit"]
            counts = audit.get("band_counts", {}).get("all", audit)
            self.assertEqual(counts["eligible_truth"], 5)
            self.assertEqual(counts["selected"], 4)
            self.assertEqual(audit["proton_selection"]["policy"], "inclusive")

    def test_truth_denominator_unchanged_and_numerator_smaller(self):
        a = tagged_arrays()
        for cls in (Accumulator, FinalMapAccumulator):
            results = {}
            for mode in ("inclusive", "required"):
                acc = cls(sample(), cfg(mode))
                acc.consume(a)
                results[mode] = acc.finish()
            old, new = (results[m] for m in ("inclusive", "required"))
            for result in (old, new):
                audit = result["audit"]
                counts = audit.get("band_counts", {}).get("all", audit)
                self.assertEqual(counts["eligible_truth"], 5)
                self.assertEqual(counts["selected_before_proton"], 4)
            audit = new["audit"]
            counts = audit.get("band_counts", {}).get("all", audit)
            self.assertEqual(counts["selected_after_proton"], 2)
            if cls is FinalMapAccumulator:
                for before, after in zip(old["rows"], new["rows"]):
                    self.assertEqual(before["denominator"], after["denominator"])
                    self.assertEqual(before["selected"], after["selected_before_proton"])
                    self.assertLessEqual(after["selected"], before["selected"])
            else:
                for before, after in zip(old["tables"]["efficiency"], new["tables"]["efficiency"]):
                    self.assertEqual(before["denominator"], after["denominator"])
                    self.assertLessEqual(after["selected"], before["selected"])
                cuts = [r for r in new["tables"]["cutflow"] if r["q2_band"] == "all" and r["stage"] == "proton_tag"]
                self.assertEqual(sum(r["entered"] for r in cuts), 4)
                self.assertEqual(sum(r["passed"] for r in cuts), 2)
                self.assertTrue(all(v == "PASS" for v in audit["count_invariants"].values()))

    def test_tagged_chunking_and_weight_independence(self):
        a = tagged_arrays()
        for cls in (Accumulator, FinalMapAccumulator):
            whole = cls(sample(), cfg(None))
            whole.consume(a)
            segmented = cls(sample(), cfg(None))
            for sl in (slice(0,2), slice(2,5)):
                part = {k:v[sl] for k,v in a.items()}
                part["weight"] = np.full(len(part["true_mass"]), np.nan)
                segmented.consume(part)
            left, right = whole.finish(), segmented.finish()
            field = "rows" if cls is FinalMapAccumulator else "tables"
            self.assertEqual(left[field], right[field])


if __name__ == "__main__":
    unittest.main()
