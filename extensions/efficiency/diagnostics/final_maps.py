"""Four truth-coordinate efficiency maps after the final performance selection.

This additive module leaves the earlier diagnostic engine unchanged. Its selection
predicate is checked against the frozen original producer in the unit tests.
"""
import copy

import numpy as np

from .core import MASSES, REQUIRED_BRANCHES, _edges, validate_arrays, wilson
from .proton import policy, policy_metadata, proton_mask


PROTON_MASS_GEV = 0.9382720813
DEFAULT_Q2_EDGES = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4,
                    .001, .01, .03, .1, .3, 1.]
MAP_AXES = {"pt_q2": ("true_pT", "true_Q2"),
            "pt_y": ("true_pT", "true_y"),
            "xv_q2": ("true_xV", "true_Q2"),
            "xb_q2": ("true_xB", "true_Q2")}


def _performance_selection_masks(a, sample, settings):
    """Literal All States performance cuts, including the original EPZ predicate."""
    ee = float(settings["electron_energy"])
    lo, hi = settings["w_range"]
    truth = (np.isfinite(a["true_mass"])
             & (np.abs(a["true_mass"] - MASSES[sample["state"]]) <= settings["true_mass_window"])
             & np.isfinite(a["true_Q2"]) & (a["true_Q2"] >= 0)
             & (a["true_Q2"] <= sample["q2_max"])
             & np.isfinite(a["true_W"]) & (a["true_W"] >= lo)
             & (a["true_W"] <= hi))
    selected = (truth & np.isfinite(a["reco_mass"])
                & (a["reco_mass"] >= 8) & (a["reco_mass"] <= 12)
                & np.isfinite(a["reco_W"]) & (a["reco_W"] >= lo)
                & (a["reco_W"] <= hi)
                & np.isfinite(a["sigma_epz"])
                & (np.abs(a["sigma_epz"] - 2 * ee) <= .2 * ee))
    return truth, selected


def final_selection_masks(a, sample, settings):
    truth, pretag = _performance_selection_masks(a, sample, settings)
    return truth, pretag & proton_mask(a, settings, len(a["true_mass"]))


class FinalMapAccumulator:
    """Streaming integer counts; no weights, flux or luminosity enter the maps."""

    def __init__(self, sample, settings):
        self.sample, self.settings = copy.deepcopy(sample), copy.deepcopy(settings)
        self.proton_policy = policy(settings)
        if sample.get("state") not in MASSES:
            raise ValueError("Unsupported Upsilon state")
        for key in ("key", "id_mode"):
            if not sample.get(key):
                raise ValueError("Sample lacks " + key)
        for value, name in [(sample["q2_max"], "q2_max"),
                            (settings["electron_energy"], "electron_energy"),
                            (settings["true_mass_window"], "true_mass_window")]:
            if not np.isfinite(value) or value <= 0:
                raise ValueError(name + " must be finite and positive")
        self.w_range = _edges(settings["w_range"], "w_range")
        if len(self.w_range) != 2 or self.w_range[0] <= PROTON_MASS_GEV:
            raise ValueError("w_range requires two endpoints above the proton mass")
        self.coverage = float(settings["coverage"])
        wilson(0, 0, self.coverage)
        self.min_denominator = int(settings["min_denominator"])
        if self.min_denominator < 1 or self.min_denominator != settings["min_denominator"]:
            raise ValueError("min_denominator must be a positive integer")
        overrides = settings.get("final_2d", {})
        if not isinstance(overrides, dict):
            raise ValueError("settings.final_2d must be an object")
        unknown = set(overrides) - {"pt_edges", "y_edges", "q2_edges", "xb_edges", "xv_edges"}
        if unknown:
            raise ValueError("Unknown final_2d settings: " + ", ".join(sorted(unknown)))
        self.edges = {
            "true_pT": _edges(overrides.get("pt_edges", settings["pt_edges"]), "final_2d.pt_edges"),
            "true_y": _edges(overrides.get("y_edges", settings["y_edges"]), "final_2d.y_edges"),
            "true_Q2": _edges(overrides.get("q2_edges", DEFAULT_Q2_EDGES), "final_2d.q2_edges"),
            "true_xB": _edges(overrides.get("xb_edges", np.geomspace(1e-14, .003, 31)), "final_2d.xb_edges"),
            "true_xV": _edges(overrides.get("xv_edges", np.geomspace(.008, .3, 19)), "final_2d.xv_edges"),
        }
        if self.edges["true_pT"][0] < 0:
            raise ValueError("Final-map pT edges must be nonnegative")
        for variable in ("true_Q2", "true_xB", "true_xV"):
            if self.edges[variable][0] <= 0:
                raise ValueError(variable + " logarithmic edges must be strictly positive")
        cap = float(sample["q2_max"])
        if cap <= self.edges["true_Q2"][0]:
            raise ValueError("q2_max must exceed the first positive Q2 plotting edge")
        self.edges["true_Q2"] = np.append(self.edges["true_Q2"][self.edges["true_Q2"] < cap], cap)
        self.hist = {}
        self.flows = {}
        for name, axes in MAP_AXES.items():
            shape = tuple(len(self.edges[v]) - 1 for v in axes)
            self.hist[name] = {"denominator": np.zeros(shape, dtype=np.int64),
                               "selected": np.zeros(shape, dtype=np.int64),
                               "selected_before_proton": np.zeros(shape, dtype=np.int64)}
            self.flows[name] = {}
        self.raw_rows = self.chunks = self.eligible_truth = self.selected = self.selected_before_proton = 0
        self.zero_q2_truth = self.zero_q2_selected = 0
        self._finished = False
        self._result = None

    def consume(self, arrays):
        if self._finished:
            raise RuntimeError("Cannot consume after finish()")
        a = validate_arrays(arrays)
        nrows = len(a["true_mass"])
        self.raw_rows += nrows
        self.chunks += 1
        truth, pretag = _performance_selection_masks(a, self.sample, self.settings)
        selected = pretag & proton_mask(arrays, self.settings, nrows)
        self.selected_before_proton += int(pretag.sum())
        self.eligible_truth += int(truth.sum())
        self.selected += int(selected.sum())
        self.zero_q2_truth += int(np.count_nonzero(truth & (a["true_Q2"] == 0)))
        self.zero_q2_selected += int(np.count_nonzero(selected & (a["true_Q2"] == 0)))
        # Evaluate derived invariants in float64 even for Float_t ROOT branches.
        q2 = np.asarray(a["true_Q2"], dtype=np.float64)
        w = np.asarray(a["true_W"], dtype=np.float64)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            denominator = w ** 2 + q2 - PROTON_MASS_GEV ** 2
            usable = np.isfinite(denominator) & (denominator > 0)
            xb = np.full(nrows, np.nan, dtype=float)
            xv = np.full(nrows, np.nan, dtype=float)
            np.divide(q2, denominator, out=xb, where=usable)
            np.divide(q2 + MASSES[self.sample["state"]] ** 2,
                      denominator, out=xv, where=usable)
        coords = dict(a, true_xB=xb, true_xV=xv)
        for name, axes in MAP_AXES.items():
            values = [coords[v] for v in axes]
            edges = [self.edges[v] for v in axes]
            valid = np.ones(nrows, dtype=bool)
            below = np.zeros(nrows, dtype=bool)
            above = np.zeros(nrows, dtype=bool)
            nonpositive_log = np.zeros(nrows, dtype=bool)
            positive_log_underflow = np.zeros(nrows, dtype=bool)
            for variable, value, edge in zip(axes, values, edges):
                valid &= np.isfinite(value) & (np.abs(value) < 900)
                below |= value < edge[0]
                above |= value > edge[-1]
                if variable in ("true_Q2", "true_xB", "true_xV"):
                    nonpositive_log |= value <= 0
                    positive_log_underflow |= (value > 0) & (value < edge[0])
            outside = below | above
            for prefix, mask, field in (("truth", truth, "denominator"),
                                        ("selected", selected, "selected"),
                                        ("before_proton", pretag, "selected_before_proton")):
                inside = mask & valid & ~outside
                if inside.any():
                    self.hist[name][field] += np.histogramdd(
                        tuple(v[inside] for v in values), bins=edges)[0].astype(np.int64)
                counts = {
                    "total": int(mask.sum()),
                    "invalid_coordinates": int(np.count_nonzero(mask & ~valid)),
                    "valid_below_any_axis": int(np.count_nonzero(mask & valid & below)),
                    "valid_above_any_axis": int(np.count_nonzero(mask & valid & above)),
                    "valid_outside_any_axis": int(np.count_nonzero(mask & valid & outside)),
                    "in_histogram": int(inside.sum()),
                    "nonpositive_log_coordinates": int(np.count_nonzero(mask & valid & nonpositive_log)),
                    "positive_log_underflow": int(np.count_nonzero(mask & valid & positive_log_underflow)),
                }
                for label, count in counts.items():
                    key = prefix + "_" + label
                    self.flows[name][key] = self.flows[name].get(key, 0) + count

    def finish(self):
        if self._result is not None:
            return copy.deepcopy(self._result)
        if self.chunks == 0:
            self.consume({name: np.array([], dtype=float) for name in REQUIRED_BRANCHES + ("reco_topology", "reco_tagged_proton_p")})
        self._finished = True
        rows = []
        for name, axes in MAP_AXES.items():
            h, flow = self.hist[name], self.flows[name]
            if np.any(h["selected"] > h["selected_before_proton"]) or np.any(h["selected_before_proton"] > h["denominator"]):
                raise AssertionError("Selected counts exceed eligible truth")
            for prefix, field, total in (("truth", "denominator", self.eligible_truth),
                                         ("selected", "selected", self.selected),
                                         ("before_proton", "selected_before_proton", self.selected_before_proton)):
                if int(h[field].sum()) != flow[prefix + "_in_histogram"]:
                    raise AssertionError("Histogram totals disagree with flow audit")
                if flow[prefix + "_total"] != total:
                    raise AssertionError("Map totals differ from full selection totals")
                if (flow[prefix + "_invalid_coordinates"] + flow[prefix + "_valid_outside_any_axis"]
                        + flow[prefix + "_in_histogram"] != total):
                    raise AssertionError("Invalid/outside/in-range partition does not close")
            ex, ey = (self.edges[v] for v in axes)
            for i, j in np.ndindex(h["denominator"].shape):
                n, k = int(h["denominator"][i, j]), int(h["selected"][i, j])
                eff, low, high = wilson(k, n, self.coverage)
                rows.append(dict(sample=self.sample["key"], state=self.sample["state"],
                                 id_mode=self.sample["id_mode"], proton_selection=self.proton_policy, map=name,
                                 x_variable=axes[0], y_variable=axes[1],
                                 x_low=float(ex[i]), x_high=float(ex[i+1]),
                                 y_low=float(ey[j]), y_high=float(ey[j+1]),
                                 denominator=n, selected=k, selected_before_proton=int(h["selected_before_proton"][i,j]), selected_after_proton=k, efficiency=eff,
                                 wilson_low=low, wilson_high=high,
                                 sparse=bool(0 < n < self.min_denominator)))
        ee = float(self.settings["electron_energy"])
        mass = MASSES[self.sample["state"]]
        width = float(self.settings["true_mass_window"])
        audit = {
            "sample": copy.deepcopy(self.sample), "raw_rows": self.raw_rows, "chunks": self.chunks,
            "eligible_truth": self.eligible_truth, "selected": self.selected,
            "selected_before_proton": self.selected_before_proton, "selected_after_proton": self.selected,
            "rejected_by_proton": self.selected_before_proton-self.selected,
            "proton_selection": policy_metadata(self.settings),
            "zero_q2_eligible_truth": self.zero_q2_truth, "zero_q2_selected": self.zero_q2_selected,
            "status": "COMPLETE" if self.raw_rows else "EMPTY_INPUT",
            "map_names": list(MAP_AXES),
            "edges": {name: values.tolist() for name, values in self.edges.items()},
            "truth_selection": {"mass_GeV": [mass-width, mass+width],
                                "Q2_GeV2": [0., float(self.sample["q2_max"])],
                                "W_GeV": self.w_range.tolist()},
            "reco_selection": {"mass_GeV": [8., 12.], "W_GeV": self.w_range.tolist(),
                               "sigma_epz_GeV": [1.8*ee, 2.2*ee],
                               "sigma_epz_predicate": "abs(sigma_epz - 2*Ee) <= 0.2*Ee"},
            "x_definitions": {"xB": "true_Q2 / (true_W**2 + true_Q2 - proton_mass_GeV**2)",
                              "xV": "(true_Q2 + nominal_state_mass_GeV**2) / (true_W**2 + true_Q2 - proton_mass_GeV**2)",
                              "proton_mass_GeV": PROTON_MASS_GEV, "nominal_state_mass_GeV": mass,
                              "xV_interpretation": "Exclusive heavy-vector-meson small-|t| proxy; t is neglected, not measured Bjorken x.",
                              "mass_convention": "Fixed nominal state mass, not reconstructed or fluctuating true_mass."},
            "histogram_flows": [{"map": name, **self.flows[name]} for name in MAP_AXES],
            "histogram_convention": "Left-closed/right-open except final upper edge included; Q2 upper edge equals the sample cap.",
            "flow_convention": "Invalid coordinates, valid outside any axis, and in-histogram partition eligible truth/selected counts independently for each map. Below-any and above-any can overlap. Nonpositive-log and positive-log-underflow counts are diagnostic subsets, not additional partition terms.",
            "log_axis_convention": "No clipping or pseudocounts. Q2=0 stays in the full eligible denominator and pT-rapidity map where coordinates allow, but is outside logarithmic Q2 maps. Positive values below a plotting edge are reported as underflow.",
            "coordinate_convention": "Generator coordinates for both numerator and denominator. Stored finite coordinates require abs(value)<900. Negative pT is not repaired and lies below the nonnegative pT plotting range. x maps do not impose pT/rapidity validity or range restrictions.",
            "efficiency_definition": "Unweighted selected eligible stored UpsilonTree rows / eligible stored rows in the same truth-coordinate cell after the full performance selection.",
            "denominator_caveat": "Denominator is eligible stored tuple rows. Missing producer rows are not recovered or independently audited here.",
            "interval_definition": "Central Wilson score intervals for unweighted finite-MC binomial counts; zero successes retain a positive upper bound; zero denominator is undefined.",
            "coverage": self.coverage, "sparse_threshold": self.min_denominator,
            "weights_used": False, "luminosity_used": False, "flux_used": False,
            "count_invariants": {"selected_subset": "PASS", "histogram_flow_reconciliation": "PASS",
                                 "invalid_outside_inside_partition": "PASS", "map_total_reconciliation": "PASS"},
        }
        self._result = {"rows": rows, "audit": audit}
        return copy.deepcopy(self._result)
