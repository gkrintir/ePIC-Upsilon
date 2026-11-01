"""Streaming, unweighted tuple-efficiency diagnostics.

The final selection deliberately matches the All States performance producer.
This module never reads event weights, rescales luminosity, repairs truth values,
or infers unrecorded producer/reconstruction decisions.  Numerators are subsets
of the same eligible stored rows used by their denominators.
"""
from statistics import NormalDist

import numpy as np

from .proton import policy, policy_metadata, proton_mask


MASSES = {"1S": 9.46030, "2S": 10.02326, "3S": 10.35520}
REQUIRED_BRANCHES = (
    "true_mass", "true_Q2", "true_W", "true_y", "true_pT",
    "reco_mass", "reco_W", "sigma_epz",
)
CUTS = ("reco_mass", "reco_W", "sigma_epz")
FAILURE_REASONS = ("nonfinite", "minus_999_sentinel", "finite_lower_rejection",
                   "finite_upper_rejection", "pass")


def wilson(selected, denominator, coverage=0.6826894921370859):
    """Return efficiency and central Wilson score limits for integer MC counts.

    An empty denominator has no efficiency or interval. Zero selected events
    with a positive denominator retain their nonzero upper confidence limit.
    """
    n, k = int(denominator), int(selected)
    if n != denominator or k != selected or n < 0 or not 0 <= k <= n:
        raise ValueError("Wilson counts must be integers with 0 <= selected <= denominator")
    if not np.isfinite(coverage) or not 0 < coverage < 1:
        raise ValueError("coverage must lie strictly between zero and one")
    if n == 0:
        return None, None, None
    z = NormalDist().inv_cdf((1.0 + float(coverage)) / 2.0)
    p = k / n
    divisor = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / divisor
    half = z * np.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / divisor
    return float(p), float(max(0.0, center - half)), float(min(1.0, center + half))


def _edges(values, name):
    result = np.asarray(values, dtype=float)
    if result.ndim != 1 or result.size < 2 or not np.all(np.isfinite(result)) or not np.all(np.diff(result) > 0):
        raise ValueError("%s must contain finite, strictly increasing edges" % name)
    return result


def validate_arrays(arrays):
    missing = [name for name in REQUIRED_BRANCHES if name not in arrays]
    if missing:
        raise ValueError("Missing required branches: " + ", ".join(missing))
    result = {}
    sizes = set()
    for name in REQUIRED_BRANCHES:
        value = np.asarray(arrays[name])
        if value.ndim != 1 or value.dtype.kind not in "iuf":
            raise ValueError("%s must be a one-dimensional numeric branch, got shape=%s dtype=%s" % (name, value.shape, value.dtype))
        result[name] = value
        sizes.add(len(value))
    if len(sizes) != 1:
        raise ValueError("Required branches do not have equal row counts")
    return result


class Accumulator:
    """Accumulate all diagnostics without retaining the event arrays in memory."""

    def __init__(self, sample, settings):
        self.sample = dict(sample)
        self.settings = dict(settings)
        self.proton_policy = policy(settings)
        for name in ("key", "state", "id_mode", "q2_max"):
            if name not in sample:
                raise ValueError("Sample lacks %s" % name)
        if sample["state"] not in MASSES:
            raise ValueError("Unknown Upsilon state: %s" % sample["state"])
        self.q2_max = float(sample["q2_max"])
        self.split = float(settings["q2_split"])
        self.ee = float(settings["electron_energy"])
        self.mass_window = float(settings["true_mass_window"])
        self.coverage = float(settings["coverage"])
        self.min_denominator = int(settings["min_denominator"])
        self.w_range = np.asarray(settings["w_range"], dtype=float)
        if self.w_range.shape != (2,) or not np.all(np.isfinite(self.w_range)) or self.w_range[0] >= self.w_range[1]:
            raise ValueError("w_range must contain two ordered finite values")
        if not all(np.isfinite(x) and x > 0 for x in (self.q2_max, self.split, self.ee, self.mass_window)):
            raise ValueError("Q2 maxima/split, electron energy and truth mass window must be positive")
        if self.min_denominator < 1 or self.min_denominator != settings["min_denominator"]:
            raise ValueError("min_denominator must be a positive integer")
        wilson(0, 0, self.coverage)
        self.edges = {
            "true_pT": _edges(settings["pt_edges"], "pt_edges"),
            "true_W": _edges(settings["w_edges"], "w_edges"),
            "true_y": _edges(settings["y_edges"], "y_edges"),
            "true_Q2": _edges(settings["q2_edges"], "q2_edges"),
        }
        self.tail_edges = _edges(settings["tail_edges"], "tail_edges")
        self.strata_w = _edges(settings["strata_w_edges"], "strata_w_edges")
        self.strata_y = _edges(settings["strata_y_edges"], "strata_y_edges")
        self.bands = ("all", "low", "high") if self.q2_max > self.split else ("all", "low")
        self.hist = {}
        self.flows = {}
        self.band_counts = {b: {"eligible_truth": 0, "selected": 0, "selected_before_proton": 0, "selected_after_proton": 0} for b in self.bands}
        self.raw_rows = 0
        self.chunks = 0
        self.invalid_input = {n: {"nonfinite": 0, "finite_abs_ge_900": 0, "negative_finite": 0, "exact_minus_999": 0} for n in REQUIRED_BRANCHES}
        self.truth_rejections = {n: 0 for n in ("mass", "q2", "w")}
        self._finished = False

    def _entry(self, band, **values):
        return dict(sample=str(self.sample["key"]), state=self.sample["state"],
                    id_mode=self.sample["id_mode"], q2_band=band, proton_selection=self.proton_policy, **values)

    def _estimate(self, k, n):
        estimate, low, high = wilson(int(k), int(n), self.coverage)
        return dict(efficiency=estimate, wilson_low=low, wilson_high=high,
                    sparse=bool(n < self.min_denominator))

    @staticmethod
    def _histogram(values, edges, mask):
        shape = tuple(len(e) - 1 for e in edges)
        if not np.any(mask):
            return np.zeros(shape, dtype=np.int64)
        if len(values) == 1:
            return np.histogram(values[0][mask], bins=edges[0])[0].astype(np.int64)
        return np.histogramdd(tuple(v[mask] for v in values), bins=edges)[0].astype(np.int64)

    def _add_pair(self, band, name, values, edges, truth, selected, record_flow=True):
        valid = np.ones(len(truth), dtype=bool)
        below = np.zeros(len(truth), dtype=bool)
        above = np.zeros(len(truth), dtype=bool)
        for value, edge in zip(values, edges):
            valid &= np.isfinite(value) & (np.abs(value) < 900)
            below |= value < edge[0]
            above |= value > edge[-1]
        key = (band, name)
        if key not in self.hist:
            shape = tuple(len(e) - 1 for e in edges)
            self.hist[key] = {"denominator": np.zeros(shape, dtype=np.int64),
                              "selected": np.zeros(shape, dtype=np.int64)}
        self.hist[key]["denominator"] += self._histogram(values, edges, truth & valid)
        self.hist[key]["selected"] += self._histogram(values, edges, selected & valid)
        if record_flow:
            flow = self.flows.setdefault(key, {})
            for prefix, mask in (("truth", truth), ("selected", selected)):
                counts = {
                    "total": int(np.count_nonzero(mask)),
                    "invalid_coordinates": int(np.count_nonzero(mask & ~valid)),
                    "valid_below_any_axis": int(np.count_nonzero(mask & valid & below)),
                    "valid_above_any_axis": int(np.count_nonzero(mask & valid & above)),
                    "valid_outside_any_axis": int(np.count_nonzero(mask & valid & (below | above))),
                    "in_histogram": int(np.count_nonzero(mask & valid & ~below & ~above)),
                }
                for suffix, count in counts.items():
                    label = prefix + "_" + suffix
                    flow[label] = flow.get(label, 0) + count

    def consume(self, arrays):
        if self._finished:
            raise RuntimeError("Cannot consume more rows after finish()")
        a = validate_arrays(arrays)
        size = len(a["true_mass"])
        proton_pass = proton_mask(arrays, self.settings, size)
        self.raw_rows += size
        self.chunks += 1
        for name, value in a.items():
            finite = np.isfinite(value)
            counter = self.invalid_input[name]
            counter["nonfinite"] += int(np.count_nonzero(~finite))
            counter["finite_abs_ge_900"] += int(np.count_nonzero(finite & (np.abs(value) >= 900)))
            counter["negative_finite"] += int(np.count_nonzero(finite & (value < 0)))
            counter["exact_minus_999"] += int(np.count_nonzero(value == -999))
        mass_ok = np.isfinite(a["true_mass"]) & (np.abs(a["true_mass"] - MASSES[self.sample["state"]]) <= self.mass_window)
        q2_ok = np.isfinite(a["true_Q2"]) & (a["true_Q2"] >= 0) & (a["true_Q2"] <= self.q2_max)
        w_ok = np.isfinite(a["true_W"]) & (a["true_W"] >= self.w_range[0]) & (a["true_W"] <= self.w_range[1])
        truth_all = mass_ok & q2_ok & w_ok
        for name, mask in (("mass", mass_ok), ("q2", q2_ok), ("w", w_ok)):
            self.truth_rejections[name] += int(np.count_nonzero(~mask))
        cuts = {
            "reco_mass": np.isfinite(a["reco_mass"]) & (a["reco_mass"] >= 8) & (a["reco_mass"] <= 12),
            "reco_W": np.isfinite(a["reco_W"]) & (a["reco_W"] >= self.w_range[0]) & (a["reco_W"] <= self.w_range[1]),
            "sigma_epz": np.isfinite(a["sigma_epz"]) & (np.abs(a["sigma_epz"] - 2 * self.ee) <= 0.2 * self.ee),
        }
        centers = {"reco_mass": 10.0, "reco_W": float(np.mean(self.w_range)), "sigma_epz": 2.0 * self.ee}
        reasons = {}
        for cut in CUTS:
            finite = np.isfinite(a[cut])
            sentinel = finite & (a[cut] == -999)
            rejected = finite & ~sentinel & ~cuts[cut]
            reasons[cut] = {
                "nonfinite": ~finite,
                "minus_999_sentinel": sentinel,
                "finite_lower_rejection": rejected & (a[cut] <= centers[cut]),
                "finite_upper_rejection": rejected & (a[cut] > centers[cut]),
                "pass": cuts[cut],
            }
        pretag_selection = cuts["reco_mass"] & cuts["reco_W"] & cuts["sigma_epz"]
        selection = pretag_selection & proton_pass
        masks = {"all": truth_all, "low": truth_all & (a["true_Q2"] <= self.split)}
        if "high" in self.bands:
            masks["high"] = truth_all & (a["true_Q2"] > self.split)
        for band in self.bands:
            truth = masks[band]
            selected = truth & selection
            pretag = truth & pretag_selection
            self.band_counts[band]["selected_before_proton"] += int(np.count_nonzero(pretag))
            self.band_counts[band]["selected_after_proton"] += int(np.count_nonzero(selected))
            self.band_counts[band]["eligible_truth"] += int(np.count_nonzero(truth))
            self.band_counts[band]["selected"] += int(np.count_nonzero(selected))
            for variable in ("true_pT", "true_W", "true_y"):
                self._add_pair(band, "efficiency/" + variable, [a[variable]], [self.edges[variable]], truth, selected)
                self._add_pair(band, "inclusive_efficiency/" + variable, [a[variable]], [self.edges[variable]], truth, pretag)
            entered = truth
            self._add_pair(band, "cutflow/truth", [a["true_pT"]], [self.edges["true_pT"]], truth, truth, False)
            for stage in CUTS:
                passed = entered & cuts[stage]
                self._add_pair(band, "cutflow/" + stage, [a["true_pT"]], [self.edges["true_pT"]], entered, passed, False)
                entered = passed
                other = truth & proton_pass
                for cut in CUTS:
                    if cut != stage:
                        other &= cuts[cut]
                self._add_pair(band, "nminus1/" + stage, [a["true_pT"]], [self.edges["true_pT"]], other, selected, False)
                for reason in FAILURE_REASONS:
                    self._add_pair(band, "failure_reasons/" + stage + "/" + reason,
                                   [a["true_pT"]], [self.edges["true_pT"]], truth,
                                   truth & reasons[stage][reason], False)
            self._add_pair(band, "cutflow/proton_tag", [a["true_pT"]], [self.edges["true_pT"]], entered, selected, False)
            self._add_pair(band, "proton_selection", [a["true_pT"]], [self.edges["true_pT"]], pretag, selected, False)
            for variable in ("true_W", "true_y", "true_Q2"):
                values = [a["true_pT"], a[variable]]
                self._add_pair(band, "maps/" + variable, values, [self.edges["true_pT"], self.edges[variable]], truth, selected)
                self._add_pair(band, "composition/" + variable, values, [self.tail_edges, self.edges[variable]], truth, selected)
            self._add_pair(band, "strata", [a["true_pT"], a["true_W"], a["true_y"]],
                           [self.edges["true_pT"], self.strata_w, self.strata_y], truth, selected)

    def _check_invariants(self):
        checks = []
        for (band, name), h in self.hist.items():
            if np.any(h["selected"] > h["denominator"]) or np.any(h["selected"] < 0):
                raise AssertionError("Selected counts are not a subset: %s/%s" % (band, name))
            flow = self.flows.get((band, name))
            if flow:
                for prefix, key in (("truth", "denominator"), ("selected", "selected")):
                    if int(h[key].sum()) != flow[prefix + "_in_histogram"]:
                        raise AssertionError("Histogram/flow disagreement: %s/%s" % (band, name))
                    total = flow[prefix + "_invalid_coordinates"] + flow[prefix + "_valid_outside_any_axis"] + flow[prefix + "_in_histogram"]
                    if total != flow[prefix + "_total"]:
                        raise AssertionError("Invalid/overflow partition disagreement")
            if band == "all":
                for kind in ("denominator", "selected"):
                    partition = self.hist[("low", name)][kind].copy()
                    if "high" in self.bands:
                        partition += self.hist[("high", name)][kind]
                    if not np.array_equal(h[kind], partition):
                        raise AssertionError("Q2 partition disagreement: %s/%s" % (name, kind))
        for kind in ("eligible_truth", "selected", "selected_before_proton", "selected_after_proton"):
            partition = self.band_counts["low"][kind] + (self.band_counts["high"][kind] if "high" in self.bands else 0)
            if self.band_counts["all"][kind] != partition:
                raise AssertionError("Aggregate Q2 partition disagreement")
        for band in self.bands:
            previous = self.hist[(band, "cutflow/truth")]["selected"]
            for stage in CUTS + ("proton_tag",):
                h = self.hist[(band, "cutflow/" + stage)]
                if not np.array_equal(previous, h["denominator"]):
                    raise AssertionError("Cutflow stages are not nested")
                previous = h["selected"]
            final = self.hist[(band, "efficiency/true_pT")]["selected"]
            if not np.array_equal(previous, final):
                raise AssertionError("Cutflow final selection differs from efficiency numerator")
            for stage in CUTS:
                if not np.array_equal(self.hist[(band, "nminus1/" + stage)]["selected"], final):
                    raise AssertionError("N-minus-one final numerator disagreement")
                partition = sum((self.hist[(band, "failure_reasons/" + stage + "/" + reason)]["selected"]
                                 for reason in FAILURE_REASONS), np.zeros_like(final))
                if not np.array_equal(partition, self.hist[(band, "efficiency/true_pT")]["denominator"]):
                    raise AssertionError("Stored-field failure categories do not partition eligible truth")
        checks.extend(["selected_counts_are_subsets", "histograms_reconcile_with_flow_counts",
                       "invalid_and_overflow_partition_reconciles", "all_equals_low_plus_high",
                       "sequential_cutflow_reconciles", "nminus1_final_numerators_reconcile",
                       "stored_field_failure_categories_partition_truth"])
        return {name: "PASS" for name in checks}

    def finish(self):
        # Empty files still receive all output rows and an explicit empty audit.
        if not self.hist:
            self.consume({name: np.array([], dtype=float) for name in REQUIRED_BRANCHES + ("reco_topology", "reco_tagged_proton_p")})
        self._finished = True
        invariants = self._check_invariants()
        tables = {name: [] for name in ("efficiency", "inclusive_efficiency", "proton_selection", "cutflow", "nminus1", "failure_reasons", "maps", "composition", "strata", "q2_mixture")}
        for band in self.bands:
            for variable in ("true_pT", "true_W", "true_y"):
                h, edges = self.hist[(band, "efficiency/" + variable)], self.edges[variable]
                for i, (n, k) in enumerate(zip(h["denominator"], h["selected"])):
                    tables["efficiency"].append(self._entry(band, variable=variable, low=float(edges[i]), high=float(edges[i+1]),
                        denominator=int(n), selected=int(k), **self._estimate(k, n)))
            for variable in ("true_pT", "true_W", "true_y"):
                h, edges = self.hist[(band, "inclusive_efficiency/" + variable)], self.edges[variable]
                for i, (n, k) in enumerate(zip(h["denominator"], h["selected"])):
                    tables["inclusive_efficiency"].append(self._entry(band, variable=variable, low=float(edges[i]), high=float(edges[i+1]), denominator=int(n), selected=int(k), comparison_category="inclusive_before_proton_selection", **self._estimate(k, n)))
            h = self.hist[(band, "proton_selection")]
            for i, (n, k) in enumerate(zip(h["denominator"], h["selected"])):
                tables["proton_selection"].append(self._entry(band, pt_low=float(self.edges["true_pT"][i]), pt_high=float(self.edges["true_pT"][i+1]), eligible_truth=int(self.hist[(band, "cutflow/truth")]["denominator"][i]), selected_before_proton=int(n), selected_after_proton=int(k), rejected_by_proton=int(n-k), tag_retention=float(k/n) if n else None, yield_convention="unweighted MC rows"))
                tables["nminus1"].append(self._entry(band, pt_low=float(self.edges["true_pT"][i]), pt_high=float(self.edges["true_pT"][i+1]), cut="proton_tag", denominator=int(n), selected=int(k), failed_without_cut=int(n-k), **self._estimate(k, n)))
            base = self.hist[(band, "cutflow/truth")]["denominator"]
            for stage in ("truth",) + CUTS + ("proton_tag",):
                h = self.hist[(band, "cutflow/" + stage)]
                for i, (n, k) in enumerate(zip(h["denominator"], h["selected"])):
                    tables["cutflow"].append(self._entry(band, pt_low=float(self.edges["true_pT"][i]), pt_high=float(self.edges["true_pT"][i+1]),
                        stage=stage, baseline=int(base[i]), entered=int(n), passed=int(k), failed=int(n-k),
                        step_efficiency=float(k/n) if n else None, **self._estimate(k, base[i])))
            for stage in CUTS:
                h = self.hist[(band, "nminus1/" + stage)]
                for i, (n, k) in enumerate(zip(h["denominator"], h["selected"])):
                    tables["nminus1"].append(self._entry(band, pt_low=float(self.edges["true_pT"][i]), pt_high=float(self.edges["true_pT"][i+1]),
                        cut=stage, denominator=int(n), selected=int(k), failed_without_cut=int(n-k), **self._estimate(k, n)))
                for reason in FAILURE_REASONS:
                    h = self.hist[(band, "failure_reasons/" + stage + "/" + reason)]
                    for i, (n, k) in enumerate(zip(h["denominator"], h["selected"])):
                        fraction, low, high = wilson(int(k), int(n), self.coverage)
                        tables["failure_reasons"].append(self._entry(band,
                            pt_low=float(self.edges["true_pT"][i]), pt_high=float(self.edges["true_pT"][i+1]),
                            cut=stage, reason=reason, baseline=int(n), denominator=int(n), count=int(k),
                            fraction=fraction, wilson_low=low, wilson_high=high,
                            sparse=bool(n < self.min_denominator)))
            for variable in ("true_W", "true_y", "true_Q2"):
                h = self.hist[(band, "maps/" + variable)]
                for i, j in np.ndindex(h["denominator"].shape):
                    n, k = int(h["denominator"][i, j]), int(h["selected"][i, j])
                    tables["maps"].append(self._entry(band, x_variable="true_pT", y_variable=variable,
                        x_low=float(self.edges["true_pT"][i]), x_high=float(self.edges["true_pT"][i+1]),
                        y_low=float(self.edges[variable][j]), y_high=float(self.edges[variable][j+1]),
                        denominator=n, selected=k, **self._estimate(k, n)))
                h = self.hist[(band, "composition/" + variable)]
                den_norm, sel_norm = h["denominator"].sum(axis=1), h["selected"].sum(axis=1)
                for i, j in np.ndindex(h["denominator"].shape):
                    n, k = int(h["denominator"][i, j]), int(h["selected"][i, j])
                    tables["composition"].append(self._entry(band, pt_low=float(self.tail_edges[i]), pt_high=float(self.tail_edges[i+1]),
                        variable=variable, low=float(self.edges[variable][j]), high=float(self.edges[variable][j+1]),
                        denominator=n, selected=k, truth_fraction=float(n/den_norm[i]) if den_norm[i] else None,
                        selected_fraction=float(k/sel_norm[i]) if sel_norm[i] else None,
                        truth_normalization_count=int(den_norm[i]), selected_normalization_count=int(sel_norm[i]),
                        normalization="valid in-range variable counts within pT region"))
            h = self.hist[(band, "strata")]
            for i, j, l in np.ndindex(h["denominator"].shape):
                n, k = int(h["denominator"][i, j, l]), int(h["selected"][i, j, l])
                tables["strata"].append(self._entry(band, pt_low=float(self.edges["true_pT"][i]), pt_high=float(self.edges["true_pT"][i+1]),
                    w_low=float(self.strata_w[j]), w_high=float(self.strata_w[j+1]), y_low=float(self.strata_y[l]), y_high=float(self.strata_y[l+1]),
                    denominator=n, selected=k, **self._estimate(k, n)))
        if "high" in self.bands:
            for variable in ("true_pT", "true_W", "true_y"):
                histograms = {band: self.hist[(band, "efficiency/" + variable)] for band in self.bands}
                edges = self.edges[variable]
                for i in range(len(edges) - 1):
                    n = {band: int(h["denominator"][i]) for band, h in histograms.items()}
                    k = {band: int(h["selected"][i]) for band, h in histograms.items()}
                    eff = {band: k[band] / n[band] if n[band] else None for band in self.bands}
                    low_contribution = k["low"] / n["all"] if n["all"] else None
                    high_contribution = k["high"] / n["all"] if n["all"] else None
                    residual = eff["all"] - low_contribution - high_contribution if n["all"] else None
                    if residual is not None and abs(residual) > 5e-15:
                        raise AssertionError("Q2 mixture decomposition does not close")
                    tables["q2_mixture"].append(self._entry("all", variable=variable,
                        low=float(edges[i]), high=float(edges[i+1]),
                        n_all=n["all"], n_low=n["low"], n_high=n["high"],
                        k_all=k["all"], k_low=k["low"], k_high=k["high"],
                        high_truth_fraction=n["high"] / n["all"] if n["all"] else None,
                        high_selected_fraction=k["high"] / k["all"] if k["all"] else None,
                        eff_all=eff["all"], eff_low=eff["low"], eff_high=eff["high"],
                        low_contribution=low_contribution, high_contribution=high_contribution,
                        mixture_residual=residual))
            invariants["q2_mixture_decomposition_reconciles"] = "PASS"
        ranges = {"all": {"low": 0.0, "high": self.q2_max, "lower_inclusive": True, "upper_inclusive": True},
                  "low": {"low": 0.0, "high": min(self.split, self.q2_max), "lower_inclusive": True, "upper_inclusive": True}}
        if "high" in self.bands:
            ranges["high"] = {"low": self.split, "high": self.q2_max, "lower_inclusive": False, "upper_inclusive": True}
        audit = {
            "sample": {k: str(v) if k == "path" else v for k, v in self.sample.items()},
            "proton_selection": policy_metadata(self.settings),
            "raw_rows": self.raw_rows, "chunks": self.chunks, "band_counts": self.band_counts,
            "q2_band_definitions_GeV2": ranges, "count_invariants": invariants,
            "input_value_counts": self.invalid_input,
            "truth_failed_requirements_independent_overlapping": self.truth_rejections,
            "histogram_flows": [{"q2_band": band, "histogram": name, **flow} for (band, name), flow in sorted(self.flows.items())],
            "histogram_convention": "Left-closed/right-open bins except final upper edge included. Q2 band boundaries are separately explicit.",
            "flow_convention": "Invalid, valid outside, and in-histogram partition the eligible count. Below-any and above-any can overlap in multidimensional maps.",
            "truth_selection": {"mass_window_GeV": [MASSES[self.sample["state"]]-self.mass_window, MASSES[self.sample["state"]]+self.mass_window],
                                "Q2_GeV2": [0.0, self.q2_max], "W_GeV": self.w_range.tolist()},
            "reco_selection": {"mass_GeV": [8.0, 12.0], "W_GeV": self.w_range.tolist(), "sigma_epz_GeV": [1.8*self.ee, 2.2*self.ee]},
            "cutflow_definition": "Sequential cuts ordered reco_mass, reco_W, sigma_epz, proton_tag. The final stage applies the explicit reconstructed-proton policy. These are not a reconstruction-algorithm cutflow.",
            "nminus1_definition": "For each original stored-field cut, denominator passes truth, the proton policy and all other stored-field cuts; failed_without_cut is the number failing that omitted cut.",
            "failure_reasons_definition": "For each stored field independently, nonfinite / exact -999 / other finite lower rejection / other finite upper rejection / pass partition all eligible truth. Pass uses the original cut predicate; other finite rejected values are split around the cut center. Missing or sentinel mass is not evidence for a mass-window or ID-specific loss. These categories are not sequential producer stages.",
            "efficiency_definition": "Unweighted selected eligible stored UpsilonTree rows / eligible stored rows in the same truth-variable bin.",
            "denominator_caveat": "No independently verified full generated-event denominator; absent producer rows are not recovered.",
            "interval_definition": "Central Wilson score interval for finite unweighted MC binomial counts; no luminosity or projected-yield scaling.",
            "q2_mixture_definition": "For samples extending above the Q2 split, eff_all = k_low/n_all + k_high/n_all. Contributions describe the stored sample composition, not an isolated causal detector effect. No difference significance is inferred.",
            "coverage": self.coverage, "weights_used": False, "luminosity_used": False,
            "sparse_threshold": self.min_denominator,
            "sparse_definition": "Denominator below the descriptive threshold; bins are retained, never repaired or merged automatically.",
            "status": "COMPLETE" if self.raw_rows else "EMPTY_INPUT",
        }
        return {"tables": tables, "audit": audit}
