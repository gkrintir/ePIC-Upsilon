"""Read-only yield and photon-flux accounting for the Upsilon W extraction.

Only NumPy and the Python standard library are required. This module performs no
I/O, weight normalization, fit, event generation, or uncertainty reinterpretation.
It reproduces the current full-sample bin-by-bin central-value selections. The
real-photon eSTARlight function is a current-source benchmark, not certification
of a particular production campaign or its finite-Q2 generator calculation.
"""
from functools import lru_cache
import os

import numpy as np

ALPHA_EM = 1.0 / 137.035999084
M_E = 0.00051099895
M_P = 0.9382720813
MODEL_UPSILON_MASS = 9.46030
MODEL_SOURCE = "https://eic-code-browser.sdcc.bnl.gov/lxr/source/estarlight/src/photonNucleusCrossSection.cpp"
REQUIRED_BRANCHES = (
    "weight", "true_mass", "true_W", "reco_mass", "reco_W", "sigma_epz", "true_Q2",
)
CUTFLOW_STAGES = (
    "all", "truthmass", "truthmassQ2", "truthmassQ2W",
    "reco_mass", "reco_epz", "reco_Q2", "reco_W_before_proton", "reco_W",
)


def proton_policy(settings=None):
    """Read the declared category without depending on the framework selector."""
    policy = (settings or {}).get("proton_selection", os.environ.get("YPSILON_PROTON_SELECTION", "required"))
    if isinstance(policy, dict):
        policy = policy.get("mode")
    if policy not in ("required", "inclusive"):
        raise ValueError("YPSILON_PROTON_SELECTION must be required or inclusive")
    return policy


def required_branches(settings=None):
    return REQUIRED_BRANCHES + (() if proton_policy(settings) == "inclusive" else
                                ("reco_topology", "reco_tagged_proton_p"))


def independent_proton_mask(data, policy):
    """Independently check the producer topology/momentum contract, then select.

    The producer writes topology 1/2 for tagged valid lepton combinations,
    3/4 for untagged combinations and 0 before a valid combination exists.
    No reconstructed observable is used as a substitute for a proton tag.
    """
    rows = len(data["weight"])
    if policy == "inclusive":
        return np.ones(rows, dtype=bool)
    topology = np.asarray(data["reco_topology"], dtype=float)
    momentum = np.asarray(data["reco_tagged_proton_p"], dtype=float)
    if topology.shape != (rows,) or momentum.shape != (rows,):
        raise ValueError("Proton branches must be aligned scalar arrays")
    if np.any(~np.isfinite(topology) | ~np.isin(topology, [0., 1., 2., 3., 4.])):
        raise ValueError("Invalid reco_topology: expected integer producer codes 0 through 4")
    tagged = (topology == 1.) | (topology == 2.)
    valid = np.where(tagged, np.isfinite(momentum) & (momentum > 0.), momentum == -999.)
    if not np.all(valid):
        raise ValueError("Proton topology/momentum contract mismatch: 1/2 require finite p>0; 0/3/4 require p=-999")
    return tagged


def proton_metadata(policy):
    return {"mode": policy, "implementation": "independent normalization auditor",
            "tag_definition": "reco_topology in {1,2} and finite reco_tagged_proton_p > 0",
            "producer_contract": "topology 0/3/4 have p=-999; topology 1/2 have finite p>0",
            "truth_denominator": "unchanged: all eligible truth rows, including untagged rows",
            "generated_normalization": "unchanged: all generated tuple rows"}


def _number(value):
    """Return a JSON-safe float; unavailable arithmetic is represented by null."""
    value = float(value)
    return value if np.isfinite(value) else None


def _ratio(numerator, denominator):
    return _number(numerator / denominator) if denominator > 0 else None


def _trapezoid(values, coordinates):
    integrate = getattr(np, "trapezoid", None)
    if integrate is None:
        integrate = np.trapz
    return float(integrate(values, coordinates))


def _positive_finite(name, value, allow_zero=False):
    value = float(value)
    if not np.isfinite(value) or (value < 0 if allow_zero else value <= 0):
        raise ValueError(name + " must be finite and " + ("nonnegative" if allow_zero else "positive"))
    return value


def validate_settings(settings):
    """Validate supplied numerical settings without changing their conventions."""
    names = (
        "electron_energy", "proton_energy", "q2_cut", "q2_max", "w_edges", "lumi_invfb",
        "br_ee", "mass", "mass_sigma", "mass_window_sigma", "true_mass_window",
        "sigma_epz_center", "sigma_epz_halfwidth", "min_eff",
    )
    missing = [name for name in names if name not in settings]
    if missing:
        raise ValueError("Missing physics settings: " + ", ".join(missing))
    result = dict(settings)
    result["proton_selection"] = proton_policy(result)
    result.setdefault("state", "1S")
    if result["state"] not in ("1S", "2S", "3S"):
        raise ValueError("state must be 1S, 2S or 3S")
    for name in ("electron_energy", "proton_energy", "q2_max", "lumi_invfb", "br_ee", "mass", "mass_sigma"):
        result[name] = _positive_finite(name, result[name])
    for name in ("mass_window_sigma", "true_mass_window", "sigma_epz_halfwidth", "min_eff"):
        result[name] = _positive_finite(name, result[name], allow_zero=True)
    if result["br_ee"] > 1:
        raise ValueError("br_ee cannot exceed one")
    if result["q2_cut"] is not None:
        result["q2_cut"] = _positive_finite("q2_cut", result["q2_cut"], allow_zero=True)
    result["sigma_epz_center"] = float(result["sigma_epz_center"])
    if not np.isfinite(result["sigma_epz_center"]):
        raise ValueError("sigma_epz_center must be finite")
    edges = np.asarray(result["w_edges"], dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or not np.all(np.isfinite(edges)) or not np.all(np.diff(edges) > 0):
        raise ValueError("w_edges must be finite and strictly increasing")
    if edges[0] <= M_P:
        raise ValueError("w_edges must start above the proton mass for a finite log-y photon-flux integral")
    result["w_edges"] = edges.tolist()
    return result


def photon_flux_dW(w, electron_energy, proton_energy, q2_max):
    """Analytically Q2-integrated electron EPA, including d y/d W.

    This is the exact expression in extract_gamma_p_cross_sections_bin_by_bin.py.
    Scalar and array W inputs are accepted; nonphysical y or Q2 ranges give zero.
    """
    values = np.asarray(w, dtype=float)
    denominator = 4.0 * electron_energy * proton_energy
    y = (values * values - M_P * M_P) / denominator
    output = np.zeros_like(values, dtype=float)
    valid = np.isfinite(values) & (y > 0.0) & (y < 1.0)
    if np.any(valid):
        yy = y[valid]
        qmin = M_E * M_E * yy * yy / (1.0 - yy)
        supported = (qmin > 0.0) & (qmin < q2_max)
        density = np.zeros_like(yy)
        yy_ok, qmin_ok = yy[supported], qmin[supported]
        density[supported] = (ALPHA_EM / (2.0 * np.pi)) * (
            ((1.0 + (1.0 - yy_ok) ** 2) / yy_ok) * np.log(q2_max / qmin_ok)
            - (2.0 * (1.0 - yy_ok) / yy_ok) * (1.0 - qmin_ok / q2_max)
        )
        output[valid] = density * (2.0 * values[valid] / denominator)
    return float(output) if output.ndim == 0 else output


def integrated_flux(w_low, w_high, electron_energy, proton_energy, q2_max, n_steps=500):
    """Framework-compatible 500-point trapezoid over W, including both edges."""
    grid = np.linspace(w_low, w_high, int(n_steps))
    return _trapezoid(photon_flux_dW(grid, electron_energy, proton_energy, q2_max), grid)


@lru_cache(maxsize=8)
def _gauss_rule(n):
    return np.polynomial.legendre.leggauss(int(n))


def _log_y_bounds(w_low, w_high, electron_energy, proton_energy, q2_max):
    denominator = 4.0 * electron_energy * proton_energy
    y_low = (w_low * w_low - M_P * M_P) / denominator
    y_high = (w_high * w_high - M_P * M_P) / denominator
    # Stable solution of m_e^2*y^2/(1-y) = Q2max.
    y_cap = 2.0 / (1.0 + np.sqrt(1.0 + 4.0 * M_E * M_E / q2_max))
    y_high = min(y_high, y_cap, np.nextafter(1.0, 0.0))
    if y_low <= 0:
        raise ValueError("Independent flux requires W_low above the proton mass")
    if y_high <= y_low:
        return None
    return np.log(y_low), np.log(y_high)


def independent_flux(w_low, w_high, electron_energy, proton_energy, q2_max, order=128):
    """Independent double Gauss quadrature of differential EPA in log y/log Q2.

    Integrates alpha/(2*pi) * [ (1+(1-y)^2)/(y*Q2)
    - 2*m_e^2*y/Q2^2 ] before Q2 integration. Thus this checks the analytic
    Q2 integral and the W-to-y Jacobian independently of photon_flux_dW().
    """
    bounds = _log_y_bounds(w_low, w_high, electron_energy, proton_energy, q2_max)
    if bounds is None:
        return 0.0
    nodes, weights = _gauss_rule(order)
    low, high = bounds
    u = 0.5 * (high + low) + 0.5 * (high - low) * nodes
    y = np.exp(u)
    qmin = M_E * M_E * y * y / (1.0 - y)
    qlog_low, qlog_high = np.log(qmin), np.log(q2_max)
    qhalf = 0.5 * (qlog_high - qlog_low)
    v = 0.5 * (qlog_high + qlog_low)[:, None] + qhalf[:, None] * nodes[None, :]
    q2 = np.exp(v)
    # The log-coordinate Jacobians cancel y*Q2 in the first term.
    density = ALPHA_EM / (2.0 * np.pi) * (
        (1.0 + (1.0 - y) ** 2)[:, None] - 2.0 * M_E * M_E * y[:, None] ** 2 / q2
    )
    q_integral = qhalf * np.sum(density * weights[None, :], axis=1)
    return float(0.5 * (high - low) * np.sum(q_integral * weights))


def estarlight_real_photon_pb(w):
    """Current-source Upsilon(1S) real-photon benchmark in pb, not a campaign fit.

    photonNucleusCrossSection.cpp lines1048--1055 gives fm^2 as
    1e-10 * 6.4 * W^0.74 * (1-(m_Upsilon+m_p)^2/W^2)^2.
    One fm^2 is 1e10 pb. This function sets the physical threshold explicitly.
    It omits the finite-Q2 production response and does not validate generator
    revision, card cuts, or the supplied total ep normalization cross section.
    """
    values = np.asarray(w, dtype=float)
    output = np.zeros_like(values)
    threshold = MODEL_UPSILON_MASS + M_P
    valid = np.isfinite(values) & (values > threshold)
    ww = values[valid]
    output[valid] = 6.4 * np.power(ww, 0.74) * np.square(1.0 - threshold * threshold / (ww * ww))
    return float(output) if output.ndim == 0 else output


def model_weighted_integral(w_low, w_high, electron_energy, proton_energy, q2_max, order=256):
    """Return integral flux*sigma(Q2=0) in pb and its flux-weighted average.

    The numerator and denominator use the same 256-point Gauss rule in log y.
    A threshold-straddling interval is split before integrating the numerator.
    """
    bounds = _log_y_bounds(w_low, w_high, electron_energy, proton_energy, q2_max)
    if bounds is None:
        return 0.0, None
    nodes, weights = _gauss_rule(order)
    denominator = 4.0 * electron_energy * proton_energy

    def integrate(low, high, with_model):
        if high <= low:
            return 0.0
        u = 0.5 * (high + low) + 0.5 * (high - low) * nodes
        y = np.exp(u)
        w = np.sqrt(M_P * M_P + denominator * y)
        # dW/du = (4 Ee Ep)*y/(2 W).
        density = photon_flux_dW(w, electron_energy, proton_energy, q2_max) * denominator * y / (2.0 * w)
        if with_model:
            density *= estarlight_real_photon_pb(w)
        return float(0.5 * (high - low) * np.sum(weights * density))

    low, high = bounds
    flux = integrate(low, high, False)
    threshold_y = ((MODEL_UPSILON_MASS + M_P) ** 2 - M_P * M_P) / denominator
    threshold_log_y = np.log(threshold_y)
    numerator = integrate(max(low, threshold_log_y), high, True)
    return numerator, _ratio(numerator, flux)


class PhysicsAccumulator:
    """Accumulate exact central-value W selections over any number of chunks.

    Required branches are never optional, including true_Q2. No weight_proj is
    consumed. The reco lane does not require a truth-mass match, matching the
    framework. All cuts and the final histogram edge are inclusive.

    Nonfinite weights are counted and excluded from weighted accumulations only;
    raw counts retain their rows. This permits a diagnostic report, but the
    framework rejects such inputs. Consumers must treat those diagnostics as a
    failure and not interpret finite-only subtotals as a repaired sample.
    """
    def __init__(self, settings):
        self.settings = validate_settings(settings)
        self.required_branches = required_branches(self.settings)
        self.edges = np.asarray(self.settings["w_edges"], dtype=float)
        n = len(self.edges) - 1
        self.histograms = {name: np.zeros(n, dtype=np.int64 if name.endswith("count") else float)
                           for name in ("truth_count", "reco_count", "truth_sumw", "reco_sumw", "truth_sumw2", "reco_sumw2")}
        self.upper_edge_counts = {name: np.zeros(n, dtype=np.int64) for name in ("truth", "reco")}
        self.cutflow = {stage: {"count": 0, "sumw": 0.0, "sumw2": 0.0, "nonfinite_weight_count": 0}
                        for stage in CUTFLOW_STAGES}
        self.diagnostics = {"chunks": 0, "rows": 0, "nonfinite_weight_rows": 0,
                            "negative_weight_rows": 0, "zero_weight_rows": 0,
                            "nonfinite_by_branch": {name: 0 for name in self.required_branches},
                            "W_domain_counts": {
                                name: {key: 0 for key in ("sentinel_minus999", "negative_non_sentinel",
                                                         "nonnegative_below_upsilon_proton_threshold",
                                                         "above_approximate_ep_endpoint")}
                                for name in ("true_W", "reco_W")}}

    def update(self, arrays):
        missing = [name for name in self.required_branches if name not in arrays]
        if missing:
            raise ValueError("Missing required physics branches: " + ", ".join(missing))
        data = {name: np.asarray(arrays[name], dtype=float) for name in self.required_branches}
        nrows = len(data["weight"]) if data["weight"].ndim == 1 else -1
        if nrows < 0 or any(value.ndim != 1 or len(value) != nrows for value in data.values()):
            raise ValueError("Physics branches must be aligned scalar arrays")
        weights = data["weight"]
        finite_weight = np.isfinite(weights)
        finite_squared = np.square(weights[finite_weight])
        if not np.all(np.isfinite(finite_squared)):
            raise ValueError("Finite input weights overflow sumw2; refusing invalid accounting")
        s = self.settings
        q2 = np.ones(nrows, dtype=bool) if s["q2_cut"] is None else ((data["true_Q2"] >= 0) & (data["true_Q2"] <= s["q2_cut"]))
        truthmass = (data["true_mass"] >= s["mass"] - s["true_mass_window"]) & (data["true_mass"] <= s["mass"] + s["true_mass_window"])
        truthmassq2 = truthmass & q2
        truth = truthmassq2 & (data["true_W"] >= self.edges[0]) & (data["true_W"] <= self.edges[-1])
        mass_window = s["mass_window_sigma"] * s["mass_sigma"]
        reco_mass = (data["reco_mass"] >= s["mass"] - mass_window) & (data["reco_mass"] <= s["mass"] + mass_window)
        reco_epz = reco_mass & (np.abs(data["sigma_epz"] - s["sigma_epz_center"]) <= s["sigma_epz_halfwidth"])
        reco_q2 = reco_epz & q2
        reco_before_proton = reco_q2 & (data["reco_W"] >= self.edges[0]) & (data["reco_W"] <= self.edges[-1])
        tagged = independent_proton_mask(data, s["proton_selection"])
        reco = reco_before_proton & tagged
        masks = dict(zip(CUTFLOW_STAGES, (np.ones(nrows, dtype=bool), truthmass, truthmassq2, truth,
                                        reco_mass, reco_epz, reco_q2, reco_before_proton, reco)))
        for stage, mask in masks.items():
            target = self.cutflow[stage]
            selected_weight = weights[mask & finite_weight]
            target["count"] += int(mask.sum())
            target["sumw"] += float(np.sum(selected_weight, dtype=float))
            target["sumw2"] += float(np.sum(selected_weight * selected_weight, dtype=float))
            target["nonfinite_weight_count"] += int(np.sum(mask & ~finite_weight))
        for prefix, values, mask in (("truth", data["true_W"], truth), ("reco", data["reco_W"], reco)):
            # Match statistical_components._indices and bincount, including
            # values equal to the final upper edge in the final bin.
            n_bins = len(self.edges) - 1
            indices = np.searchsorted(self.edges, values, side="right") - 1
            indices[values == self.edges[-1]] = n_bins - 1
            self.histograms[prefix + "_count"] += np.bincount(indices[mask], minlength=n_bins)
            valid = mask & finite_weight
            self.histograms[prefix + "_sumw"] += np.bincount(indices[valid], weights=weights[valid], minlength=n_bins)
            self.histograms[prefix + "_sumw2"] += np.bincount(indices[valid], weights=weights[valid] ** 2, minlength=n_bins)
            # These counts use all corresponding non-bin cuts. An event on
            # an internal upper edge belongs to the next bin but is recorded
            # here to diagnose cross-report differences in final-bin roles.
            for i, upper in enumerate(self.edges[1:]):
                self.upper_edge_counts[prefix][i] += int(np.sum(mask & (values == upper)))
        self.diagnostics["chunks"] += 1
        self.diagnostics["rows"] += nrows
        self.diagnostics["nonfinite_weight_rows"] += int(np.sum(~finite_weight))
        self.diagnostics["negative_weight_rows"] += int(np.sum(finite_weight & (weights < 0)))
        self.diagnostics["zero_weight_rows"] += int(np.sum(finite_weight & (weights == 0)))
        for name, values in data.items():
            self.diagnostics["nonfinite_by_branch"][name] += int(np.sum(~np.isfinite(values)))
        threshold = s["mass"] + M_P
        endpoint = np.sqrt(M_P * M_P + 4.0 * s["electron_energy"] * s["proton_energy"])
        for name in ("true_W", "reco_W"):
            values = data[name]
            finite = np.isfinite(values)
            counts = self.diagnostics["W_domain_counts"][name]
            counts["sentinel_minus999"] += int(np.sum(values == -999.0))
            counts["negative_non_sentinel"] += int(np.sum(finite & (values < 0) & (values != -999.0)))
            counts["nonnegative_below_upsilon_proton_threshold"] += int(np.sum(finite & (values >= 0) & (values < threshold)))
            counts["above_approximate_ep_endpoint"] += int(np.sum(finite & (values > endpoint)))
        return self

    def result(self):
        s = self.settings
        model_applicable = s["state"] == "1S" and np.isclose(s["mass"], MODEL_UPSILON_MASS, rtol=0., atol=1e-8)
        all_rows = self.cutflow["all"]["count"]
        all_weight = self.cutflow["all"]["sumw"]
        rows = []
        for i, (low, high) in enumerate(zip(self.edges[:-1], self.edges[1:])):
            flux = integrated_flux(low, high, s["electron_energy"], s["proton_energy"], s["q2_max"])
            check = independent_flux(low, high, s["electron_energy"], s["proton_energy"], s["q2_max"])
            model_integral, model_average = (model_weighted_integral(low, high, s["electron_energy"], s["proton_energy"], s["q2_max"])
                                             if model_applicable else (None, None))
            truth = float(self.histograms["truth_sumw"][i])
            reco = float(self.histograms["reco_sumw"][i])
            efficiency = reco / truth if truth > 0 else np.nan
            valid = truth > 0 and reco > 0 and np.isfinite(efficiency) and efficiency > 0
            corrected = reco / efficiency if valid else np.nan
            denominator = 1000.0 * s["lumi_invfb"] * s["br_ee"] * flux
            factor = 1.0 / denominator if denominator > 0 else np.nan
            sigma = corrected * factor
            row = {"W_low": float(low), "W_high": float(high), "W_center": float((low + high) / 2.0)}
            row.update({key: int(value[i]) if key.endswith("count") else _number(value[i]) for key, value in self.histograms.items()})
            row.update({"eff_reco_over_true": _number(efficiency), "corrected_visible": _number(corrected),
                        "upper_edge_inclusive": i == len(self.edges) - 2,
                        "truth_count_exact_upper_edge": int(self.upper_edge_counts["truth"][i]),
                        "reco_count_exact_upper_edge": int(self.upper_edge_counts["reco"][i]),
                        "sigma_gamma_p_pb": _number(sigma), "sigma_gamma_p_truth_closure_pb": _number(truth * factor),
                        "closure_unfold_over_truth": _ratio(corrected, truth),
                        "flux_int": flux, "independent_flux": check,
                        "flux_relative_difference": _ratio(flux - check, check),
                        "model_average_pb": model_average, "model_center_pb": estarlight_real_photon_pb((low + high) / 2.0) if model_applicable else None,
                        "model_weighted_integral_pb": model_integral,
                        "truth_fraction_all_rows": _ratio(self.histograms["truth_count"][i], all_rows),
                        "truth_weight_fraction_all_rows": _ratio(truth, all_weight),
                        "reco_fraction_all_rows": _ratio(self.histograms["reco_count"][i], all_rows),
                        "plotted": bool(np.isfinite(sigma) and sigma > 0 and efficiency > s["min_eff"])})
            rows.append(row)
        cutflow = {stage: {**values, "fraction_of_all_rows": _ratio(values["count"], all_rows),
                           "weight_fraction_of_all_rows": _ratio(values["sumw"], all_weight)}
                   for stage, values in self.cutflow.items()}
        diagnostics = {**self.diagnostics,
                       "W_domain_definition": {
                           "scope": "All tuple rows before mass/Q2/W/reco cuts; diagnostics do not change selections",
                           "upsilon_proton_threshold_GeV": s["mass"] + M_P,
                           "approximate_ep_endpoint_GeV": float(np.sqrt(M_P * M_P + 4.0 * s["electron_energy"] * s["proton_energy"])),
                           "endpoint_convention": "sqrt(m_p^2 + 4 Ee Ep), as in the framework EPA; finite-mass producer kinematics can differ slightly",
                           "reco_interpretation": "Reconstruction can migrate beyond physical truth support; these counts are diagnostic, not proof of corrupt events"},
                       "framework_would_reject_nonfinite_weights": bool(self.diagnostics["nonfinite_weight_rows"]),
                       "weighted_values_are_finite_weight_subtotals": bool(self.diagnostics["nonfinite_weight_rows"]),
                       "correctable_bins": sum(row["corrected_visible"] is not None for row in rows),
                       "plotted_bins": sum(row["plotted"] for row in rows)}
        return {"bins": rows, "cutflow": cutflow, "diagnostics": diagnostics,
                "settings": {name: s[name] for name in ("state", "electron_energy", "proton_energy", "q2_cut", "q2_max", "w_edges", "lumi_invfb", "br_ee", "mass", "mass_sigma", "mass_window_sigma", "true_mass_window", "sigma_epz_center", "sigma_epz_halfwidth", "min_eff", "proton_selection")},
                "proton_selection": proton_metadata(s["proton_selection"]),
                "conventions": {
                    "selection": "Independent truth and reco lanes; declared proton category applies only to reco; inclusive mass, Q2 and W boundaries; final W bin includes its upper edge",
                    "weights": "Stored weight only; no weight_proj, sample renormalization, luminosity rescaling or branching-factor inference",
                    "correction": "R/(R/T)=T is a same-MC algebraic closure, not an independent efficiency validation; R/T may exceed one with migrations",
                    "units": "sigma_gamma_p_pb = corrected_visible/(1000*lumi_invfb*br_ee*flux_int)",
                    "finite_handling": "Raw counts include nonfinite-weight rows; weighted subtotals exclude those weights and are invalid for certification when diagnostics flag them",
                    "flux": "Framework 500-point trapezoid in W; independent 128x128 Gauss integration in log(y), log(Q2)",
                    "model": "6.4*W^0.74*(1-(m_Upsilon+m_p)^2/W^2)^2 pb above threshold, Q2=0 Upsilon(1S) current-source benchmark",
                    "model_status": "AVAILABLE" if model_applicable else "NOT_APPLICABLE",
                    "model_source": MODEL_SOURCE,
                    "model_scope": "Not a full finite-Q2 generator calculation; exact campaign revision and total-rate provenance are not certified",
                }}


def physics_summary(arrays, settings):
    """Convenience wrapper for a single in-memory sample."""
    return PhysicsAccumulator(settings).update(arrays).result()
