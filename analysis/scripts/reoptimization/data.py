"""Read-only, generated-denominator-preserving inputs for reoptimization.

Only prepared pure-state tuples are accepted. ``weight`` is the visible
dielectron expected yield at the preparation luminosity. No selection enters
normalization, and invalid reconstructed values are inefficiencies rather than
reasons to discard otherwise valid generated truth rows.
"""

import hashlib
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata

import numpy as np
from scipy.integrate import quad


STATE_DEFAULTS = {
    "1S": {"mass": 9.46030, "sigma": 0.053, "br_ee": 0.0238,
           "visible_xsec_fb": 210.136, "q2_max": 1.0},
    "2S": {"mass": 10.02326, "sigma": 0.056, "br_ee": 0.0191,
           "visible_xsec_fb": 53.133, "q2_max": 0.01},
    "3S": {"mass": 10.35520, "sigma": 0.059, "br_ee": 0.0218,
           "visible_xsec_fb": 42.969, "q2_max": 0.01},
}
ALPHA_EM = 1.0 / 137.035999084
M_E = 0.00051099895
M_P = 0.9382720813
NORMALIZATION_RTOL = 1e-7
NORMALIZATION_ATOL = 1e-10
REQUIRED_BRANCHES = (
    "weight", "true_mass", "true_W", "true_Q2", "true_pT", "true_t",
    "reco_mass", "reco_W", "reco_pT", "sigma_epz",
)
KNOWN_PROFILES = {
    "jul26-estarlight-939b11a", "jul26-purestates-estarlight-939b11a",
}
TRUTH_SCHEMA_POLICY = "dvmp_coherent_minus999_v1"
TRUTH_BLOCK = ("true_mass", "true_W", "true_Q2", "true_pT", "true_t")


def _available_truth(arrays):
    """Recognize the producer's all-event missing-truth row, without repair.

    DVMP_Upsilon_Analysis(3).C resets this entire block to -999 each event,
    fills it only inside if (foundGenEvent), and writes every event once.
    A true_t-only invalid value remains an observable-specific issue. Other
    partial blocks and unrecognized negative mass/W/pT values are not accepted.
    """
    missing = np.logical_and.reduce([arrays[name] == -999.0 for name in TRUTH_BLOCK])
    main_names = ("true_mass", "true_W", "true_Q2", "true_pT")
    partial = np.logical_or.reduce([arrays[name] == -999.0 for name in main_names]) & ~missing
    if np.any(partial):
        examples = np.flatnonzero(partial)[:5].tolist()
        raise ValueError("Inconsistent partial -999 generated truth block: %d rows, examples %s. "
                         "Expected all of %s to be -999 for a missing generated combination; "
                         "audit producer/schema rather than repairing values."
                         % (int(partial.sum()), examples, ", ".join(TRUTH_BLOCK)))
    for name in main_names:
        if not np.all(np.isfinite(arrays[name])):
            raise ValueError("Nonfinite generated truth in %s. Repair/audit the input rather than "
                             "silently dropping truth denominator rows." % name)
    for name in ("true_mass", "true_W", "true_pT"):
        invalid = (arrays[name] < 0) & ~missing
        if np.any(invalid):
            values = np.unique(arrays[name][invalid])[:5].tolist()
            raise ValueError("Negative generated truth in %s outside the coherent -999 producer "
                             "sentinel: %d rows, example values %s. An explicit producer/schema "
                             "audit is required; no absolute-value repair is applied."
                             % (name, int(invalid.sum()), values))
    if np.all(missing):
        raise ValueError("No available generated truth: every row has the coherent -999 producer sentinel")
    return ~missing


def _number(value, name, positive=False):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(name + " must be a finite number, not a boolean")
    try:
        result = float(value)
    except (ValueError, TypeError):
        raise ValueError(name + " must be a finite number") from None
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValueError(name + " must be finite" + (" and positive" if positive else ""))
    return result


def state_settings(state_config):
    """Resolve explicit state settings without overriding deliberate values."""
    state = state_config.get("state")
    if state not in STATE_DEFAULTS:
        raise ValueError("state_config.state must be one of 1S, 2S, 3S")
    result = dict(STATE_DEFAULTS[state])
    result.update(state_config)
    for key in ("mass", "sigma", "br_ee", "visible_xsec_fb", "q2_max"):
        result[key] = _number(result[key], "state_config." + key, positive=True)
    if result["br_ee"] >= 1:
        raise ValueError("state_config.br_ee must lie strictly between 0 and 1")
    return result


def _settings(config, state_config):
    state = state_settings(state_config)
    energy_e = _number(config.get("electron_energy", 9.0), "electron_energy", True)
    energy_p = _number(config.get("proton_energy", 275.0), "proton_energy", True)
    target_lumi = _number(config.get("target_lumi_invfb", 2.5), "target_lumi_invfb", True)
    w_range = config.get("w_range", [20.0, 100.0])
    if not isinstance(w_range, (list, tuple)) or len(w_range) != 2:
        raise ValueError("w_range must contain exactly two bounds")
    w_min, w_max = [_number(x, "w_range", True) for x in w_range]
    if w_min >= w_max:
        raise ValueError("w_range bounds must increase")
    if w_max <= M_P or w_min >= math.sqrt(4 * energy_e * energy_p + M_P * M_P):
        raise ValueError("w_range contains no physical photon-flux support")
    return state, energy_e, energy_p, target_lumi, w_min, w_max


def photon_flux_dW(w, electron_energy, proton_energy, q2_max):
    """Equivalent-photon dN/dW used by the existing differential extraction."""
    denominator = 4.0 * electron_energy * proton_energy
    y = (w * w - M_P * M_P) / denominator
    if y <= 0.0 or y >= 1.0:
        return 0.0
    q2_min = M_E * M_E * y * y / (1.0 - y)
    if q2_min <= 0.0 or q2_min >= q2_max:
        return 0.0
    dndy = (ALPHA_EM / (2.0 * math.pi)) * (
        ((1.0 + (1.0 - y) ** 2) / y) * math.log(q2_max / q2_min)
        - 2.0 * (1.0 - y) / y * (1.0 - q2_min / q2_max)
    )
    return max(0.0, dndy * 2.0 * w / denominator)


def integrated_flux(config, state_config):
    """Integrate the state-specific analytic flux, clipping physical support.

    Adaptive quadrature is independent of the legacy fixed-grid trapezoid.
    This is a numerical integration cross-check, not a generator validation.
    """
    state, ee, ep, _, lower, upper = _settings(config, state_config)
    qmax = state["q2_max"]
    # Positive root of m_e^2*y^2 + Qmax^2*y - Qmax^2 = 0,
    # in a rationalized form that avoids catastrophic cancellation.
    y_max = 2.0 / (1.0 + math.sqrt(1.0 + 4.0 * M_E * M_E / qmax))
    physical_upper = math.sqrt(M_P * M_P + 4.0 * ee * ep * y_max)
    upper = min(upper, physical_upper)
    lower = max(lower, M_P)
    if upper <= lower:
        raise ValueError("Selection has no physical photon-flux support")
    # The lower W boundary must be above the proton mass for a finite integral.
    if lower <= M_P:
        raise ValueError("Photon-flux integration needs W_min above the proton mass")
    value, error = quad(photon_flux_dW, lower, upper,
                        args=(ee, ep, qmax), epsabs=1e-12, epsrel=1e-10, limit=200)
    if not math.isfinite(value) or value <= 0 or error > max(1e-9, 1e-6 * value):
        raise ValueError("Photon-flux quadrature failed its numerical accuracy check")
    return float(value)


def cross_section_factor(config, state_config):
    """Convert visible truth yield to inclusive flux-averaged sigma in pb.

    Divide additionally by the variable-bin width for a differential spectrum.
    """
    state, _, _, luminosity, _, _ = _settings(config, state_config)
    return 1.0 / (1000.0 * luminosity * state["br_ee"]
                  * integrated_flux(config, state_config))


def observable_cross_section_scaling(config, state_config, observable, edges):
    """Return explicit per-bin count-to-cross-section factors in pb units.

    W is a photon-flux-weighted gamma-p cross section averaged within each W
    interval, so no additional division by Delta W is made. Rapidity is the
    Upsilon rapidity, not DIS inelasticity; its output is inclusive ep dσ/dy
    within the configured W/Q2 region, without an arbitrary flux division.
    """
    edges = np.asarray(edges, dtype=float)
    if (edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all()
            or np.any(np.diff(edges) <= 0)):
        raise ValueError("Cross-section edges must be finite and increasing")
    state, _, _, luminosity, w_min, w_max = _settings(config, state_config)
    widths = np.diff(edges)
    base = 1.0 / (1000.0 * luminosity * state["br_ee"])
    result = dict(observable=observable, bin_edges=edges.tolist(),
                  luminosity_invfb=luminosity, branching_fraction_ee=state["br_ee"],
                  conversion_from_fb_to_pb=0.001,
                  normalization="inclusive Upsilon; input weights are visible dielectron yields")
    if observable == "w":
        if edges[0] < w_min - 1e-12 or edges[-1] > w_max + 1e-12:
            raise ValueError("W extraction bins must lie within the configured truth W range")
        fluxes = [integrated_flux(dict(config, w_range=[float(a), float(b)]), state)
                  for a, b in zip(edges[:-1], edges[1:])]
        result.update(quantity="sigma_gamma_p_w", quantity_label="Flux-averaged gamma-p cross section",
                      units="pb", observable_units="GeV", bin_flux=fluxes,
                      factors=(base / np.asarray(fluxes)).tolist(),
                      formula="N_i / (1000 * L * B_ee * integral_bin(dN_gamma/dW dW)); no Delta W division",
                      photon_flux_applied=True, bin_width_applied=False)
    elif observable == "y":
        flux = integrated_flux(config, state)
        result.update(quantity="dsigma_ep_dy", quantity_label="Inclusive ep differential cross section",
                      units="pb", observable_units="dimensionless", bin_flux=None,
                      factors=(base / widths).tolist(),
                      formula="N_i / (1000 * L * B_ee * Delta y_i); no photon-flux division",
                      photon_flux_applied=False, bin_width_applied=True,
                      alternate_factors=(base / widths / flux).tolist(),
                      alternate_quantity="flux_normalized_ep_rapidity",
                      alternate_label="(1/Phi) d sigma_ep / d y",
                      alternate_units="pb", alternate_flux_integral=flux,
                      alternate_interpretation="Legacy flux-normalized ep rapidity spectrum; not a universal gamma-p rapidity distribution")
    elif observable in ("pt2", "abs_t_via_pt2"):
        flux = integrated_flux(config, state)
        result.update(quantity="dsigma_gamma_p_dpt2" if observable == "pt2" else "dsigma_gamma_p_dabs_t",
                      quantity_label="Flux-averaged gamma-p differential cross section",
                      units="pb/GeV^2", observable_units="GeV^2", flux_integral=flux,
                      bin_flux=None, factors=(base / flux / widths).tolist(),
                      formula="N_i / (1000 * L * B_ee * integrated_flux * Delta x_i)",
                      photon_flux_applied=True, bin_width_applied=True)
    else:
        raise ValueError("Unknown cross-section observable: " + str(observable))
    result["x_units"] = result["observable_units"]
    result["ylabel"] = {"w": r"$\langle\sigma_{\gamma p}\rangle_{\Phi}$ [pb]",
                        "y": r"$d\sigma_{ep}/dy_{\Upsilon}$ [pb]",
                        "pt2": r"$d\sigma_{\gamma p}/dp_T^2$ [pb/GeV$^2$]",
                        "abs_t_via_pt2": r"$d\sigma_{\gamma p}/d|t|$ [pb/GeV$^2$]"}[observable]
    if observable == "y":
        result["alternate_ylabel"] = r"$(1/\Phi)\,d\sigma_{ep}/dy_{\Upsilon}$ [pb]"
    return result


def _same_number(metadata, key, expected, exact_count=False):
    if key not in metadata:
        raise ValueError("Prepared ROOT metadata is missing " + key)
    actual = _number(metadata[key], "framework_normalization." + key)
    if not math.isclose(actual, expected,
                        rel_tol=0.0 if exact_count else NORMALIZATION_RTOL,
                        abs_tol=1e-6 if exact_count else NORMALIZATION_ATOL):
        raise ValueError("Prepared ROOT metadata %s=%s disagrees with expected %s"
                         % (key, actual, expected))


def _sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_sample(root_path, state_config, config):
    """Load a normalized state tuple without modifying files or normalization.

    Normalization tolerance is relative 1e-7 plus absolute 1e-10 for weights
    and totals, and absolute 1e-6 for generated counts. All returned arrays have
    the original number and ordering of rows. The deterministic split is made
    before all analysis cuts and is later exposure-corrected by its configured
    probability, never by a selected-event normalization.
    """
    state, ee, ep, target_lumi, w_min, w_max = _settings(config, state_config)
    source_path = Path(root_path).expanduser().resolve()
    if not source_path.is_file():
        raise FileNotFoundError(
            "Prepared ROOT input not found: %s. Set --framework-root or the state run/input "
            "path to the existing allstates-projected-L10 run containing "
            "artifacts/inputs/9x275_normalized.root; do not substitute a raw merged tuple."
            % source_path)
    try:
        import uproot
    except ImportError:
        raise ImportError("Reading the prepared tuple requires uproot>=5. Install the "
                          "reoptimization requirements in the analysis Python environment.") from None

    tree_name = config.get("tree", "UpsilonTree")
    if not isinstance(tree_name, str) or not tree_name:
        raise ValueError("tree must be a nonempty string")
    requested_observables = config.get("requested_observables", ("pt2", "abs_t_via_pt2"))
    if not isinstance(requested_observables, (list, tuple)) or any(
            x not in ("pt2", "abs_t_via_pt2", "w", "y") for x in requested_observables):
        raise ValueError("requested_observables must be a list of supported observable names")
    unavailable_observables = {}
    proton_policy = proton_selection_metadata()
    if "proton_selection" in config and config["proton_selection"] != proton_policy:
        raise ValueError("Configured proton selection differs from the active selection policy")
    required_branches = REQUIRED_BRANCHES + tuple(proton_branches())
    train_fraction = _number(config.get("train_fraction", 0.5), "train_fraction")
    if not 0.0 < train_fraction < 1.0:
        raise ValueError("train_fraction must lie strictly between 0 and 1")
    seed = config.get("seed", 12345)
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a nonnegative integer")
    truth_window = _number(config.get("truth_mass_window", 0.2), "truth_mass_window", True)
    reco_width = _number(config.get("reco_mass_window_sigma", 3.0), "reco_mass_window_sigma", True)
    epz_center = _number(config.get("sigma_epz_center", 2.0 * ee), "sigma_epz_center")
    epz_halfwidth = _number(config.get("sigma_epz_halfwidth", 0.2 * ee), "sigma_epz_halfwidth", True)
    min_state_fraction = _number(config.get("min_state_fraction", 0.95), "min_state_fraction")
    min_available_fraction = _number(config.get("min_available_truth_fraction", 0.95),
                                     "min_available_truth_fraction")
    max_q2_outliers = _number(config.get("max_q2_outlier_fraction", 0.001), "max_q2_outlier_fraction")
    if (not 0.0 < min_state_fraction <= 1 or not 0 <= max_q2_outliers < 1
            or not 0.0 < min_available_fraction <= 1):
        raise ValueError("Invalid state-population/Q2-outlier thresholds")

    with uproot.open(source_path) as source:
        if tree_name not in source:
            raise ValueError("Prepared ROOT input has no tree " + tree_name)
        tree = source[tree_name]
        count = int(tree.num_entries)
        if count <= 1:
            raise ValueError("At least two generated rows are required for independent splits")
        missing = sorted(set(required_branches) - set(tree.keys()))
        if missing:
            raise ValueError("Prepared ROOT input is missing required branches: " + ", ".join(missing))
        if "y" in requested_observables:
            missing_y = sorted({"true_y", "reco_y"} - set(tree.keys()))
            if missing_y:
                unavailable_observables["y"] = "Missing rapidity branches: " + ", ".join(missing_y)
            else:
                required_branches += ("true_y", "reco_y")
        if "framework_normalization" not in source:
            raise ValueError("Input has no framework_normalization metadata. Use the separately "
                             "prepared normalized ROOT, not the raw merged file.")
        try:
            normalization = json.loads(str(source["framework_normalization"]))
        except (TypeError, ValueError):
            raise ValueError("framework_normalization is not readable JSON") from None
        if not isinstance(normalization, dict):
            raise ValueError("framework_normalization must contain a JSON object")
        profile = normalization.get("normalization_profile")
        expected_profile = state.get("normalization_profile", config.get("normalization_profile"))
        if expected_profile is not None:
            if not isinstance(expected_profile, str) or not expected_profile or profile != expected_profile:
                raise ValueError("Prepared normalization profile does not match configured normalization_profile")
        elif profile not in KNOWN_PROFILES:
            raise ValueError("Unknown normalization profile. Set an explicit matching normalization_profile "
                             "only for an intentionally supplied alternative normalization.")
        if profile == "jul26-estarlight-939b11a" and state["state"] != "1S":
            raise ValueError("The jul26-estarlight-939b11a profile is 1S-only; use pure-state prepared inputs")
        if normalization.get("normalization_scope") != (
                "All generated events, before truth-state, Q2, W, and reconstruction cuts"):
            raise ValueError("Prepared normalization_scope is not the complete generated sample")
        source_lumi = _number(normalization.get("luminosity_invfb"), "source luminosity", True)
        expected_lumi = state.get("expected_source_lumi_invfb", config.get("expected_source_lumi_invfb"))
        if expected_lumi is not None and not math.isclose(source_lumi, _number(expected_lumi, "expected_source_lumi_invfb", True), rel_tol=1e-10):
            raise ValueError("Prepared source luminosity disagrees with expected_source_lumi_invfb")
        total = state["visible_xsec_fb"] * source_lumi
        per_event = total / count
        for key, expected, integer in (
            ("entries", count, True), ("generated_events_h_EventCount", count, True),
            ("cross_section_visible_ee_fb", state["visible_xsec_fb"], False),
            ("generator_branching_fraction_ee", state["br_ee"], False),
            ("target_sum_weight", total, False), ("corrected_event_weight", per_event, False),
        ):
            _same_number(normalization, key, expected, integer)
        if "written_sum_weight" in normalization:
            _same_number(normalization, "written_sum_weight", total)
        if "h_EventCount" not in source:
            raise ValueError("Prepared ROOT input is missing h_EventCount generated denominator")
        histogram = source["h_EventCount"]
        try:
            ordinary = np.asarray(histogram.values(flow=False), dtype=float)
            flow = np.asarray(histogram.values(flow=True), dtype=float)
        except (AttributeError, TypeError, ValueError):
            raise ValueError("h_EventCount is not a readable one-dimensional count histogram") from None
        if ordinary.ndim != 1 or flow.ndim != 1 or flow.size != ordinary.size + 2:
            raise ValueError("h_EventCount has unrecognized dimension/flow geometry")
        if (not np.all(np.isfinite(flow)) or np.any(flow < 0)
                or np.any(flow[[0, -1]] != 0)
                or not math.isclose(float(ordinary.sum()), count, rel_tol=0, abs_tol=1e-6)):
            raise ValueError("h_EventCount must equal all tuple rows with finite nonnegative contents and zero flows")
        arrays = {}
        for name in required_branches:
            try:
                values = np.asarray(tree[name].array(library="np"), dtype=np.float64)
            except (ValueError, TypeError):
                if name in ("true_y", "reco_y"):
                    unavailable_observables["y"] = "Rapidity branch %s is not a scalar numeric array" % name
                    continue
                raise ValueError("Branch %s is not a scalar numeric array" % name) from None
            if values.ndim != 1 or values.shape != (count,):
                if name in ("true_y", "reco_y"):
                    unavailable_observables["y"] = "Rapidity branch %s does not have one scalar per generated row" % name
                    continue
                raise ValueError("Branch %s does not have one scalar per generated row" % name)
            arrays[name] = values

    weight = arrays["weight"]
    if not np.all(np.isfinite(weight) & (weight > 0)):
        raise ValueError("Ordinary weight must be finite and positive for every generated row")
    if not np.allclose(weight, per_event, rtol=NORMALIZATION_RTOL, atol=NORMALIZATION_ATOL):
        raise ValueError("Ordinary weights are not uniform at the full-generated expected normalization")
    if not math.isclose(float(np.sum(weight, dtype=np.float64)), total,
                        rel_tol=NORMALIZATION_RTOL, abs_tol=NORMALIZATION_ATOL):
        raise ValueError("Full-generated sum(weight) disagrees with visible cross section times source luminosity")

    available_truth = _available_truth(arrays)
    missing_truth = ~available_truth
    available_count = int(available_truth.sum())
    available_fraction = float(np.mean(available_truth))
    if available_fraction < min_available_fraction:
        raise ValueError("The available generated truth fraction %.6g is below min_available_truth_fraction=%.6g "
                         "(%d coherent -999 rows of %d). Audit producer truth finding before changing "
                         "this explicit input-quality threshold; generated normalization must remain unchanged."
                         % (available_fraction, min_available_fraction, int(missing_truth.sum()), count))
    state_mask = available_truth & (np.abs(arrays["true_mass"] - state["mass"]) <= truth_window)
    state_fraction = float(state_mask.sum() / available_count)
    if state_fraction < min_state_fraction:
        raise ValueError("Tuple mass population is inconsistent with configured pure state " + state["state"])
    q2_outlier = available_truth & ((arrays["true_Q2"] < 0) | (arrays["true_Q2"] > state["q2_max"] * (1 + 1e-6)))
    q2_outlier_fraction = float(q2_outlier.sum() / available_count)
    if q2_outlier_fraction > max_q2_outliers:
        raise ValueError("Truth Q2 population disagrees with the configured generated range/rate. "
                         "Do not assign a new Q2 acceptance without the corresponding input normalization.")
    truth_selected = (state_mask & (arrays["true_W"] >= w_min) & (arrays["true_W"] <= w_max)
                      & (arrays["true_Q2"] >= 0) & (arrays["true_Q2"] <= state["q2_max"]))
    reco_finite = np.ones(count, dtype=bool)
    for name in ("reco_mass", "reco_W", "reco_pT", "sigma_epz"):
        reco_finite &= np.isfinite(arrays[name])
    reco_valid = reco_finite & (arrays["reco_pT"] >= 0)
    reco_only_selected = (reco_valid
                          & (np.abs(arrays["reco_mass"] - state["mass"]) <= reco_width * state["sigma"])
                          & (np.abs(arrays["sigma_epz"] - epz_center) <= epz_halfwidth)
                          & (arrays["reco_W"] >= w_min) & (arrays["reco_W"] <= w_max))
    reco_only_before_proton = reco_only_selected.copy()
    proton_selected = proton_tag_mask(arrays)
    reco_only_selected &= proton_selected
    reco_selected = truth_selected & reco_only_selected
    if not np.any(truth_selected):
        raise ValueError("No truth events survive the configured fiducial selection")
    split_seed = int(seed) + list(STATE_DEFAULTS).index(state["state"]) * 1000003
    is_train = np.random.default_rng(split_seed).random(count) < train_fraction
    if not np.any(is_train & truth_selected) or not np.any(~is_train & truth_selected):
        raise ValueError("Training and held-out splits each need selected truth support")
    target_weight = weight * (target_lumi / source_lumi)
    with np.errstate(invalid="ignore", over="ignore"):
        truth_pt2 = np.full(count, np.nan)
        truth_pt2[available_truth] = arrays["true_pT"][available_truth] ** 2
        reco_pt2 = arrays["reco_pT"] ** 2
    if not np.all(np.isfinite(truth_pt2[available_truth])):
        raise ValueError("Squaring true_pT produced nonfinite generated values")
    # A finite but extreme malformed reco value can overflow on squaring.
    # Preserve the generated row and classify this as reconstruction failure.
    reco_selected &= np.isfinite(reco_pt2)
    reco_only_selected &= np.isfinite(reco_pt2)
    reco_only_before_proton &= np.isfinite(reco_pt2)
    reco_before_proton = truth_selected & reco_only_before_proton
    proton_cutflow = {
        "scope": "Shared reconstructed selection before observable-specific validity cuts; truth-fiducial numerator separately recorded",
        "reco_only_before_rows": int(reco_only_before_proton.sum()),
        "reco_only_after_rows": int(reco_only_selected.sum()),
        "reco_only_before_target_yield": float(target_weight[reco_only_before_proton].sum()),
        "reco_only_after_target_yield": float(target_weight[reco_only_selected].sum()),
        "truth_fiducial_reco_before_rows": int(reco_before_proton.sum()),
        "truth_fiducial_reco_after_rows": int(reco_selected.sum()),
        "truth_fiducial_reco_before_target_yield": float(target_weight[reco_before_proton].sum()),
        "truth_fiducial_reco_after_target_yield": float(target_weight[reco_selected].sum()),
        "truth_fiducial_rejected_rows": int((reco_before_proton & ~proton_selected).sum()),
        "truth_fiducial_rejected_target_yield": float(target_weight[reco_before_proton & ~proton_selected].sum()),
        "truth_denominator_unchanged": True, "generated_normalization_unchanged": True,
    }
    invalid_true_t = ~np.isfinite(arrays["true_t"]) | (arrays["true_t"] < 0)
    warnings = []
    if np.any(missing_truth):
        warnings.append("Coherent -999 producer truth rows retained in full generated normalization and "
                        "excluded from the truth-fiducial study without renormalization. Missing truth "
                        "does not establish that these events lie outside physical phase space; "
                        "any bias from producer truth-finding failures is not covered by statistical toys.")
    if q2_outlier.any():
        warnings.append("Small truth-Q2 outlier population retained in generated normalization and excluded by selection")
    if np.any(~reco_valid):
        warnings.append("Invalid reconstructed values retained as generated truth rows and classified as reconstruction failures")
    if np.any(invalid_true_t & available_truth):
        warnings.append("Invalid true_t retained without repair; abs_t_via_pt2 must reject invalid selected truth, while pt2 is independent of this branch")
    has_y_arrays = all(name in arrays for name in ("true_y", "reco_y"))
    if has_y_arrays:
        invalid_true_y = ~np.isfinite(arrays["true_y"]) | (arrays["true_y"] == -999.0)
        invalid_reco_y = ~np.isfinite(arrays["reco_y"]) | (arrays["reco_y"] == -999.0)
        if np.any(invalid_true_y & truth_selected):
            unavailable_observables["y"] = ("Invalid selected generated Upsilon rapidity: %d rows; no sentinel/absolute-value repair applied"
                % int(np.sum(invalid_true_y & truth_selected)))
            warnings.append("Invalid selected true_y retained without repair; rapidity extraction must reject this input, while other observables remain independent")
        if np.any(invalid_reco_y & reco_selected):
            warnings.append("Invalid reconstructed rapidity classified as reconstruction failure for the rapidity observable only")
    metadata = {
        "input": str(source_path), "tree": tree_name, "state": state["state"],
        "source_luminosity_invfb": source_lumi, "target_luminosity_invfb": target_lumi,
        "visible_cross_section_fb": state["visible_xsec_fb"], "branching_fraction_ee": state["br_ee"],
        "electron_energy_GeV": ee, "proton_energy_GeV": ep,
        "normalization_profile": profile, "normalization_metadata": normalization,
        "proton_selection": proton_policy, "proton_selection_cutflow": proton_cutflow,
        "normalization_tolerance": {"relative": NORMALIZATION_RTOL, "absolute": NORMALIZATION_ATOL,
                                    "generated_count_absolute": 1e-6},
        "rows": count, "all_generated_source_weight_sum": float(weight.sum()),
        "all_generated_target_weight_sum": float(target_weight.sum()),
        "target_expected_visible_generated_yield": state["visible_xsec_fb"] * target_lumi,
        "truth_schema_policy": TRUTH_SCHEMA_POLICY,
        "truth_schema_evidence": "DVMP_Upsilon_Analysis(3).C: reset block at lines 449-454; "
                                 "fill inside foundGenEvent at 531-541; unconditional outTree->Fill at 796. "
                                 "This checks the tuple pattern, not the production executable identity.",
        "available_truth_rows": available_count,
        "available_truth_fraction": available_fraction,
        "min_available_truth_fraction": min_available_fraction,
        "missing_truth_rows": int(missing_truth.sum()),
        "missing_truth_fraction": float(np.mean(missing_truth)),
        "missing_truth_target_yield": float(target_weight[missing_truth].sum()),
        "missing_truth_source_weight_sum": float(weight[missing_truth].sum()),
        "missing_truth_example_row_indices": np.flatnonzero(missing_truth)[:10].tolist(),
        "missing_truth_reco_only_selected_rows": int((missing_truth & reco_only_selected).sum()),
        "missing_truth_reco_only_selected_target_yield": float(target_weight[missing_truth & reco_only_selected].sum()),
        "truth_diagnostic_fraction_denominator": "Rows with available truth, excluding only coherent -999 blocks",
        "pure_state_fraction_all_generated": float(np.mean(state_mask)),
        "truth_selected_rows": int(truth_selected.sum()), "reco_selected_rows": int(reco_selected.sum()),
        "truth_selected_target_yield": float(target_weight[truth_selected].sum()),
        "reco_selected_target_yield": float(target_weight[reco_selected].sum()),
        "pure_state_fraction": state_fraction, "q2_outlier_fraction": q2_outlier_fraction,
        "q2_outlier_rows": int(q2_outlier.sum()), "q2_outlier_fraction_all_generated": float(np.mean(q2_outlier)),
        "nonfinite_reco_rows": int((~reco_finite).sum()),
        "negative_reco_pt_rows": int((arrays["reco_pT"] < 0).sum()),
        "invalid_true_t_rows": int(invalid_true_t.sum()),
        "invalid_true_t_fiducial_rows": int((invalid_true_t & truth_selected).sum()),
        "reco_only_selected_rows": int(reco_only_selected.sum()),
        "reco_only_excluded_by_truth_fiducial_rows": int((reco_only_selected & ~truth_selected).sum()),
        "reco_only_excluded_by_truth_fiducial_target_yield": float(target_weight[reco_only_selected & ~truth_selected].sum()),
        "reco_only_excluded_by_truth_W_rows": int((reco_only_selected & ((arrays["true_W"] < w_min) | (arrays["true_W"] > w_max))).sum()),
        "reco_only_excluded_by_truth_Q2_rows": int((reco_only_selected & ((arrays["true_Q2"] < 0) | (arrays["true_Q2"] > state["q2_max"]))).sum()),
        "reco_only_excluded_by_truth_mass_rows": int((reco_only_selected & ~state_mask).sum()),
        "truth_exclusion_counts_overlap": True,
        "truth_selection": {"mass_GeV": [state["mass"] - truth_window, state["mass"] + truth_window],
                            "W_GeV": [w_min, w_max], "Q2_GeV2": [0.0, state["q2_max"]]},
        "reco_selection": {"mass_GeV": [state["mass"] - reco_width * state["sigma"],
                                        state["mass"] + reco_width * state["sigma"]],
                           "W_GeV": [w_min, w_max], "sigma_epz_GeV": [epz_center - epz_halfwidth, epz_center + epz_halfwidth],
                           "requires_truth_fiducial": True, "fixed_mass_sigma_GeV": state["sigma"],
                           "proton_selection": proton_policy},
        "ordinary_weight_used": "weight", "weight_proj_used": False,
        "requested_observables": list(requested_observables),
        "unavailable_observables": unavailable_observables,
        "normalization_scope": "All generated rows before any analysis cut; target/source luminosity scaling only",
        "analysis_scope": "Conditional truth-fiducial differential study. Out-of-W/Q2/truth-state feed-in is not a detector-only acceptance model. Outside-variable-range feed-in is handled separately by response.py.",
        "production_provenance": "Assigned-rate consistency checked; EOS production executable is not authenticated",
        "train_fraction": train_fraction, "split_seed": split_seed,
        "train_truth_rows": int((is_train & truth_selected).sum()),
        "heldout_truth_rows": int((~is_train & truth_selected).sum()),
        "split_scaling": "Use configured split probability, not selected row count, to restore target exposure",
        "warnings": warnings,
    }
    if config.get("hash_inputs", False):
        metadata["input_sha256"] = _sha256(source_path)
    result = {"truth_pt2": truth_pt2, "truth_abs_t": arrays["true_t"], "reco_pt2": reco_pt2,
            "truth_w": arrays["true_W"], "reco_w": arrays["reco_W"],
            "truth_selected": truth_selected, "reco_selected": reco_selected,
            "weight": target_weight, "is_train": is_train, "metadata": metadata,
            "source_lumi": source_lumi, "target_lumi": target_lumi, "train_fraction": train_fraction}
    if has_y_arrays:
        result.update(truth_y=arrays["true_y"], reco_y=arrays["reco_y"],
                      reco_selected_y=reco_selected & ~invalid_reco_y)
        metadata.update(invalid_true_y_fiducial_rows=int((invalid_true_y & truth_selected).sum()),
                        invalid_reco_y_selected_rows=int((invalid_reco_y & reco_selected).sum()),
                        rapidity_definition="Upsilon laboratory rapidity from true_y/reco_y; negative physical rapidity is allowed; -999 is invalid")
    return result
