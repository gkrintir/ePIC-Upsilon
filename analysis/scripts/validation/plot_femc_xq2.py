#!/usr/bin/env python3
"""Audit and plot fEMC selected-yield losses in truth x and Q2.

Input tuples are read only. The visible rate already includes the dielectron
branching fraction. Every generated row remains in its sample's denominator.
Loss maps are central MC predictions; unknown cross-sample covariance is not
replaced by a fictitious binomial or independent-Poisson uncertainty.
"""
import argparse

# Explicit reconstructed proton category; truth rows remain unchanged.
import sys as _proton_sys
from pathlib import Path as _ProtonPath
_proton_sys.path.insert(0, str(_ProtonPath(__file__).resolve().parents[2] / "workflow"))
from proton_selection import proton_tag_mask, proton_branches, proton_selection_metadata
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.stats import beta
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle


MASSES = {"1S": 9.46030, "2S": 10.02326, "3S": 10.35520}
PROTON_MASS = 0.9382720813
COVERAGE = math.erf(1.0 / math.sqrt(2.0))
REQUIRED = ("weight", "true_mass", "true_W", "true_Q2", "reco_mass", "reco_W", "sigma_epz") + tuple(proton_branches())
SAMPLE_KEYS = ("nominal", "remove_inner", "remove_outer")
X_LABELS = {"xb": r"True Bjorken $x_B$", "xv": r"True $x_V$ ($t\simeq0$ proxy)"}


def clean_json(value):
    if isinstance(value, dict):
        return {str(k): clean_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean_json(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def write_json(path, payload):
    path.write_text(json.dumps(clean_json(payload), indent=2, allow_nan=False) + "\n")


def coordinates(q2, w, mass):
    q2, w = np.asarray(q2, dtype=float), np.asarray(w, dtype=float)
    denominator = w * w + q2 - PROTON_MASS * PROTON_MASS
    valid = np.isfinite(denominator) & (denominator > 0)
    xb = np.divide(q2, denominator, out=np.full_like(q2, np.nan), where=valid)
    xv = np.divide(q2 + mass * mass, denominator, out=np.full_like(q2, np.nan), where=valid)
    return xb, xv


def selection_masks(arrays, config):
    mass = MASSES[config["state"]]
    q2, w = arrays["true_Q2"], arrays["true_W"]
    truth = (np.isfinite(arrays["true_mass"])
             & (np.abs(arrays["true_mass"] - mass) <= config["truth_mass_window"])
             & np.isfinite(q2) & (q2 >= 0) & (q2 <= config["q2_max"])
             & np.isfinite(w))
    if config["selection_scope"] == "analysis":
        truth &= (w >= config["w_min"]) & (w <= config["w_max"])
    else:
        truth &= w >= mass + PROTON_MASS
    reco = (truth & np.isfinite(arrays["reco_mass"])
            & (arrays["reco_mass"] >= 8) & (arrays["reco_mass"] <= 12)
            & np.isfinite(arrays["reco_W"])
            & np.isfinite(arrays["sigma_epz"])
            & (np.abs(arrays["sigma_epz"] - 2 * config["electron_energy"])
               <= 0.2 * config["electron_energy"]))
    if config["selection_scope"] == "analysis":
        reco &= (arrays["reco_W"] >= config["w_min"]) & (arrays["reco_W"] <= config["w_max"])
    else:
        reco &= arrays["reco_W"] > 0
    reco &= proton_tag_mask(arrays)
    return truth, reco


def clopper_pearson(success, total, coverage=COVERAGE):
    success, total = np.broadcast_arrays(np.asarray(success, dtype=float), np.asarray(total, dtype=float))
    if np.any(success < 0) or np.any(total < success) or np.any(success != np.floor(success)) or np.any(total != np.floor(total)):
        raise ValueError("Clopper-Pearson requires integer 0 <= selected <= truth counts")
    alpha = 1.0 - coverage
    low, high = np.full(success.shape, np.nan), np.full(success.shape, np.nan)
    positive = total > 0
    low[positive & (success == 0)] = 0.
    high[positive & (success == total)] = 1.
    use = positive & (success > 0)
    low[use] = beta.ppf(alpha / 2, success[use], total[use] - success[use] + 1)
    use = positive & (success < total)
    high[use] = beta.ppf(1 - alpha / 2, success[use] + 1, total[use] - success[use])
    return low, high


def comparison_metrics(nominal_counts, scenario_counts, nominal_ngen, scenario_ngen):
    baseline = np.asarray(nominal_counts, dtype=float) / nominal_ngen
    scenario = np.asarray(scenario_counts, dtype=float) / scenario_ngen
    ratio = np.divide(scenario, baseline, out=np.full_like(baseline, np.nan), where=baseline > 0)
    return 100. * (1. - ratio)


def validate_edges(value, label):
    edges = np.asarray(value, dtype=float)
    if edges.ndim != 1 or len(edges) < 2 or np.any(~np.isfinite(edges)) or np.any(edges <= 0) or np.any(np.diff(edges) <= 0):
        raise ValueError(label + " must contain positive, finite, strictly increasing log-axis edges")
    return edges


def validate_config(config):
    result = dict(config)
    result.setdefault("tree", "UpsilonTree")
    result.setdefault("selection_scope", "analysis")
    result.setdefault("min_nominal_mc", 25)
    result.setdefault("truth_mass_window", .2)
    result.setdefault("hash_inputs", False)
    if result.get("state") not in MASSES:
        raise ValueError("state must be 1S, 2S or 3S")
    if result["selection_scope"] not in ("analysis", "full"):
        raise ValueError("selection_scope must be analysis or full")
    for key in ("electron_energy", "proton_energy", "visible_xsec_fb", "lumi_invfb", "q2_max", "truth_mass_window"):
        value = float(result[key])
        if not math.isfinite(value) or value <= 0:
            raise ValueError(key + " must be positive and finite")
        result[key] = value
    result["w_min"], result["w_max"] = float(result["w_min"]), float(result["w_max"])
    if not math.isfinite(result["w_min"]) or not math.isfinite(result["w_max"]) or not 0 < result["w_min"] < result["w_max"]:
        raise ValueError("Require positive ordered w_min,w_max")
    if not isinstance(result["min_nominal_mc"], int) or isinstance(result["min_nominal_mc"], bool) or result["min_nominal_mc"] < 1:
        raise ValueError("min_nominal_mc must be a positive integer")
    samples = result.get("samples", [])
    if len(samples) != 3 or sorted(s.get("key", "") for s in samples) != sorted(SAMPLE_KEYS):
        raise ValueError("Exactly nominal, remove_inner and remove_outer samples are required")
    ordered = {s["key"]: dict(s) for s in samples}
    paths = []
    for key in SAMPLE_KEYS:
        sample = ordered[key]
        path = Path(sample["path"]).expanduser().resolve()
        if not path.is_file():
            raise ValueError("Missing input: " + str(path))
        sample["path"] = str(path)
        sample.setdefault("label", key.replace("_", " "))
        paths.append(path)
    if len(set(paths)) != 3:
        raise ValueError("Nominal and removal inputs must be three distinct files")
    result["samples"] = [ordered[k] for k in SAMPLE_KEYS]
    for key in ("q2_edges", "xb_edges", "xv_edges"):
        if result.get(key) is not None:
            result[key] = validate_edges(result[key], key).tolist()
    return result


def event_count(source, tree, path):
    nrows = int(tree.num_entries)
    if nrows <= 0:
        raise ValueError(str(path) + ": empty tuple")
    if "h_EventCount" not in source:
        raise ValueError(str(path) + ": missing h_EventCount; cannot establish all-generated normalization")
    regular = np.asarray(source["h_EventCount"].values(flow=False), dtype=float)
    flows = np.asarray(source["h_EventCount"].values(flow=True), dtype=float)
    if regular.ndim != 1 or regular.size < 1 or flows.ndim != 1 or flows.size != regular.size + 2:
        raise ValueError(str(path) + ": h_EventCount must be a one-dimensional event counter")
    if np.any(~np.isfinite(flows)) or np.any(flows < 0) or flows[0] != 0 or flows[-1] != 0:
        raise ValueError(str(path) + ": invalid h_EventCount values or nonzero flows")
    ngen = float(regular.sum())
    if not math.isfinite(ngen) or ngen <= 0 or not math.isclose(ngen, nrows, rel_tol=0, abs_tol=1e-6):
        raise ValueError(str(path) + ": generated count does not equal all tuple rows")
    return nrows


def chunks(tree):
    for values in tree.iterate(expressions=list(REQUIRED), step_size="64 MB", library="np"):
        if isinstance(values, np.ndarray) and values.dtype.names:
            values = {key: values[key] for key in values.dtype.names}
        arrays = {key: np.asarray(values[key], dtype=float) for key in REQUIRED}
        if any(a.ndim != 1 for a in arrays.values()):
            raise ValueError("All required branches must be scalar per-event quantities")
        yield arrays


def checked_tree(source, config, path):
    if config["tree"] not in source:
        raise ValueError(str(path) + ": missing tree " + config["tree"])
    tree = source[config["tree"]]
    missing = sorted(set(REQUIRED) - set(tree.keys()))
    if missing:
        raise ValueError(str(path) + ": missing branches " + ", ".join(missing))
    return tree


def discover_edges(config):
    nominal = config["samples"][0]
    min_q, min_xb, max_xb = float("inf"), float("inf"), 0.
    if config.get("q2_edges") is None or config.get("xb_edges") is None:
        with uproot.open(nominal["path"]) as source:
            tree = checked_tree(source, config, nominal["path"])
            event_count(source, tree, nominal["path"])
            for arrays in chunks(tree):
                truth, _ = selection_masks(arrays, config)
                xb, _ = coordinates(arrays["true_Q2"], arrays["true_W"], MASSES[config["state"]])
                positive = truth & (arrays["true_Q2"] > 0) & np.isfinite(xb) & (xb > 0)
                if np.any(positive):
                    min_q = min(min_q, float(arrays["true_Q2"][positive].min()))
                    min_xb = min(min_xb, float(xb[positive].min()))
                    max_xb = max(max_xb, float(xb[positive].max()))
        if not math.isfinite(min_q) or not math.isfinite(min_xb):
            raise ValueError("Nominal truth has no positive Q2/xB support for automatic logarithmic binning; supply explicit edges to audit a zero-Q2-only sample")
    if config.get("q2_edges") is None:
        lower = 10. ** math.floor(math.log10(min_q))
        upper = config["q2_max"]
        if lower >= upper:
            lower = upper / 10.
        qedges = np.geomspace(lower, upper, 7)
    else:
        qedges = validate_edges(config["q2_edges"], "q2_edges")
    if config.get("xb_edges") is None:
        lower = 10. ** math.floor(math.log10(min_xb))
        threshold = config["w_min"] if config["selection_scope"] == "analysis" else MASSES[config["state"]] + PROTON_MASS
        expected_max = config["q2_max"] / (threshold ** 2 + config["q2_max"] - PROTON_MASS ** 2)
        upper = min(1., max(expected_max, max_xb))
        if lower >= upper:
            lower = upper / 10.
        xbedges = np.geomspace(lower, upper, 8)
    else:
        xbedges = validate_edges(config["xb_edges"], "xb_edges")
    xvedges = validate_edges(config.get("xv_edges") or [.005, .01, .02, .05, .1, .2, .4, .8, 1.], "xv_edges")
    return {"q2": qedges, "xb": xbedges, "xv": xvedges}


def update_flows(flows, x, q, mask, xedges, qedges):
    finite = np.isfinite(x) & np.isfinite(q)
    in_grid = finite & (x >= xedges[0]) & (x <= xedges[-1]) & (q >= qedges[0]) & (q <= qedges[-1])
    checks = {"rows": mask, "in_grid": mask & in_grid, "outside_grid": mask & ~in_grid,
              "q2_zero": mask & (q == 0), "q2_underflow": mask & finite & (q < qedges[0]),
              "q2_overflow": mask & finite & (q > qedges[-1]),
              "x_underflow": mask & finite & (x < xedges[0]),
              "x_overflow": mask & finite & (x > xedges[-1]),
              "nonfinite_coordinate": mask & ~finite}
    for key, use in checks.items():
        flows[key] = flows.get(key, 0) + int(np.count_nonzero(use))


def load_histograms(config, sample, edges):
    path = Path(sample["path"])
    hist = {coord: {kind: np.zeros((len(edges[coord])-1, len(edges["q2"])-1), dtype=np.int64)
                    for kind in ("truth", "selected")} for coord in ("xb", "xv")}
    flow = {coord: {kind: {} for kind in ("truth", "selected")} for coord in hist}
    state_rows = invalid_q2 = above_q2 = truth_rows = selected_rows = rows_seen = 0
    sum_raw_weight = 0.
    weight_ref = None
    weight_min, weight_max = float("inf"), -float("inf")
    q2_min, q2_max = float("inf"), -float("inf")
    with uproot.open(path) as source:
        tree = checked_tree(source, config, path)
        ngen = event_count(source, tree, path)
        for arrays in chunks(tree):
            weight = arrays["weight"]
            if np.any(~np.isfinite(weight)) or np.any(weight <= 0):
                raise ValueError(str(path) + ": generated flat normalization requires finite positive base weight on every row")
            if weight_ref is None:
                weight_ref = float(weight[0])
            if not np.allclose(weight, weight_ref, rtol=1e-8, atol=0):
                raise ValueError(str(path) + ": nonuniform base weight; flat generated normalization is unsupported for weighted-event shapes")
            weight_min, weight_max = min(weight_min, float(weight.min())), max(weight_max, float(weight.max()))
            sum_raw_weight += float(weight.sum())
            rows_seen += len(weight)
            mass = arrays["true_mass"]
            q2 = arrays["true_Q2"]
            state_rows += int(np.count_nonzero(np.isfinite(mass) & (np.abs(mass-MASSES[config["state"]]) <= config["truth_mass_window"])))
            invalid_q2 += int(np.count_nonzero(~np.isfinite(q2) | (q2 < 0)))
            above_q2 += int(np.count_nonzero(np.isfinite(q2) & (q2 > config["q2_max"] * (1+1e-6))))
            if np.any(np.isfinite(q2)):
                finite_q2 = q2[np.isfinite(q2)]
                q2_min, q2_max = min(q2_min, float(finite_q2.min())), max(q2_max, float(finite_q2.max()))
            truth, selected = selection_masks(arrays, config)
            truth_rows += int(truth.sum())
            selected_rows += int(selected.sum())
            xb, xv = coordinates(q2, arrays["true_W"], MASSES[config["state"]])
            for coord, x in (("xb", xb), ("xv", xv)):
                for kind, mask in (("truth", truth), ("selected", selected)):
                    good = mask & np.isfinite(x) & np.isfinite(q2)
                    hist[coord][kind] += np.histogram2d(x[good], q2[good], bins=(edges[coord], edges["q2"]))[0].astype(np.int64)
                    update_flows(flow[coord][kind], x, q2, mask, edges[coord], edges["q2"])
    if rows_seen != ngen:
        raise ValueError(str(path) + ": iteration row count disagrees with generated denominator")
    if state_rows / ngen < .95:
        raise ValueError(str(path) + ": selected-state truth-mass population is below 95%; check the input state")
    if (invalid_q2 + above_q2) / ngen > .001:
        raise ValueError(str(path) + ": more than 0.1% invalid/out-of-profile truth Q2; check the generation phase space")
    event_weight = config["visible_xsec_fb"] * config["lumi_invfb"] / ngen
    metadata = {"key": sample["key"], "label": sample["label"], "path": str(path),
                "input_rows": ngen, "ngen": ngen, "event_weight": event_weight,
                "state_truth_rows": state_rows, "truth_rows": truth_rows, "selected_rows": selected_rows,
                "truth_fraction_of_generated": truth_rows / ngen,
                "selected_fraction_of_generated": selected_rows / ngen,
                "projected_selected_yield": selected_rows * event_weight,
                "projected_all_generated_yield": ngen * event_weight,
                "visible_xsec_fb": config["visible_xsec_fb"], "lumi_invfb": config["lumi_invfb"],
                "raw_weight_min": weight_min, "raw_weight_max": weight_max, "raw_sum_weight": sum_raw_weight,
                "raw_weights_used_for_yield": False, "weight_proj_used": False,
                "invalid_q2_rows": invalid_q2, "above_profile_q2_rows": above_q2,
                "true_q2_min": q2_min, "true_q2_max": q2_max, "flows": flow}
    stat = path.stat()
    metadata["input_fingerprint"] = {"size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if config["hash_inputs"]:
        h = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(8*1024*1024), b""):
                h.update(block)
        metadata["input_fingerprint"]["sha256"] = h.hexdigest()
    return {"metadata": metadata, "histograms": hist}


def make_rows(samples, coord, edges, config):
    baseline = samples[0]
    btruth = baseline["histograms"][coord]["truth"]
    breco = baseline["histograms"][coord]["selected"]
    beff = np.divide(breco, btruth, out=np.full(btruth.shape, np.nan), where=btruth > 0)
    rows = []
    for sample in samples:
        truth = sample["histograms"][coord]["truth"]
        reco = sample["histograms"][coord]["selected"]
        metadata = sample["metadata"]
        eff = np.divide(reco, truth, out=np.full(truth.shape, np.nan), where=truth > 0)
        low, high = clopper_pearson(reco, truth)
        loss = comparison_metrics(breco, reco, baseline["metadata"]["ngen"], metadata["ngen"])
        eff_ratio_loss = 100 * (1 - np.divide(eff, beff, out=np.full(eff.shape, np.nan), where=beff > 0))
        for ix, iq in np.ndindex(truth.shape):
            flags = []
            if breco[ix, iq] == 0:
                flags.append("UNDEFINED_NOMINAL_ZERO")
            elif breco[ix, iq] < config["min_nominal_mc"]:
                flags.append("LOW_NOMINAL_MC")
            if truth[ix, iq] == 0:
                flags.append("NO_SCENARIO_TRUTH")
            elif reco[ix, iq] == 0:
                flags.append("ZERO_SCENARIO_SELECTED")
            if np.isfinite(loss[ix, iq]) and loss[ix, iq] < 0:
                flags.append("GAIN")
            rows.append({"coordinate": coord, "sample": metadata["key"], "label": metadata["label"],
                         "x_low": edges[coord][ix], "x_high": edges[coord][ix+1],
                         "q2_low": edges["q2"][iq], "q2_high": edges["q2"][iq+1],
                         "truth_mc": int(truth[ix, iq]), "selected_mc": int(reco[ix, iq]),
                         "ngen": metadata["ngen"], "event_weight": metadata["event_weight"],
                         "projected_selected_yield": reco[ix, iq]*metadata["event_weight"],
                         "projected_truth_yield": truth[ix, iq]*metadata["event_weight"],
                         "efficiency": eff[ix, iq], "efficiency_cp_low": low[ix, iq], "efficiency_cp_high": high[ix, iq],
                         "efficiency_interval_coverage": COVERAGE, "baseline_selected_mc": int(breco[ix, iq]),
                         "loss_percent": loss[ix, iq], "efficiency_ratio_loss_percent": eff_ratio_loss[ix, iq],
                         "loss_interval": "NOT_EVALUATED_UNKNOWN_INTERSAMPLE_COVARIANCE", "flag": "|".join(flags) or "OK"})
    return rows


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ("" if isinstance(v, (float, np.floating)) and not np.isfinite(v) else v) for k, v in row.items()})


def describe_cuts(config):
    if config["selection_scope"] == "analysis":
        wtext = rf"${config['w_min']:g}\leq W_{{\rm true,reco}}\leq {config['w_max']:g}$ GeV"
    else:
        wtext = r"$W_{\rm true}\geq M_\Upsilon+m_p$; $W_{\rm reco}>0$"
    return (wtext + rf", $Q^2_{{\rm true}}\leq {config['q2_max']:g}$ GeV$^2$; "
            r"$8\leq m_{ee}^{\rm reco}\leq12$ GeV; "
            r"$|\Sigma(E-p_z)-2E_e|\leq0.2E_e$")


def setup_axis(ax, coord, edges):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(edges[coord][0], edges[coord][-1])
    ax.set_ylim(edges["q2"][0], edges["q2"][-1])
    ax.set_xlabel(X_LABELS[coord])
    ax.set_ylabel(r"True $Q^2$ [GeV$^2$]")
    ax.set_facecolor("#eeeeee")
    ax.tick_params(which="both", direction="in", top=True, right=True)


def annotate(ax, values, coord, edges, percent=False, low=None):
    for ix, iq in np.ndindex(values.shape):
        v = values[ix, iq]
        if not np.isfinite(v) or (not percent and v <= 0):
            continue
        x = math.sqrt(edges[coord][ix] * edges[coord][ix+1])
        y = math.sqrt(edges["q2"][iq] * edges["q2"][iq+1])
        if percent:
            label = ("%.0f%%" % v) if abs(v) >= 10 else ("%.1f%%" % v)
        else:
            label = "%.2g" % v
        if low is not None and low[ix, iq]:
            label += "*"
        ax.text(x, y, label, ha="center", va="center", fontsize=7,
                bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": .75})


def loss_scale(losses):
    finite = np.concatenate([a[np.isfinite(a)] for a in losses])
    minimum = float(finite.min()) if finite.size else 0.
    if minimum < 0:
        return TwoSlopeNorm(vmin=min(-1., minimum), vcenter=0., vmax=100.), "RdYlBu_r"
    return Normalize(vmin=0., vmax=100.), "viridis"


def draw_loss(ax, loss, baseline_counts, coord, edges, config, scale, cmap):
    setup_axis(ax, coord, edges)
    mesh = ax.pcolormesh(edges[coord], edges["q2"], np.ma.masked_invalid(loss.T),
                         norm=scale, cmap=cmap, shading="flat", edgecolors="white", linewidth=.35)
    low = (baseline_counts > 0) & (baseline_counts < config["min_nominal_mc"])
    for ix, iq in zip(*np.nonzero(low)):
        ax.add_patch(Rectangle((edges[coord][ix], edges["q2"][iq]),
                              edges[coord][ix+1]-edges[coord][ix], edges["q2"][iq+1]-edges["q2"][iq],
                              facecolor="none", edgecolor="#555555", linewidth=0, hatch="///", alpha=.3))
    annotate(ax, loss, coord, edges, percent=True, low=low)
    return mesh


def save_figure(fig, output, stem):
    paths = []
    for suffix in ("png", "pdf"):
        path = output / (stem + "." + suffix)
        fig.savefig(path, dpi=180, facecolor="white", bbox_inches="tight")
        paths.append(path.name)
    plt.close(fig)
    return paths


def render_maps(samples, coord, edges, config, output):
    baseline = samples[0]
    counts = baseline["histograms"][coord]["selected"]
    nominal_yield = counts * baseline["metadata"]["event_weight"]
    losses = [comparison_metrics(counts, sample["histograms"][coord]["selected"],
                                baseline["metadata"]["ngen"], sample["metadata"]["ngen"])
              for sample in samples[1:]]
    scale, cmap = loss_scale(losses)
    title = rf"$\Upsilon({config['state']})\to e^+e^-$  |  ${config['electron_energy']:g}\times{config['proton_energy']:g}$ GeV"
    foot = (rf"Visible selected yields at ${config['lumi_invfb']:g}$ fb$^{{-1}}$. "
            r"Loss $=100(1-N_{\rm scenario}/N_{\rm nominal})$; negative values are gains."
            "\nGray loss bins have zero nominal selected MC. Hatched/* bins have fewer than "
            + str(config["min_nominal_mc"]) + " nominal selected MC events. Loss-ratio uncertainties are not evaluated.")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12,
                         "axes.labelsize": 12, "savefig.facecolor": "white"})
    fig, axes = plt.subplots(1, 3, figsize=(18.5, 6.3))
    setup_axis(axes[0], coord, edges)
    positive = nominal_yield[nominal_yield > 0]
    if positive.size:
        lo, hi = float(positive.min()), float(positive.max())
        mesh = axes[0].pcolormesh(edges[coord], edges["q2"], np.ma.masked_less_equal(nominal_yield.T, 0),
                                  norm=LogNorm(vmin=lo, vmax=max(hi, lo*1.01)), cmap="viridis",
                                  shading="flat", edgecolors="white", linewidth=.35)
        fig.colorbar(mesh, ax=axes[0], pad=.02, label="Expected selected events / bin")
    else:
        axes[0].text(.5, .5, "No nominal selected events", transform=axes[0].transAxes, ha="center")
    annotate(axes[0], nominal_yield, coord, edges)
    axes[0].set_title("Nominal projected yield")
    for ax, loss, sample in zip(axes[1:], losses, samples[1:]):
        mesh = draw_loss(ax, loss, counts, coord, edges, config, scale, cmap)
        ax.set_title(sample["metadata"]["label"])
        fig.colorbar(mesh, ax=ax, pad=.02, label="Relative selected-yield loss [%]")
    fig.suptitle(title, fontsize=17, y=.98)
    fig.text(.5, .915, describe_cuts(config), ha="center", fontsize=10)
    fig.text(.5, .02, foot, ha="center", fontsize=9)
    fig.subplots_adjust(top=.83, bottom=.22, left=.06, right=.98, wspace=.28)
    paths = save_figure(fig, output, coord + "_q2_comparison")
    for loss, sample in zip(losses, samples[1:]):
        fig, ax = plt.subplots(figsize=(8.8, 7.3))
        mesh = draw_loss(ax, loss, counts, coord, edges, config, scale, cmap)
        fig.colorbar(mesh, ax=ax, pad=.025, label="Relative selected-yield loss [%]")
        ax.set_title(title + "\n" + sample["metadata"]["label"], fontsize=15, pad=16)
        fig.text(.5, .105, describe_cuts(config), ha="center", fontsize=8)
        fig.text(.5, .025, foot, ha="center", fontsize=7)
        fig.subplots_adjust(top=.86, bottom=.24, left=.12, right=.92)
        paths += save_figure(fig, output, coord + "_q2_loss_" + sample["metadata"]["key"])
    return paths


def run(config, out_dir):
    config = validate_config(config)
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=True)
    edges = discover_edges(config)
    samples = [load_histograms(config, sample, edges) for sample in config["samples"]]
    warnings, comparisons, fraction_comparisons = [], [], []
    baseline = samples[0]
    for sample in samples:
        meta = sample["metadata"]
        if meta["invalid_q2_rows"] or meta["above_profile_q2_rows"]:
            warnings.append(meta["key"] + ": small Q2 outlier population retained in normalization denominator but excluded by selection")
        for coord in ("xb", "xv"):
            flow = meta["flows"][coord]["selected"]
            if flow["outside_grid"]:
                warnings.append("%s %s: %d selected MC rows outside plotted log grid; see flow audit" % (meta["key"], coord, flow["outside_grid"]))
        if meta["key"] == "nominal":
            continue
        nominal_fraction = baseline["metadata"]["truth_fraction_of_generated"]
        scenario_fraction = meta["truth_fraction_of_generated"]
        fraction_ratio = scenario_fraction / nominal_fraction if nominal_fraction > 0 else float("nan")
        fraction_comparisons.append({
            "sample": meta["key"],
            "nominal_truth_fraction_of_generated": nominal_fraction,
            "scenario_truth_fraction_of_generated": scenario_fraction,
            "scenario_over_nominal_truth_fraction": fraction_ratio,
            "relative_discrepancy_threshold": .05,
            "interpretation": "Descriptive generated-population fraction check, not a p-value; unknown sample correlation",
        })
        if not np.isfinite(fraction_ratio) or abs(fraction_ratio - 1.) > .05:
            warnings.append(meta["key"] + ": truth selection fraction T/Ngen differs from nominal by more than 5% or nominal fraction is zero; detector-only interpretation is not established")
        for coord in ("xb", "xv"):
            # Include all outside-grid truth as an extra category; this diagnostic
            # is not a hypothesis test and assumes no intersample independence.
            a = np.append(baseline["histograms"][coord]["truth"].ravel(), baseline["metadata"]["flows"][coord]["truth"]["outside_grid"])
            b = np.append(sample["histograms"][coord]["truth"].ravel(), meta["flows"][coord]["truth"]["outside_grid"])
            if a.sum() and b.sum():
                distance = .5 * np.abs(a/a.sum() - b/b.sum()).sum()
            else:
                distance = float("nan")
            comparisons.append({"sample": meta["key"], "coordinate": coord, "truth_shape_total_variation": distance,
                                "threshold": .05, "interpretation": "Descriptive shape diagnostic, not a p-value; unknown sample correlation"})
            if not np.isfinite(distance) or distance > .05:
                warnings.append(meta["key"] + " " + coord + ": truth distributions differ or have no support; detector-only interpretation is not established")
    files, all_rows = [], []
    for coord in ("xb", "xv"):
        csv_path = output / (coord + "_q2_bins.csv")
        rows = make_rows(samples, coord, edges, config)
        all_rows.extend(rows)
        write_csv(csv_path, rows)
        files.append(csv_path.name)
        files += render_maps(samples, coord, edges, config, output)
    write_csv(output / "bins.csv", all_rows)
    files.append("bins.csv")
    audit = {"status": "WARN" if warnings else "PASS", "config": config, "samples": [s["metadata"] for s in samples],
             "edges": edges, "outputs": files, "warnings": warnings, "truth_comparisons": comparisons,
             "truth_fraction_comparisons": fraction_comparisons,
             "coordinate_definitions": {"xb": "Q2/(W_true^2+Q2-m_p^2)", "xv": "(Q2+M_Upsilon^2)/(W_true^2+Q2-m_p^2)",
                                        "xv_interpretation": "t=0 exclusive heavy-meson proxy, not Bjorken x or an exact gluon momentum fraction",
                                        "mass_GeV": MASSES[config["state"]], "proton_mass_GeV": PROTON_MASS},
             "normalization": "All generated rows: sigma_visible_fb * lumi_invfb / each sample Ngen. No BR division and no post-selection renormalization.",
             "proton_selection": proton_selection_metadata(),
             "selection": {"scope": config["selection_scope"], "truth_mass_half_window_GeV": config["truth_mass_window"],
                           "truth_q2_GeV2": [0, config["q2_max"]], "reco_mass_GeV": [8, 12],
                           "sigma_epz_half_width_GeV": .2*config["electron_energy"],
                           "true_reco_W_GeV": [config["w_min"], config["w_max"]] if config["selection_scope"] == "analysis" else None,
                           "full_scope_W": "truth W >= M_Upsilon+m_p; reco W > 0; no analysis W window"},
             "statistical_interpretation": {"yield": "Expected visible selected signal counts; no data realization or projected error bars",
                                           "efficiency_intervals": "Exact Clopper-Pearson intervals from selected/truth integer MC counts",
                                           "efficiency_interval_coverage": COVERAGE,
                                           "loss_intervals": "Not evaluated: intersample event identities and covariance are unknown",
                                           "low_mc": "Hatching flags low nominal selected MC, not low projected data yield",
                                           "efficiency_ratio_loss": "CSV-only diagnostic with each bin's truth population standardized; differs from projected-yield loss"},
             "production_comparability": "Same state, beam, generation phase space and electron-ID production must be independently established; supplied paths and truth shapes do not authenticate producer settings",
             "log_flow_note": "All positive-axis under/overflow and Q2=0 counts are audited. Axis flow categories can overlap; outside_grid is their non-overlapping union."}
    write_json(output / "audit.json", audit)
    print("%s: %s" % (audit["status"], output / "audit.json"))
    for sample in samples:
        m = sample["metadata"]
        print("%s: Ngen=%d, selected MC=%d, projected yield=%.8g" % (m["key"], m["ngen"], m["selected_rows"], m["projected_selected_yield"]))
    return audit


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        config = json.loads(args.config.read_text())
        run(config, args.out_dir)
        return 0
    except (OSError, ValueError, KeyError, TypeError) as error:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        write_json(args.out_dir / "audit.json", {"status": "FAIL", "error": str(error), "config_path": str(args.config)})
        print("ERROR: " + str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
