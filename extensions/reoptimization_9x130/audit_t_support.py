#!/usr/bin/env python3
"""Audit the 9x130 fixed-bin true-|t| response without fitting or altering inputs.

Reproduces the existing pT2-estimator selection and random split exactly, then
compares it with a pre-proton-tag selection and a physical reco-pT requirement.
Run this external script with --framework-root pointing at the tagged framework.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

# The diagnostic must not create import caches inside the supplied framework.
sys.dont_write_bytecode = True

import numpy as np
import uproot

STATES = {
    "1S": {"mass": 9.46030, "sigma": 0.053, "q2_max": 1.0, "visible_xsec_fb": 119.685},
    "2S": {"mass": 10.02326, "sigma": 0.056, "q2_max": 0.01, "visible_xsec_fb": 29.297},
    "3S": {"mass": 10.35520, "sigma": 0.059, "q2_max": 0.01, "visible_xsec_fb": 23.447},
}
BASE_BRANCHES = ("weight", "true_mass", "true_t", "true_W", "reco_mass",
                 "reco_W", "sigma_epz", "reco_pT", "true_Q2")
TAG_BRANCHES = ("reco_topology", "reco_tagged_proton_p")
EDGES = np.array([0., .2, .5, 1., 5.])


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValueError("Cannot load framework helper: " + str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def helpers(framework_root, output):
    root = Path(framework_root).expanduser().resolve()
    policy_module = load_module(root / "workflow/proton_selection.py", "_support_audit_proton")
    isolation = load_module(root / "workflow/proton_category_isolation.py", "_support_audit_isolation")
    marker = isolation.enforce_category_workspace(root, output)
    if marker is None:
        raise ValueError("This audit requires a completed isolated proton-tagged category marker")
    policy = policy_module.proton_selection_metadata()
    if policy.get("mode") != "required" or policy.get("contract") != "jul26-reco-topology-v1":
        raise ValueError("This audit requires the jul26-reco-topology-v1 proton selection in required mode")
    return root, policy_module, policy


def counts(mask, weights, target_scale):
    values = weights[mask]
    return {"rows": int(np.count_nonzero(mask)), "positive_weight_rows": int(np.count_nonzero(values > 0)),
            "source_sum_weight": float(values.sum()), "target_expected_yield": float(values.sum() * target_scale)}


def response(truth, reco, weights, selected, train, test, edges, projection_scale):
    nb = len(edges) - 1
    truth_train = np.histogram(truth[train], edges, weights=weights[train])[0]
    truth_test = np.histogram(truth[test], edges, weights=weights[test])[0] * projection_scale
    data_reco = np.histogram(reco[test & selected], edges, weights=weights[test & selected])[0] * projection_scale
    ib = np.digitize(truth, edges) - 1
    jb = np.digitize(reco, edges) - 1
    ib[truth == edges[-1]] = nb - 1
    jb[reco == edges[-1]] = nb - 1
    valid = train & selected & (ib >= 0) & (ib < nb) & (jb >= 0) & (jb < nb)
    matrix = np.zeros((nb, nb))
    np.add.at(matrix, (ib[valid], jb[valid]), weights[valid])
    probability = np.divide(matrix, truth_train[:, None], out=np.zeros_like(matrix), where=truth_train[:, None] > 0)
    missing_truth = np.flatnonzero((truth_test > 0) & (probability.sum(axis=1) <= 0))
    missing_reco = np.flatnonzero((data_reco > 0) & (probability.sum(axis=0) <= 0))
    return {"truth_train_source_sumw": truth_train.tolist(), "truth_test_target_yield": truth_test.tolist(),
            "reco_test_target_yield": data_reco.tolist(), "H_truth_rows_reco_columns": matrix.tolist(),
            "A_truth_rows_reco_columns": probability.tolist(), "efficiency": probability.sum(axis=1).tolist(),
            "unsupported_truth_bins": missing_truth.tolist(), "unsupported_reco_bins": missing_reco.tolist()}


def audit_arrays(arrays, args, policy_module, policy):
    props = STATES[args.state]
    arr = {name: np.asarray(arrays[name], dtype=float) for name in BASE_BRANCHES + TAG_BRANCHES}
    sizes = {value.shape for value in arr.values()}
    if len(sizes) != 1 or len(next(iter(sizes))) != 1:
        raise ValueError("All required branches must be scalar event arrays of equal length")
    n = len(arr["weight"])
    if not n:
        raise ValueError("Input tree has no rows")
    weight = arr["weight"]
    if not np.all(np.isfinite(weight)) or np.any(weight < 0):
        raise ValueError("Projected response requires finite nonnegative source weights")
    policy_module.proton_tag_mask(arr)  # Same schema validation as the production script.
    global_finite = np.logical_and.reduce([np.isfinite(arr[name]) for name in BASE_BRANCHES])
    mass = props["mass"]
    truth_raw = (np.isfinite(arr["true_mass"]) & (arr["true_mass"] >= mass - .2) & (arr["true_mass"] <= mass + .2)
                 & np.isfinite(arr["true_Q2"]) & (arr["true_Q2"] >= 0) & (arr["true_Q2"] <= props["q2_max"])
                 & np.isfinite(arr["true_W"]) & (arr["true_W"] >= 20.) & (arr["true_W"] <= 70.)
                 & np.isfinite(arr["true_t"]))
    truth_legacy = truth_raw & global_finite
    in_legacy_range = truth_legacy & (arr["true_t"] >= EDGES[0]) & (arr["true_t"] <= EDGES[-1])
    sub = {name: value[in_legacy_range] for name, value in arr.items()}
    if len(sub["weight"]) < 2:
        raise ValueError("Fewer than two legacy truth-fiducial rows remain")
    train = np.random.default_rng(args.seed).random(len(sub["weight"])) < args.train_fraction
    test = ~train
    if not np.any(train) or not np.any(test):
        raise ValueError("The declared random split leaves an empty training or held-out sample")
    tags = policy_module.proton_tag_mask(sub)
    with np.errstate(over="ignore", invalid="ignore"):
        reco_t = sub["reco_pT"] ** 2
    mass_cut = (sub["reco_mass"] >= mass - 3. * props["sigma"]) & (sub["reco_mass"] <= mass + 3. * props["sigma"])
    epz_cut = np.abs(sub["sigma_epz"] - 18.) <= 1.8
    w_cut = (sub["reco_W"] >= 20.) & (sub["reco_W"] <= 70.)
    estimator_range = (reco_t >= EDGES[0]) & (reco_t <= EDGES[-1])
    pretag = mass_cut & epz_cut & w_cut & estimator_range
    legacy_selected = pretag & tags
    physical_selected = legacy_selected & (sub["reco_pT"] >= 0)
    stages = {
        "truth_fiducial_legacy_finite_filter": np.ones(len(train), dtype=bool),
        "proton_tag_alone": tags,
        "reco_mass": mass_cut,
        "reco_mass_then_epz": mass_cut & epz_cut,
        "reco_mass_then_epz_then_W": mass_cut & epz_cut & w_cut,
        "all_other_cuts_before_proton_tag": pretag,
        "legacy_final_with_proton_tag": legacy_selected,
        "physical_reco_pT_final_with_proton_tag": physical_selected,
    }
    source_scale = args.projection_lumi / args.source_lumi
    projection_scale = args.projection_lumi / (args.source_lumi * (1. - args.train_fraction))
    split_scales = {"all": source_scale, "train": source_scale / args.train_fraction,
                    "test": projection_scale}
    cutflow = []
    for i, (low, high) in enumerate(zip(EDGES[:-1], EDGES[1:])):
        in_bin = (sub["true_t"] >= low) & ((sub["true_t"] <= high) if i == len(EDGES)-2 else (sub["true_t"] < high))
        for split, subset in {"all": np.ones(len(train), dtype=bool), "train": train, "test": test}.items():
            for stage, mask in stages.items():
                row = {"truth_bin": i, "truth_low": low, "truth_high": high, "split": split, "stage": stage}
                row.update(counts(in_bin & subset & mask, sub["weight"], split_scales[split]))
                cutflow.append(row)
    responses = {}
    for name, selected in {"legacy_tagged": legacy_selected, "without_proton_tag": pretag,
                           "physical_reco_pT_tagged": physical_selected}.items():
        responses[name] = response(sub["true_t"], reco_t, sub["weight"], selected, train, test, EDGES, projection_scale)
    diagnoses = []
    for i in responses["legacy_tagged"]["unsupported_truth_bins"]:
        low, high = EDGES[i:i+2]
        in_bin = (sub["true_t"] >= low) & ((sub["true_t"] <= high) if i == len(EDGES)-2 else (sub["true_t"] < high))
        truth_positive = int(np.count_nonzero(in_bin & train & (sub["weight"] > 0)))
        before = int(np.count_nonzero(in_bin & train & pretag & (sub["weight"] > 0)))
        if truth_positive == 0:
            reason = "No positive-weight training truth in this bin; held-out truth is positive."
        elif before == 0:
            reason = "No positive-weight training event survives the non-proton reconstruction cuts."
        else:
            reason = "The proton requirement removes all positive-weight training events surviving the other cuts."
        diagnoses.append({"truth_bin": i, "observed_finite_MC_reason": reason})
    # Outside-fit feed-in at 0<=reco pT2<=1, including true t>5 omitted by the legacy unfolding.
    with np.errstate(over="ignore", invalid="ignore"):
        reco_squared_all = arr["reco_pT"] ** 2
    reco_core = (global_finite & (arr["reco_pT"] >= 0) & (reco_squared_all <= 1.)
                 & (arr["reco_mass"] >= mass - 3. * props["sigma"])
                 & (arr["reco_mass"] <= mass + 3. * props["sigma"])
                 & (np.abs(arr["sigma_epz"] - 18.) <= 1.8)
                 & (arr["reco_W"] >= 20.) & (arr["reco_W"] <= 70.))
    all_tags = policy_module.proton_tag_mask(arr)
    tail_rows = []
    for region, region_mask in {"truth_0_to_1_inclusive": (arr["true_t"] >= 0) & (arr["true_t"] <= 1),
                                "truth_above_1_through_5": (arr["true_t"] > 1) & (arr["true_t"] <= 5),
                                "truth_above_5": arr["true_t"] > 5}.items():
        for selection, mask in {"truth_before_legacy_global_finite_filter": truth_raw,
                                "truth_after_legacy_global_finite_filter": truth_legacy,
                                "reco_core_without_proton_tag": truth_raw & reco_core,
                                "reco_core_with_proton_tag": truth_raw & reco_core & all_tags}.items():
            row = {"region": region, "selection": selection}
            row.update(counts(region_mask & mask, weight, source_scale))
            tail_rows.append(row)
    missing_truth = np.logical_and.reduce([arr[name] == -999. for name in ("true_mass", "true_t", "true_W", "true_Q2")])
    expected_sum = args.source_lumi * props["visible_xsec_fb"]
    normalization = {"source_lumi_invfb": args.source_lumi, "target_lumi_invfb": args.projection_lumi,
                     "assigned_visible_xsec_fb": props["visible_xsec_fb"], "stored_weight_sum": float(weight.sum()),
                     "assigned_expected_source_sum": expected_sum,
                     "weight_sum_matches_assigned": bool(np.isclose(weight.sum(), expected_sum, rtol=1e-6, atol=1e-8)),
                     "weights_changed": False}
    return {
        "schema_version": 1, "kind": "9x130_legacy_true_t_support_audit", "state": args.state,
        "electron_energy": 9., "proton_energy": 130., "truth_edges": EDGES.tolist(),
        "truth_variable": "positive stored true_t", "reconstructed_estimator": "reco_pT squared",
        "w_range": [20., 70.], "q2_max": props["q2_max"], "mass": mass,
        "reco_mass_sigma": props["sigma"], "proton_selection": policy,
        "seed": args.seed, "train_fraction": args.train_fraction,
        "normalization": normalization, "legacy_test_projection_scale": projection_scale,
        "input_rows": n, "legacy_global_finite_rows": int(global_finite.sum()),
        "rows_dropped_by_legacy_global_finite_filter": int((~global_finite).sum()),
        "truth_fiducial_rows_dropped_by_legacy_global_finite_filter": int((truth_raw & ~global_finite).sum()),
        "coherent_missing_truth_rows": int(missing_truth.sum()),
        "legacy_truth_range_rows": len(train), "train_rows": int(train.sum()), "test_rows": int(test.sum()),
        "negative_reco_pT_rows_in_legacy_truth": int((sub["reco_pT"] < 0).sum()),
        "negative_reco_pT_rows_passing_legacy_final": int((legacy_selected & (sub["reco_pT"] < 0)).sum()),
        "sentinel_reco_pT_minus999_rows": int((arr["reco_pT"] == -999.).sum()),
        "responses": responses, "diagnoses": diagnoses, "cutflow": cutflow, "tail_feed_in": tail_rows,
        "scope": [
            "Exact legacy finite-row filtering and RNG split are retained for failure reconstruction.",
            "The physical reco-pT comparison uses that same split; it does not redefine the truth denominator.",
            "Zero support in finite MC does not establish that the physical acceptance is exactly zero.",
            "Removing a proton requirement here is a diagnostic counterfactual, not a replacement result.",
            "The target-luminosity scaling cannot change response support.",
            "This audit does not select new bins, validate an unfolding, or change ROOT/production outputs.",
        ],
        "publication_ready": False,
    }


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--framework-root", type=Path, default=Path.cwd())
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--tree", default="UpsilonTree")
    p.add_argument("--state", choices=tuple(STATES), default="1S")
    p.add_argument("--projection-lumi", type=float, default=1.)
    p.add_argument("--source-lumi", type=float, default=10.)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--train-fraction", type=float, default=.5)
    p.add_argument("--out-dir", type=Path, required=True, help="A fresh directory; never overwrite an existing audit")
    return p


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(argv=None):
    p = parser()
    args = p.parse_args(argv)
    if not np.isfinite(args.projection_lumi) or args.projection_lumi <= 0 or not np.isfinite(args.source_lumi) or args.source_lumi <= 0:
        p.error("Luminosities must be finite and positive")
    if not 0 < args.train_fraction < 1:
        p.error("--train-fraction must lie strictly between zero and one")
    args.input = args.input.expanduser().resolve()
    args.out_dir = args.out_dir.expanduser().resolve()
    if args.out_dir.exists():
        p.error("Output directory already exists; use a fresh --out-dir. Existing audit unchanged.")
    root, policy_module, policy = helpers(args.framework_root, args.out_dir)
    with uproot.open(args.input) as source:
        tree = source[args.tree]
        missing = sorted(set(BASE_BRANCHES + TAG_BRANCHES) - set(tree.keys()))
        if missing:
            raise ValueError("Required branches missing: " + ", ".join(missing))
        arrays = {name: tree[name].array(library="np") for name in BASE_BRANCHES + TAG_BRANCHES}
    result = audit_arrays(arrays, args, policy_module, policy)
    stat = args.input.stat()
    result["input"] = {"path": str(args.input), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns,
                       "tree": args.tree}
    result["source_sha256"] = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in
                                (Path(__file__), root / "workflow/proton_selection.py",
                                 root / "scripts/unfolding/unfold_t_spectrum.py")}
    args.out_dir.mkdir(parents=True, exist_ok=False)
    (args.out_dir / "audit.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    write_csv(args.out_dir / "truth_bin_cutflow.csv", result["cutflow"])
    write_csv(args.out_dir / "tail_feed_in.csv", result["tail_feed_in"])
    summary = []
    for label, response_result in result["responses"].items():
        for i in range(len(EDGES)-1):
            summary.append({"selection": label, "truth_bin": i, "truth_low": EDGES[i], "truth_high": EDGES[i+1],
                            "training_truth_source_sumw": response_result["truth_train_source_sumw"][i],
                            "heldout_truth_target_yield": response_result["truth_test_target_yield"][i],
                            "efficiency": response_result["efficiency"][i],
                            "unsupported_truth": i in response_result["unsupported_truth_bins"]})
    write_csv(args.out_dir / "response_support.csv", summary)
    for label, response_result in result["responses"].items():
        print(label, "unsupported truth", response_result["unsupported_truth_bins"],
              "unsupported reco", response_result["unsupported_reco_bins"])
    for row in result["diagnoses"]:
        print("TRUTH BIN", row["truth_bin"], row["observed_finite_MC_reason"])
    print("Audit:", args.out_dir / "audit.json")
    print("Target luminosity:", args.projection_lumi, "fb^-1; source weights unchanged")
    if not result["normalization"]["weight_sum_matches_assigned"]:
        print("WARN stored sum(weight) differs from the assigned 9x130 rate times the declared source luminosity")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
