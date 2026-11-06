"""Shared SWAN I/O, publication style and fixed-correction yield projections.

No ROOT input is modified. Cross sections are already flux/BR corrected by
the framework; only expected reconstructed counts and their limits change
with luminosity. No second luminosity, branching-ratio or flux correction.
"""
from .electron_id import electron_id_metadata
from dataclasses import dataclass, field
from pathlib import Path
import hashlib
import json
import warnings

import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy.stats import gamma, poisson

STATES = ("1S", "2S", "3S")
COLORS = {"1S": "#0072B2", "2S": "#D55E00", "3S": "#009E73"}
MARKERS = {"1S": "o", "2S": "s", "3S": "^"}
DEFAULT_COVERAGE = 0.6826894921370859


@dataclass
class Config:
    runs: dict
    energy: str = "9x275"
    lumi_invfb: float = 2.5
    output_dir: Path = Path("figures/results")
    campaign: str = "26.07.1 Campaign"
    categories: dict = field(default_factory=lambda: {
        "1S": "new electron ID", "2S": "standard electron ID", "3S": "standard electron ID"})
    q2_max: dict = field(default_factory=lambda: {"1S": 1., "2S": .01, "3S": .01})
    coverage: float = DEFAULT_COVERAGE
    synthetic: bool = False
    require_proton_tag: bool = False  # Enabled by both tagged settings constructors.
    allow_diagnostic_results: bool = False  # Explicit opt-in for unvalidated saved inference.
    reoptimization_run: Path = None
    rebin_edges: dict = field(default_factory=dict)

    def __post_init__(self):
        self.runs = {s: Path(p).expanduser() for s, p in self.runs.items()}
        self.output_dir = Path(self.output_dir).expanduser()
        if self.reoptimization_run is not None:
            self.reoptimization_run = Path(self.reoptimization_run).expanduser()
        if set(self.runs) != set(STATES):
            raise ValueError("Set one explicit run directory for each of 1S, 2S and 3S.")
        if self.energy not in ("9x130", "9x275"):
            raise ValueError("Supported current beam keys: 9x130 and 9x275.")
        if not np.isfinite(self.lumi_invfb) or self.lumi_invfb <= 0:
            raise ValueError("Projection luminosity must be positive.")
        if not 0 < self.coverage < 1:
            raise ValueError("Coverage must lie strictly between zero and one.")


def resolve_artifacts(path):
    path = Path(path).expanduser()
    result = path if path.name == "artifacts" else path / "artifacts"
    if not result.is_dir():
        raise FileNotFoundError(f"Framework artifacts not found: {result}. Check RUN_IDS/FRAMEWORK_ROOT in the configuration cell.")
    return result


def read_csv(path, required=(), state=None):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Required framework CSV is missing: {path}. PNGs alone cannot supply the numerical inputs.")
    frame = pd.read_csv(path)
    if state is not None and "state" in frame:
        frame = frame.loc[frame.state.astype(str) == state].copy()
    missing = sorted(set(required) - set(frame.columns))
    if missing or frame.empty:
        raise ValueError(f"{path}: missing columns {missing}, selected rows={len(frame)}.")
    frame.attrs["source"] = str(path)
    frame.attrs["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return frame


def read_json(path):
    return json.loads(Path(path).read_text())


def _json_clean(obj):
    if isinstance(obj, dict):
        return {str(k): _json_clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, np.ndarray)):
        return [_json_clean(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None
    return obj


def write_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_json_clean(obj), indent=2, allow_nan=False) + "\n")
    return str(path)


def source_luminosity(frame, metadata=None):
    """Luminosity of saved weighted yields, never an assumed 10 fb^-1."""
    candidates = []
    for key in ("yield_lumi_invfb", "active_luminosity_invfb"):
        if key in frame:
            candidates.append((key, pd.to_numeric(frame[key], errors="coerce").to_numpy()))
    if not candidates and "projection_lumi_invfb" in frame:
        if "stat_mode" not in frame or set(frame.stat_mode.astype(str)) != {"projected"}:
            raise ValueError("projection_lumi_invfb alone is not sufficient for non-projected source yields.")
        candidates.append(("projection_lumi_invfb", pd.to_numeric(frame.projection_lumi_invfb).to_numpy()))
    if not candidates and metadata:
        stat = metadata.get("statistical_configuration", metadata)
        for key in ("active_luminosity_invfb", "yield_lumi_invfb"):
            if key in stat:
                candidates.append((key, np.atleast_1d(float(stat[key]))))
                break
    if not candidates:
        raise ValueError("Source yield luminosity is absent. Regenerate metadata; do not guess it from a run name.")
    values = np.concatenate([np.atleast_1d(v) for _, v in candidates])
    if not np.all(np.isfinite(values)) or np.any(values <= 0) or not np.allclose(values, values[0], rtol=1e-8, atol=0):
        raise ValueError("Saved source luminosities are invalid or inconsistent.")
    return float(values[0])


def expected_intervals(mean, coverage=DEFAULT_COVERAGE):
    """Framework median-Garwood forecast; the mean itself is not rounded.

    These are ensemble-median endpoints for a future integer Poisson count,
    not confidence limits computed by treating a fractional mean as data.
    """
    mu = np.asarray(mean, dtype=float)
    if np.any(~np.isfinite(mu)) or np.any(mu < 0) or not 0 < coverage < 1:
        raise ValueError("Finite nonnegative expected counts and 0<coverage<1 are required.")
    n = poisson.ppf(.5, mu)
    tail = (1 - coverage) / 2
    low = np.zeros_like(mu)
    positive = n > 0
    low[positive] = gamma.ppf(tail, n[positive])
    high = gamma.isf(tail, n + 1)
    return {"lower": low, "upper": high, "reference_count": n,
            "method": "median_garwood", "coverage": coverage}


def boolean(values):
    return pd.Series(values).astype(str).str.lower().isin(("1", "true", "1.0")).to_numpy()


def project_cross_section(frame, target_lumi, central_column,
                          count_column="n_reco_weighted_visible", scale_column=None,
                          coverage=DEFAULT_COVERAGE):
    """Recompute count limits with fixed corrections, retaining central sigma.

    Excludes unsupported or upstream-masked bins. A supplied finite positive
    conversion factor permits a genuine supported zero-yield bin; it is never
    invented from 0/0. This function is for bin-by-bin tables, not unfolding.
    """
    f = frame.copy()
    if not np.isfinite(target_lumi) or target_lumi <= 0:
        raise ValueError("Target luminosity must be finite and positive.")
    if "unfolding_identifiable" in f:
        raise ValueError("Use bin-by-bin tables here. Unfolded covariance/toys require a new projection run.")
    if "output_convention" in f and set(f.output_convention.astype(str)) != {"inclusive_upsilon"}:
        raise ValueError("Expected inclusive_upsilon cross sections, already divided by branching fraction.")
    source = source_luminosity(f)
    scale = target_lumi / source
    y = pd.to_numeric(f[central_column], errors="coerce").to_numpy(float)
    counts = pd.to_numeric(f[count_column], errors="coerce").to_numpy(float)
    if np.any(np.isfinite(counts) & (counts < 0)):
        raise ValueError("Negative expected signal counts cannot use this Poisson projection.")
    if scale_column is None:
        known = {"sigma_gamma_p_pb": "reco_count_to_sigma_pb",
                 "dsigma_dx_pb": "reco_count_to_dsigma_dx_pb"}
        scale_column = known.get(central_column)
    if scale_column and scale_column in f:
        source_factor = pd.to_numeric(f[scale_column], errors="coerce").to_numpy(float)
    else:
        source_factor = np.divide(y, counts, out=np.full(len(f), np.nan), where=counts > 0)
    supported = np.isfinite(source_factor) & (source_factor > 0) & np.isfinite(y) & (y >= 0) & np.isfinite(counts)
    if not np.allclose(y[supported], (counts * source_factor)[supported], rtol=2e-6, atol=1e-9):
        raise ValueError("Saved central cross section does not close with reconstructed yield times correction. Check input normalization.")
    valid = supported.copy()
    for key in ("plotted", "passes_min_eff"):
        if key in f:
            valid &= boolean(f[key])
    mean_count = counts * scale
    factor = source_factor / scale
    low, high, reference = (np.full(len(f), np.nan) for _ in range(3))
    if np.any(supported):
        ci = expected_intervals(mean_count[supported], coverage)
        low[supported] = ci["lower"] * factor[supported]
        high[supported] = ci["upper"] * factor[supported]
        reference[supported] = ci["reference_count"]
    f["mean_count"], f["mean"] = mean_count, y
    f["lower"], f["upper"], f["valid"] = low, high, valid
    f["source_lumi"] = source
    f["target_lumi_invfb"] = target_lumi
    f["count_to_cross_section"] = factor
    f["interval_reference_count"] = reference
    f["projection_variance"] = mean_count * factor**2
    f["projection_interval_method"] = "median_garwood"
    f.attrs.update(frame.attrs)
    f.attrs["projection_assumptions"] = "Signal-only independent reconstructed counts; fixed bin-by-bin corrections. No response-MC, background, luminosity, branching or systematic uncertainty."
    return f


def sample_metadata(config, state):
    base = resolve_artifacts(config.runs[state])
    path = base / "performance" / config.energy / "performance_audit.json"
    audit = read_json(path)
    proton_metadata = {}
    if config.require_proton_tag:
        from .proton_inputs import full_run_policy
        proton_metadata = full_run_policy(config, state, audit)
    if str(audit.get("state", state)) != state:
        raise ValueError(f"{path}: wrong state.")
    if str(audit.get("energy_label", "")).replace(" ", "").replace("GeV", "") not in ("", config.energy):
        raise ValueError(f"{path}: energy label conflicts with requested energy.")
    selection = audit["truth_selection"]
    q2 = float(selection["Q2_GeV2"][1])
    w_range = list(map(float, selection["W_GeV"]))
    if state in config.q2_max and not np.isclose(q2, config.q2_max[state]):
        raise ValueError(f"{state}: actual Q2 maximum {q2:g} disagrees with configured {config.q2_max[state]:g}. Update configuration only after checking the intended run.")
    na_path = base / "normalization_audit" / config.energy / "audit.json"
    flux_status = "not available in performance audit"
    normalization_status, normalization_warnings = "not available", []
    if na_path.is_file():
        na = read_json(na_path)
        normalization_status = str(na.get("status", "unknown")).upper()
        failed = [c for c in na.get("checks", []) if str(c.get("status", "")).upper() in ("FAIL", "ERROR")]
        if normalization_status in ("FAIL", "ERROR") or failed:
            names = [str(c.get("id", c.get("name", c.get("check", c.get("key", "unspecified"))))) for c in failed]
            raise ValueError(f"{state}: upstream normalization audit {na_path} has status {normalization_status}; failed checks: {names}. Resolve that failure before rendering cross sections.")
        normalization_warnings = [c for c in na.get("checks", []) if str(c.get("status", "")).upper() == "WARN"]
        if normalization_status == "WARN" or normalization_warnings:
            warnings.warn(f"{state}: normalization audit reports warnings; retained in figure metadata: {na_path}")
        settings = na.get("settings", {})
        for key in ("q2_cut", "q2_max"):
            if key in settings and not np.isclose(float(settings[key]), q2):
                raise ValueError(f"{state}: normalization audit {key} conflicts with performance Q2 selection.")
        if all(k in settings for k in ("q2_cut", "q2_max")):
            flux_status = "cross-checked against normalization audit settings"
    normalization = base / "normalization" / f"{config.energy}_input_audit.json"
    raw_hash = read_json(normalization).get("sha256") if normalization.is_file() else None
    id_metadata = electron_id_metadata(config.categories[state], raw_hash, state, config.energy)
    return {**proton_metadata, **id_metadata, "state": state, "q2_max": q2,
            "w_range": w_range, "source": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "flux_selection_check": flux_status,
            "normalization_audit_status": normalization_status,
            "normalization_audit_warnings": normalization_warnings,
            "denominator_caveat": audit.get("denominator_caveat", "Eligible stored tuple rows"),
            "source_lumi": source_luminosity(pd.DataFrame(), audit)}


def preflight(config, section="all"):
    """Print all missing required files together; never choose a different run."""
    if config.reoptimization_run is not None:
        from .reoptimization_inputs import preflight_reoptimization
        return preflight_reoptimization(config, section)
    rows, missing = [], []
    for state in STATES:
        base = config.runs[state] if config.runs[state].name == "artifacts" else config.runs[state] / "artifacts"
        perf = base / "performance" / config.energy
        files = [perf / "performance_audit.json"]
        if section in ("all", "reconstruction"):
            files += [perf / n for n in ("mass_bins.csv", "mass_fit.json", "selection_efficiency.csv")]
        if section in ("all", "w"):
            files += [base / "cross_sections" / f"{config.energy}_sigma_gamma_p_{state.lower()}{suffix}.csv" for suffix in ("", "_normalization")]
        if section in ("all", "differential"):
            files += [base / "differential" / f"{config.energy}_{state.lower()}_{state}_{v}.csv" for v in ("pT2", "y")]
            files += [base / "differential" / f"{config.energy}_{state.lower()}_metadata.csv"]
        missing += [str(p) for p in files if not p.is_file()]
        if (perf / "performance_audit.json").is_file():
            rows.append(sample_metadata(config, state))
    if missing:
        raise FileNotFoundError("Missing numerical framework outputs:\n" + "\n".join(missing))
    if len({r["q2_max"] for r in rows}) > 1 or len({r["category"] for r in rows}) > 1:
        warnings.warn("Check the per-state Q2 ranges and sample metadata when comparing the spectra.")
    columns = ["state", "q2_max", "category", "w_range", "source_lumi", "flux_selection_check", "normalization_audit_status"]
    if config.require_proton_tag:
        columns += ["proton_selection_mode"]
    return pd.DataFrame(rows)[columns]


# Source sizes chosen for the note's approximately 0.47 single-column scale.
# Keep explicit text sizes in the renderers tied to these shared roles.
FONT_SIZES = dict(axis=16.0, tick=14.5, legend=14.0, annotation=13.0,
                  panel=15.0, header=16.5, inset_tick=12.5, inset_label=13.5)


def set_style():
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
        "font.size": FONT_SIZES["annotation"], "mathtext.fontset": "dejavusans", "axes.labelsize": FONT_SIZES["axis"],
        "axes.titlesize": FONT_SIZES["panel"], "axes.linewidth": 1., "axes.grid": False,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 6, "ytick.major.size": 6, "xtick.minor.size": 3, "ytick.minor.size": 3,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.labelsize": FONT_SIZES["tick"], "ytick.labelsize": FONT_SIZES["tick"],
        "legend.frameon": False, "legend.fontsize": FONT_SIZES["legend"], "legend.handlelength": 1.6,
        "lines.linewidth": 1.6, "lines.markersize": 4.5, "errorbar.capsize": 2.5,
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 220, "text.usetex": False})


def state_label(state):
    return rf"$\Upsilon({state})$"


def annotate(ax, config, panel=None):
    energy = config.energy.replace("x", r"\times")
    ax.text(.035, .965, "ePIC Performance", transform=ax.transAxes, ha="left", va="top", weight="bold", fontsize=FONT_SIZES["header"])
    ax.text(.035, .865, config.campaign + "\n" + rf"$e+p$, ${energy}$ GeV, $\mathcal{{L}}={config.lumi_invfb:g}$ fb$^{{-1}}$",
            transform=ax.transAxes, ha="left", va="top", fontsize=FONT_SIZES["annotation"], linespacing=1.35)
    if panel:
        ax.text(.975, .965, f"({panel})", transform=ax.transAxes, ha="right", va="top", fontsize=FONT_SIZES["panel"])


def draw_intervals(ax, x, mean, lower, upper, xerr=None, **kwargs):
    """Do not clamp absolute forecast limits to force inclusion of the mean."""
    x, y, low, high = np.broadcast_arrays(*[np.asarray(v, float) for v in (x, mean, lower, upper)])
    if np.any(np.isfinite(low) & np.isfinite(high) & (high < low)):
        raise ValueError("Reversed interval endpoints.")
    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(low) & np.isfinite(high)
    if np.all((low[finite] <= y[finite]) & (y[finite] <= high[finite])):
        return ax.errorbar(x, y, yerr=[y-low, high-y], xerr=xerr, **kwargs)
    artist = ax.errorbar(x, y, xerr=xerr, **kwargs)
    color = kwargs.get("ecolor", kwargs.get("color", "black"))
    width = kwargs.get("elinewidth", kwargs.get("linewidth", 1.))
    ax.vlines(x[finite], low[finite], high[finite], color=color, lw=width, zorder=kwargs.get("zorder", 2))
    for endpoint in (low, high):
        ax.plot(x[finite], endpoint[finite], "_", color=color, ms=2*kwargs.get("capsize", 2.5))
    return artist


def save_figure(fig, config, stem):
    config.output_dir.mkdir(parents=True, exist_ok=True)
    if config.synthetic:
        fig.text(.5, -.025, "SYNTHETIC TEST INPUTS - NOT PHYSICS RESULTS", ha="center", fontsize=8, color="#a12323", weight="bold")
    paths = []
    for ext in ("pdf", "png"):
        path = config.output_dir / f"{stem}.{ext}"
        fig.savefig(path, bbox_inches="tight", pad_inches=.10)
        paths.append(str(path))
    return paths
