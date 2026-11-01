"""Four standalone maps of the unchanged final-selection efficiency.

No projected-yield scaling, interpolation, smoothing, or interval replacement is
performed here. Empty MC cells remain undefined; populated zero-efficiency cells
remain at zero on the common color scale. Numerical intervals live in the CSV.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np


MAPS = {
    "pt_q2": (r"$p_{T,\Upsilon}^{\mathrm{true}}$ [GeV]",
               r"$Q^{2}_{\mathrm{true}}$ [GeV$^2$]", False, True),
    "pt_y": (r"$p_{T,\Upsilon}^{\mathrm{true}}$ [GeV]",
              r"$y_{\Upsilon}^{\mathrm{true}}$", False, False),
    "xv_q2": (r"$x_V^{\mathrm{true}}$",
               r"$Q^{2}_{\mathrm{true}}$ [GeV$^2$]", True, True),
    "xb_q2": (r"$x_B^{\mathrm{true}}$",
               r"$Q^{2}_{\mathrm{true}}$ [GeV$^2$]", True, True),
}


def _grid(rows):
    """Validate a complete rectangular count grid without imputing cells."""
    if not rows:
        raise ValueError("No numerical rows supplied for final efficiency map")
    x = np.array(sorted({float(r[k]) for r in rows for k in ("x_low", "x_high")}))
    y = np.array(sorted({float(r[k]) for r in rows for k in ("y_low", "y_high")}))
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise ValueError("Nonfinite map edges")
    n = np.zeros((len(y) - 1, len(x) - 1), dtype=np.int64)
    k = np.zeros_like(n)
    seen = np.zeros_like(n, dtype=bool)
    for row in rows:
        ix = int(np.searchsorted(x, float(row["x_low"])))
        iy = int(np.searchsorted(y, float(row["y_low"])))
        if (ix >= len(x) - 1 or iy >= len(y) - 1
                or x[ix + 1] != float(row["x_high"])
                or y[iy + 1] != float(row["y_high"])):
            raise ValueError("Map cells must tile adjacent edges")
        if seen[iy, ix]:
            raise ValueError("Duplicate map cell")
        ni, ki = float(row["denominator"]), float(row["selected"])
        if not (np.isfinite(ni) and np.isfinite(ki)
                and ni.is_integer() and ki.is_integer() and 0 <= ki <= ni):
            raise ValueError("Map counts must be integers with 0 <= selected <= denominator")
        if ni > 0:
            stored = row.get("efficiency")
            if stored is None or not np.isfinite(float(stored)) or not np.isclose(
                    float(stored), ki / ni, rtol=0, atol=1e-12):
                raise ValueError("Stored efficiency disagrees with map counts")
        n[iy, ix], k[iy, ix], seen[iy, ix] = int(ni), int(ki), True
    if not np.all(seen):
        raise ValueError("Incomplete map grid; missing cells must be explicitly counted")
    efficiency = np.zeros_like(n, dtype=float)
    np.divide(k, n, out=efficiency, where=n > 0)
    return x, y, n, np.ma.array(100.0 * efficiency, mask=n == 0)


def render_final_maps(rows, sample, settings, output_dir, synthetic=False):
    """Write four single-panel PDF/PNG pairs, refusing existing targets.

    ``output_dir`` is the plots directory for one sample. Returned paths are
    absolute. Input counts always refer to truth coordinates in both numerator
    and denominator. The numerator is the complete final selection.
    """
    output = Path(output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    selected = [r for r in rows if r.get("sample", sample["key"]) == sample["key"]]
    paths = [output / (name + "." + extension)
             for name in MAPS for extension in ("pdf", "png")]
    existing = [str(p) for p in paths if p.exists()]
    if existing:
        raise FileExistsError("Refusing to overwrite final maps: " + ", ".join(existing))
    # Validate every map before writing any figure.
    grids = {name: _grid([r for r in selected if r["map"] == name]) for name in MAPS}
    for name, (_, _, logx, logy) in MAPS.items():
        x, y, _, _ = grids[name]
        if (logx and np.any(x <= 0)) or (logy and np.any(y <= 0)):
            raise ValueError("Logarithmic map axes require positive edges: " + name)
    style = {
        "font.family": "sans-serif", "font.size": 12,
        "axes.labelsize": 15, "xtick.labelsize": 12, "ytick.labelsize": 12,
        "axes.linewidth": 1.0, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 180,
    }
    with plt.rc_context(style):
        for name, (xlabel, ylabel, logx, logy) in MAPS.items():
            x, y, n, efficiency = grids[name]
            fig, ax = plt.subplots(figsize=(7.4, 5.8))
            fig.subplots_adjust(left=0.135, right=0.825, bottom=0.14, top=0.835)
            cmap = plt.get_cmap("viridis").copy()
            cmap.set_bad("#eeeeee")
            mesh = ax.pcolormesh(x, y, efficiency, cmap=cmap,
                                 norm=Normalize(0, 100), shading="flat", rasterized=True)
            if logx:
                ax.set_xscale("log")
            if logy:
                ax.set_yscale("log")
            ax.set_xlim(x[0], x[-1])
            ax.set_ylim(y[0], y[-1])
            ax.set_xlabel(xlabel, labelpad=8)
            ax.set_ylabel(ylabel, labelpad=8)
            ax.tick_params(which="both", top=True, right=True, length=5)
            ax.tick_params(which="minor", length=2.5)
            color_ax = fig.add_axes([0.85, 0.14, 0.028, 0.695])
            colorbar = fig.colorbar(mesh, cax=color_ax, ticks=np.arange(0, 101, 20))
            colorbar.set_label("Efficiency [%]", fontsize=14, labelpad=9)
            colorbar.ax.tick_params(labelsize=11, direction="in")
            fig.text(0.135, 0.955, "ePIC Simulation", fontsize=16, weight="bold", va="top")
            fig.text(0.825, 0.955, r"$%g\times%g$ GeV" % (
                settings.get("electron_energy", 9), settings.get("proton_energy", 275)),
                fontsize=13, ha="right", va="top")
            fig.text(0.135, 0.895, r"$\Upsilon(%s)\to e^+e^-$" % sample["state"],
                     fontsize=11, va="top")
            fig.text(0.825, 0.895, r"$Q^2\leq %g$ GeV$^2$" % sample["q2_max"],
                     fontsize=12, ha="right", va="top")
            if synthetic:
                ax.text(0.5, 0.5, "SYNTHETIC TEST INPUTS\nNOT PHYSICS RESULTS",
                        transform=ax.transAxes, ha="center", va="center", rotation=24,
                        fontsize=17, color="#c52828", weight="bold", alpha=0.8, zorder=5)
            try:
                for extension in ("pdf", "png"):
                    fig.savefig(output / (name + "." + extension))
            finally:
                plt.close(fig)
    return [str(path) for path in paths]
