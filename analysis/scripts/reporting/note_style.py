"""Shared presentation adapted from the reviewed figures in the supplied note.

The reviewed PNGs, rather than an available style source, define this layout.
New samples keep their actual energy, selections, and normalization in subtitles.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE = "#1E90FF"
ORANGE = "#FF8C00"
PURPLE = "#9467BD"
SIGNAL_BLUE = "#1F77B4"
COLORS = [BLUE, ORANGE, PURPLE]


def set_style():
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 13,
        "axes.labelsize": 15, "axes.titlesize": 15,
        "axes.titleweight": "normal", "axes.linewidth": 0.9,
        "axes.grid": True, "axes.axisbelow": True,
        "grid.color": "#b0b0b0", "grid.alpha": 0.20, "grid.linewidth": 0.8,
        "xtick.labelsize": 12, "ytick.labelsize": 12,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
        "legend.fontsize": 11, "legend.frameon": False,
        "lines.linewidth": 1.8, "lines.markersize": 4,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white", "savefig.dpi": 180,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.formatter.use_mathtext": True,
    })


def finish(fig, outdir, stem, title, subtitle, handles=None, labels=None,
           legend_ncol=3):
    """Reserve note-like header and common-legend space; save PNG and PDF."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if title:
        fig.suptitle(title, x=0.5, y=0.988, fontsize=21, fontweight="bold")
    if subtitle:
        fig.text(0.5, 0.940, subtitle, ha="center", va="top", fontsize=12,
                 linespacing=1.5)
    if handles:
        fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.012),
                   ncol=legend_ncol, frameon=False, fontsize=11,
                   columnspacing=1.4, handlelength=2.7)
    # Use a consistent page structure rather than tight bounding-box cropping.
    fig.tight_layout(rect=(0.02, 0.12 if handles else 0.055, 0.99, 0.90),
                     h_pad=2.4, w_pad=2.0)
    paths = []
    for suffix in ("png", "pdf"):
        path = outdir / (stem + "." + suffix)
        fig.savefig(path, dpi=180)
        paths.append(str(path))
    plt.close(fig)
    return paths
