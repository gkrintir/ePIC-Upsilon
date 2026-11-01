"""Frozen reference from the supplied All States performance producer.

Used only to verify selection parity, independently of the diagnostic engine.
Source payload: scripts/validation/plot_latest_sample_performance.py
Source SHA256: 584a0f982450f176e3f99d3dc0a5f0fb073c00aa6226104c80fd59b26a1ff1ec
"""
import numpy as np
MASSES = {"1S": 9.46030, "2S": 10.02326, "3S": 10.35520}

def selection_masks(arrays, state, electron_energy, q2_cut, w_min, w_max,
                    true_mass_window=0.2):
    """Select one pure-state tuple while retaining its full reconstructed mass tails."""
    mass = MASSES[state]
    truth = (np.isfinite(arrays["true_mass"])
             & (np.abs(arrays["true_mass"] - mass) <= true_mass_window)
             & np.isfinite(arrays["true_Q2"]) & (arrays["true_Q2"] >= 0)
             & (arrays["true_Q2"] <= q2_cut)
             & np.isfinite(arrays["true_W"]) & (arrays["true_W"] >= w_min)
             & (arrays["true_W"] <= w_max))
    reco = (truth & np.isfinite(arrays["reco_mass"])
            & (arrays["reco_mass"] >= 8) & (arrays["reco_mass"] <= 12)
            & np.isfinite(arrays["reco_W"]) & (arrays["reco_W"] >= w_min)
            & (arrays["reco_W"] <= w_max)
            & np.isfinite(arrays["sigma_epz"])
            & (np.abs(arrays["sigma_epz"] - 2 * electron_energy) <= 0.2 * electron_energy))
    return truth, reco
