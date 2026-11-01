# Selections and normalization

This repository organizes existing calculations. It does not silently harmonize
physics selections or turn a failed validation into an accepted measurement.

## Executable Jul26 baseline

| State | Q² maximum [GeV²] | B(ee) | Visible σ at 9×130 [fb] | Visible σ at 9×275 [fb] | Input directory |
|---|---:|---:|---:|---:|---|
| 1S | 1 | 0.0238 | 119.685 | 210.136 | output_newID |
| 2S | 0.01 | 0.0191 | 29.297 | 53.133 | output |
| 3S | 0.01 | 0.0218 | 23.447 | 42.969 | output |

These values are enforced by `analysis/workflow/normalization_contract.py` and
current state configuration. They differ from recent note prose asserting Q²<1
for all states. That prose must not be used to change the executable baseline:
matching generated samples and rates are needed before changing the excited-state
phase space.

The all-generated normalization uses ordinary `weight`. Preparation assigns a
uniform visible dielectron weight σvis × Ltuple / Ngenerated to each generated
row. The original producer's per-job merged weight sum is not treated as the
physical total rate. `weight_proj` is not substituted for this ordinary weight.
Normalized files are new artifacts; raw ROOT files remain unchanged. The photon
flux and branching fraction are applied by extraction according to its declared
cross-section convention.

Tuple normalization at 10 fb⁻¹ and projection luminosity are distinct. A 1 fb⁻¹
projection scales expected yields by 1/10, and a 2.5 fb⁻¹ projection by 2.5/10.
Use configuration and manifests to establish luminosity, never a folder name.
Training/test response splits restore exposure using their configured split
probability; they do not silently halve the intended projection.

## Selection categories

`inclusive` means **no explicit reconstructed recoil-proton requirement**. It
still describes the exclusive-signal simulation and retains the remaining
analysis cuts; it does not denote a different inclusive-production cross section.

`proton-tagged` requires `reco_topology` in {1,2} and a finite positive
`reco_tagged_proton_p`. The requirement acts on reconstructed selection only.
Generated normalization and truth denominators are unchanged. It represents the
producer's RP/off-momentum candidate, not a truth-match or exclusivity guarantee.
Keep categories in separate runtime workspaces and explicit plot profiles.

Current source validation retains coherent missing-truth -999 rows in the
all-generated normalization; they are excluded from truth-fiducial studies
without renormalizing the remaining events. Such producer losses are not covered
by statistical toys. Invalid reconstruction is classified as inefficiency.

## Observables and validation

- `pt2`: true meson pT² unfolded/fitted using reconstructed pT².
- `abs_t_via_pt2`: true positive |t| using reconstructed meson pT². It is a
  different response and slope definition.
- `w`: flux-bin-averaged σ(γp); flux is integrated separately in each W bin.
- `y`: primary inclusive-Υ ep dσ/dy, corrected for B(ee). The separately labeled
  total-flux-normalized companion is not a local γp cross section in rapidity.

The 9×130 add-on studies truth 0–1 GeV² for the two slope observables, with
reconstructed control up to 5 GeV² and modeled outside-fit feed-in. This does not
recover an unsupported truth tail extending to 5 GeV².

`COMPLETE` is an execution status. `PASS_CONDITIONAL` describes signal-only
conditional tests under the listed assumptions; no result automatically becomes
publication-ready. `NOT_VALIDATED` and `NO_FEASIBLE_SLOPE_BINNING` must remain
visible downstream. Empty result CSVs can be an intentional unavailable outcome.
Never fill them from an older run.

Source provenance, assigned-rate consistency and independent numerical photon
flux integration do not authenticate which executable produced an EOS dataset.
