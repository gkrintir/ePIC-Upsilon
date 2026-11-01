# Known issues in the supplied analysis snapshot

The files under `results/reference/` are retained exactly as evidence. They are
not automatically eligible as report inputs. The quality checker records these
known problems so that a future run can distinguish legacy findings from newly
introduced failures.

## Critical reference-output contradictions

1. The file
   `sigma_gp_1S_10x100_W20to26p5_Q2lt0p001_unfolded_withEffMC.csv` has zero
   reconstructed yield, zero efficiency, blank unfolded cross section, and
   `plotted=False` in every bin. It also begins at W = 10 GeV despite its
   W20 filename.

2. Its normalization table uses `target_visible_xsec_fb=304.793` and infers
   327.584 fb, values characteristic of the high-energy sample. The canonical
   low-energy 1S value used elsewhere is 104.793 fb.

3. `overlay_sigma_gp_1S_10x100_18x275.csv` consequently contains no usable
   low-energy points. Existing selected PNG/PDF figures predate the current
   low-energy CSV and cannot be reproduced from that input snapshot.

4. `diff_gamma_p_family_10x100_W20to26p5_metadata.csv` and
   `diff_gamma_p_1S_10x100_W20to26p5_bSlope_metadata.csv` record W min = 10
   GeV. The former also records zero raw visible cross section for 2S and 3S,
   yielding zero or NaN products.

5. The high-energy common-W eX and pT2 unfolded t files have zero efficiency
   in every high-t bin. The promoted hybrid CSV therefore has invalid
   high-energy entries. The metadata shows only 539 high-energy truth-common
   events before the reconstruction selection.

6. `unfold_t_1S_10x100_eX_zeusBins_normalization.csv` records `method=pT2`, in
   conflict with the filename.

7. A legacy t-method product named W20to26p5 records W max = 126.5 GeV.

## Method and convention risks

- The family and differential scripts labeled `fixed_unfolded` in the original
  archive are same-sample bin-by-bin efficiency corrections. Their unity
  closure is structurally expected and is not an independent unfolding test.
- Standalone t unfolding normalizes after its W and Q2 selections. Several
  other workflows normalize by state before the W selection. The workflow
  makes the intended normalization convention explicit per task.
- Multiple ROOT aliases appear across source and metadata. Six supplied files
  now have verified checksums and explicit roles in
  `config/verified_root_inputs.toml`. The two high-energy tuples share the same
  kinematics but different weights; the updated low-energy tuple is a subset
  of the legacy tuple with different weights. Never substitute an alias
  silently.
- The HERA b-overlay scripts contain approximate hand-entered comparison
  anchors. Only the digitized ZEUS inputs in `data/reference/` have a preserved
  tabular source in this package.
- No supplied output records input checksums or a complete software environment.

## Notebook risks

- Cached mass-fit output contains non-positive-definite covariance warnings,
  call-limit termination, abnormal termination, and negative fitted component
  normalizations. The displayed numerical fit values are unvalidated.
- A two-energy notebook cell reads the low-energy file while labeling or
  normalizing it as 18x275. That cell is excluded from the default pipeline.
- A decay-efficiency cell has mismatched input/label state and cached outputs;
  it is not automated by default.
- One t-migration notebook creates a literal `~/t_binning_study` directory
  because it does not expand the user path. The new runner rejects literal
  tilde output components.
- Some notebook outputs are cached from earlier source or input states. Source
  text and displayed results must be treated separately.

## Report-specific provenance gaps

The original and follow-up report decks state that their reconstruction began
from `DVMP_JPsi_Analysis.C`, but none of the retained or promoted Python code
calls that macro. The unavailable lower layer is the historical conversion
from raw ePIC output to `UpsilonTree`. Once the tuple exists, all ten supplied
TH2 histograms are now reconstructed exactly by encoded NumPy rules. Some TH1
objects still exist only as prefilled histograms, including mass resolution,
electron E/p, electron isolation, and `h_CutFlow`.

The saved Notebook 07 source is not the source state that made the archived
campaign PNGs. It normalizes the histograms, changes labels, ignores TH2, and
overcounts skipped empty plots. The promoted campaign renderer restores the
verified raw-count state.

No current producer was found for `Fit_Inclusive_Upsilon_DSCB.png/.pdf`, the
raw-versus-afterburner eSTARlight dashboard in the Performance report, or the
manual cross-section table. Current notebook code recovers the mass-fit model
and can create newer fit figures, but it does not export the exact original
report filename/layout. These distinctions are explicit in
`REPORT_CODE_MAP.md`, `report_output_map.csv`, and
`UPLOADED_ROOT_CROSSCHECK.md`.
