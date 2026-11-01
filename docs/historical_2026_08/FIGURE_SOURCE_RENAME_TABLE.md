# Report figure to renamed source table

The reports do not use consistent figure numbers, so each figure is identified
by report and PDF page. `C#` means a zero-indexed notebook code cell. Original
files are preserved under `archive/original_code/`; renamed files are the
canonical framework paths.

Coverage labels:

- **Tuple-direct**: the retained code reads `UpsilonTree` and fills the plotted
  data itself.
- **Tuple-reconstructed**: the recovered selection, branches, axes, and fills
  reproduce a stored histogram bin-for-bin, including flow bins.
- **Pre-binned**: the code renders an already-filled ROOT histogram.
- **Partial**: relevant ingredients survive but the exact visible panel is not
  explicitly reproduced.
- **External/manual/missing**: no local producer was identified.

None of the retained or promoted code calls or loads `DVMP_JPsi_Analysis.C`.
The uploaded tuples and encoded selections are the runtime inputs. See
`UPLOADED_ROOT_CROSSCHECK.md` for the six-file validation.

## Original report

| Figure | Old script/notebook | New script/notebook | Coverage after ROOT cross-check |
|---|---|---|---|
| p. 2, DSCB example | - | - | External figure |
| p. 4, generated Upsilon mass | `MassFit-Copy4.ipynb` C7; `MassFit-Copy8.ipynb` C0 | `notebooks/03_event_selection_truth_and_bt_studies.ipynb` C7; `scripts/plotting/compare_report_campaign_histograms.py` | Tuple-direct; historical campaign TH1 inputs identified |
| p. 5, generated Upsilon rapidity and pT | `MassFit-Copy4.ipynb` C7; `MassFit-Copy8.ipynb` C0 | notebook 03 C7; campaign comparison script | Tuple-direct; historical campaign TH1 inputs identified |
| p. 6, W vs Upsilon rapidity | `MassFit-Copy4.ipynb` C7 | `scripts/validation/reconstruct_report_histograms_from_tuple.py`; `scripts/plotting/compare_report_campaign_histograms.py` | Tuple-reconstructed exactly as `GenLevel/h2_Gen_W_vs_Rapidity` |
| p. 6, generated \|t\| vs W | `MassFit-Copy4.ipynb` C7 | reconstruction script; campaign comparison script; notebook 03 C7 | Tuple-reconstructed exactly as `GenLevel/h2_Gen_t_vs_W` |
| p. 7, lepton pT correlation and eta-vs-pT | `MassFit-Copy4.ipynb` C7 | reconstruction script; campaign comparison script; notebook 03 C7 | Three generator-level TH2 fills recovered exactly |
| p. 8, scattered-electron energy and eta | `MassFit-Copy4.ipynb` C7; `MassFit-Copy8.ipynb` C0 | notebook 03 C7; campaign comparison script | Tuple-direct; archived campaign PNG inputs identified |
| p. 9, truth/reco lepton maps and efficiency | `MassFit-Copy4.ipynb` C7; `MassFit-Copy2.ipynb` C2/C5 | reconstruction script; campaign comparison script | Truth and reco maps tuple-reconstructed exactly; binwise efficiency is derived with an empty-bin caveat |
| p. 10, generated vs reconstructed scattered-electron eta | `MassFit-Copy4.ipynb` C7 and C0/C1 | notebook 03 C7 and C0/C1 | Tuple-direct ingredients; exact overlay remains absent |
| p. 11, electron E/p | `MassFit-Copy8.ipynb` C0 | campaign comparison script | Pre-binned renderer; archived PNG reproduced exactly |
| p. 12, reconstructed mass, W, and Upsilon rapidity | `MassFit-Copy4.ipynb` C0/C1 | notebook 03 C0/C1 | Tuple-direct producer; uploaded tuple reproduces dashboard exactly |
| p. 12, generated/reconstructed W and rapidity overlays | `MassFit-Copy4.ipynb` C7 and C0/C1 | notebook 03 C7 and C0/C1 | Tuple-direct ingredients; exact overlay layout remains absent |
| p. 12, W-vs-rapidity inset | `MassFit-Copy4.ipynb` C7 | reconstruction script; campaign comparison script | Tuple-reconstructed exactly as `GenLevel/h2_Gen_W_vs_Rapidity` |
| p. 13, mass resolution by rapidity | `MassFit-Copy8.ipynb` C0 | campaign comparison script | Pre-binned; all nonempty archived campaign PNGs reproduced exactly |
| p. 14, inclusive fit without electron ID | `MassFit.ipynb` C1, closest analogue | `notebooks/01_mass_fits_inclusive_and_rapidity.ipynb` C1 | Uploaded low tuple reproduces image; fit is numerically invalid |
| p. 14, inclusive fit with electron ID | `MassFit.ipynb` C0 | notebook 01 C0 | Supplied campaign histograms do not reproduce archived fit; exact historical mass input remains missing |
| p. 17, rapidity-bin mass fits | `MassFit.ipynb` C2 | notebook 01 C2 | Uploaded updated high tuple reproduces image; fit covariances are invalid |

## Follow-up report

| Figure | Old script/notebook | New script/notebook | Coverage after ROOT cross-check |
|---|---|---|---|
| p. 3, low-Q2 tagger tracks | `MassFit-Copy4.ipynb` C0/C1; `MassFit-Copy8.ipynb` C0 | notebook 03 C0/C1; campaign comparison script | Tuple-direct distribution; archived campaign PNG reproduced exactly |
| p. 4, central-vs-tagger tracks | `plot_central_tagger_correlation.py` | `scripts/validation/reconstruct_report_histograms_from_tuple.py`; `scripts/plotting/compare_report_campaign_histograms.py`; `scripts/validation/plot_track_multiplicity_correlation.py` | Raw all-event TH2 tuple-reconstructed exactly; selected 1S views also rerun |
| p. 5, central tracks | `MassFit-Copy4.ipynb` C0/C1; `MassFit-Copy8.ipynb` C0 | notebook 03 C0/C1; campaign comparison script | Tuple-direct distribution; archived campaign PNG reproduced exactly |
| p. 6, tagged-proton momentum | `MassFit-Copy8.ipynb` C0 | campaign comparison script | Pre-binned; archived PNG reproduced exactly |
| p. 7, reconstructed \|t\| for 80<W<100 GeV | `MassFit-Copy8.ipynb` C0; `compare_reco_t_distributions_two_energies_swan.py` | campaign comparison script; `scripts/plotting/plot_reco_t_estimators_by_energy.py` | Exact stored-key renderer plus pixel-exact tuple-script rerun under historical settings |
| p. 8, electron isolation fraction | `MassFit-Copy8.ipynb` C0 | campaign comparison script | Pre-binned; archived PNG reproduced exactly |
| p. 9, generated W in acceptance | `MassFit-Copy8.ipynb` C0; `MassFit-Copy4.ipynb` C7 | campaign comparison script; notebook 03 C7 | Accepted-W TH1 survives; exact event-level acceptance expression remains unresolved |

## Performance report

| Figure | Old script/notebook | New script/notebook | Coverage after ROOT cross-check |
|---|---|---|---|
| p. 2, mass-resolution campaign comparison | `MassFit-Copy8.ipynb` C0 | campaign comparison script | Exact inputs found; all archived nonempty TH1 PNGs reproduce byte-for-byte |
| p. 3, track campaign comparison | `MassFit-Copy8.ipynb` C0 | campaign comparison script | Exact inputs found; archived PNGs reproduce byte-for-byte |
| p. 4, eSTARlight vs afterburner | - | - | No retained producer |
| p. 5, cross-section and branching-ratio table | - | - | Manual table; no generator found |
| p. 6, event-selection cutflow | `MassFit-Copy4.ipynb` C8 | notebook 03 C8 | Reads stored `h_CutFlow`; tuple files are identified |
| p. 7, four 1S efficiency panels | `MassFit-Copy2.ipynb` C0 | `notebooks/02_efficiency_and_two_energy_studies.ipynb` C0 | Legacy high tuple reproduces all state outputs byte-for-byte |
| p. 8, scaled-yield dashboard | `MassFit-Copy4.ipynb` C0; revised C1 | notebook 03 C0/C1 | Both legacy and updated tuple inputs identified; retained outputs reproduce exactly |
| p. 9, ZEUS spectrum | - | - | External figure |
| p. 10, Upsilon-family gamma-p cross sections | `MassFit-Copy7.ipynb` C0 | `notebooks/06_gamma_p_cross_section_family.ipynb` C0 | Legacy high tuple reproduces retained output byte-for-byte |
