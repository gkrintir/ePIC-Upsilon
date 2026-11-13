# ePIC Υ dielectron projections after ECAL rematching

The six merged tuples cover Υ(1S), Υ(2S) and Υ(3S) at 9 × 130 and 9 × 275 GeV. They use the downstream ECAL rematching for EICrecon issue 3011 and full fEMC. The high-energy 1S file is from fresh-01. The other five files are from rest-01. The original projection and the corrected projection use the same physics selections, visible cross sections and beam-specific kinematic windows. The numerical tables distinguish expected selected events from yields extracted in the mass fit.

## Expected proton-tagged signal

| Beam | Luminosity | Matching | Υ(1S) | Υ(2S) | Υ(3S) | Selected sum |
|---|---:|---|---:|---:|---:|---:|
| 9 × 275 GeV | 2.5 fb⁻¹ | Stored | 66.72 | 21.00 | 17.60 | 105.32 |
| 9 × 275 GeV | 2.5 fb⁻¹ | Rematched | 91.56 | 27.34 | 22.69 | 141.58 |
| 9 × 130 GeV | 10 fb⁻¹ | Stored | 107.96 | 33.70 | 27.70 | 169.36 |
| 9 × 130 GeV | 10 fb⁻¹ | Rematched | 148.56 | 43.79 | 35.74 | 228.08 |

The normalizing weight is the visible dielectron cross section times luminosity divided by the actual number of generated rows in the file. The mass, Q² and W truth cuts determine the denominator for the reported selection efficiency. W spans 20–70 GeV at 9 × 130 and 20–100 GeV at 9 × 275, in truth and reconstruction. The 1S truth selection uses Q² up to 1 GeV². The 2S and 3S selections use Q² up to 0.01 GeV². The reconstruction requires a dielectron mass of 8–12 GeV, |Σ(E−p_z)−18| ≤ 1.8 GeV and a positive-momentum reconstructed recoil proton satisfying the framework topology requirement. The summary tables also include the selection with the proton optional.

At 9 × 275 and 2.5 fb⁻¹ the simultaneous fit yields 91.79, 27.08 and 22.71 signal events with core widths of 54.44, 60.26 and 63.20 MeV. The original simultaneous fit gives 66.91, 20.78 and 17.64 events. All mass-fit uncertainty intervals are signal-only profile estimates. Backgrounds and the observed pure-template shape mismatch lie outside these intervals. Low-count coverage remains unvalidated.

## Numerical checks

* For 1S at 9 × 275 the old and new samples have the same 1,137,600 truth records, including multiplicities. Among 1,126,936 eligible truth events the tagged count rises from 144,473 to 198,265. Two of the other five sample pairs have different generated populations. The per-file normalization and a separate comparison of common truth records account for these differences. See audit/selection_audit_README.md.
* At the earlier opposite-charge candidate stage the expected 1S yield at 9 × 275 and 100 fb⁻¹ rises from 9,183 to 12,646, 26.9% below the ATHENA tracking-weighted estimate of 17,290 on this beam configuration. The published ATHENA Fig. 26 mass spectrum uses 18 × 275. The comparison at the candidate stage precedes scattered-electron, mass and proton selections.
* The mass spectra and pure-state calibrations entering the simultaneous fit are checked against each of the six native state runs. The 528 displayed cells of the 1S efficiency maps were checked directly against the corrected ROOT tree. The independent map check also records overflow.
* The fixed-bin state analyses each complete 26 tasks for 1S. Each excited-state run completes 24 of 26. A truth-populated highest-W bin has zero reconstruction response for 2S and 3S at both energies. The corresponding fixed-bin W extraction and dependent plot are withheld. One display-only correction rounds an error-bar length of −3.39 × 10⁻²¹ to zero. The counts and statistical limits retain their calculated values.
* The separate adaptive optimization examined 24 state/energy/observable combinations. Four passed conditional validation, 14 failed one or more checks, and six had no feasible slope binning. The four conditional passes are the three pT² slopes at 9 × 130 and the 1S pT² slope at 9 × 275. For the latter, b = 2.46 with a projected interval of 1.93–3.01 GeV⁻². The W and rapidity projections are diagnostics. Adaptive W boundaries can differ from the historical display. Use the saved bin edges before comparing rates.
* For the trace-anomaly projection only the Υ(1S) 22 ≤ W < 28 GeV window enters the forward-intercept sensitivity. The left panel reproduces published ATHENA Fig. 40 right by digitization. The ePIC mass-fraction result assumes a fixed polarizability and the threshold relation f = 0.224√A while profiling the slope and migration normalizations. At 9 × 275 and 100 fb⁻¹ the rematched conditional interval spans 0.114–0.312. At 2.5 fb⁻¹ a finite two-sided interval cannot be obtained at that energy. Both energies and luminosities are in trace275/fig40/both_energies_forward_forecasts.csv.

## Files

* ePIC_Upsilon_rematched_figures.pdf contains the main figures and a one-page summary.
* audit/ contains six-sample comparisons, paired truth tests, the Fig. 26 window check, binned spectra and the associated scripts.
* massfits/ contains the simultaneous-fit spectra, calibration inputs, per-state central values and profiled intervals for both versions and energies.
* framework/ contains the six state manifests, validations, selected counts, fixed-bin outputs, reoptimization tables and conditionally validated plots. The status of each observable appears in native_reoptimization_status.csv. Absolute source paths in manifests identify the run as executed.
* paper_figures/ contains the figures, plotted bins and audits for the corrected paper setting. The diagnostic W and differential plots are under spectra/. The displayed bins combine selected native bins as recorded in the per-figure audit.
* trace130/ and trace275/ contain the near-threshold studies, original reference tables and regenerated baseline checks.
* sources/verified_inputs.json identifies each merged ROOT file by URL, size, Adler-32 and SHA-256. The archive omits the ROOT files. sources/fetch.py stages the six downloads using their CERNBox metadata. The associated analysis repository is https://github.com/gkrintir/ePIC-Upsilon at 072c1fdc1c642a5cf21242b0363f24bcd89f84d5. The original input definitions are in the repository's inputs/paper_inputs.csv.
