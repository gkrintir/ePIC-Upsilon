# Read-only Upsilon normalization audit

This package audits the original merged ROOT **as stored**, compares it with the separate prepared ROOT, and recomputes the full-sample bin-by-bin W cross section for `--state 1S`, `2S` or `3S`. Each ROOT input must contain a single generated state. It writes reports and CSVs only. It does not change weights, create ROOT files, rerun the framework, or run eSTARlight.

The normalized ROOT was created by `prepare_latest_sample.py`: every `weight` is replaced by the assigned visible cross section times tuple luminosity divided by the generated-event count. `weight_proj` is changed to preserve its original ratio to `weight`. The original merged file remains separate. This audit exposes both sets of weights explicitly.

## Install and run one energy

The audit is included in the framework under `tools/normalization_audit`. Python 3.9 or newer works; ROOT, SciPy, pandas and tomllib are not required.

```bash
cd /afs/cern.ch/work/g/gkrintir/private/Ypsilon-Fit-Framework
python3 -m pip install -r tools/normalization_audit/requirements.txt

# Set this to the existing run whose plots you are checking.
RUN275=results/generated/runs/YOUR_9x275_RUN_ID

python3 tools/normalization_audit/audit_upsilon_normalization.py \
  --state 1S \
  --merged /eos/cms/store/group/phys_top/gkrintir/upsilon_data/Upsilon1S/9x275/hiAcc/output/merged_Upsilon_final.root \
  --run-dir "$RUN275" \
  --electron-energy 9 --proton-energy 275 \
  --assigned-visible-xsec-fb 210.136 --lumi-invfb 10 --br-ee 0.0238 \
  --w-edges 20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100 \
  --out-dir results/normalization-audit/9x275 --hash-inputs
```

The current July 1S baseline is **119.685 fb at 9x130 and 210.136 fb at 9x275**, with **B_ee=0.0238**. Both visible rates were reproduced by 100,000-event reference-generator runs using the unchanged published cards. They are external inputs to this audit, not inferred from ROOT weights.

The state defaults and July reference visible rates are:

| State | Mass [GeV] | Reco width [GeV] | B_ee | Event/flux Q² maximum [GeV²] | 9x130 [fb] | 9x275 [fb] |
|---|---:|---:|---:|---:|---:|---:|
| 1S | 9.46030 | 0.053 | 0.0238 | 1 | 119.685 | 210.136 |
| 2S | 10.02326 | 0.056 | 0.0191 | 0.01 | 29.297 | 53.133 |
| 3S | 10.35520 | 0.059 | 0.0218 | 0.01 | 23.447 | 42.969 |

`--state` defaults to `1S` for existing commands. Explicit mass, width, branching fraction or Q² arguments override its defaults. `--assigned-visible-xsec-fb` remains a required external input; selecting a state does not infer its normalization from the ROOT file. The tighter Q² settings for the excited states match their July cards and must also be used in the framework extraction.

For example, to audit an existing 9x275 2S run, use the same command with `--state 2S`, the `Upsilon2S/.../output/merged_Upsilon_final.root` path, its own run directory, and `--assigned-visible-xsec-fb 53.133 --br-ee 0.0191 --q2-cut 0.01 --q2-max 0.01`. Use a distinct output directory for each state and energy. The all-state framework runner supplies these settings explicitly.

Use a **new framework run prepared with these corrected rates and branching fractions**. Auditing an older 9x275 normalized file that still sums to 3275.84 events against 210.136 fb should fail. For historical reproduction, explicitly supply the old assignment and `--br-ee 0.0234`, using the direct CLI. The old 327.584 fb value is associated with 18x275, not the July 9x275 configuration.

Reference revisions are cards `f62965a1612fc734dbab18e141f4809740d6b234` and eSTARlight `939b11a24499398392d959db81c7502aeec91046`. This reference reproduction does not authenticate the generator executable used for existing EOS production.

`--run-dir` looks for these exact files beneath that run:

| Artifact | Relative path for 9x275 |
|---|---|
| Prepared ROOT | `artifacts/inputs/9x275_normalized.root` |
| W CSV | `artifacts/cross_sections/9x275_sigma_gamma_p_<state>.csv`, where `<state>` is `1s`, `2s` or `3s` |
| Preparation audit | `artifacts/normalization/9x275_input_audit.json` |

It does not select a different run automatically. Each state has its own run directory, so the prepared ROOT filename remains unchanged. Missing prepared ROOT is a failed requested comparison. Missing optional CSV/audit is reported. To use a scenario nested under `artifacts/scenarios/keep_yneg/`, supply `--normalized`, `--csv`, and `--normalization-audit` explicitly. Each argument overrides discovery. The original merged path must match that scenario: `output`, `output_newID`, `output_keep_yneg`, `output_keep_xpos`, `output_remove_inner`, or `output_remove_outer` as appropriate.

You can omit `--run-dir` and supply `--normalized /full/path/9x275_normalized.root` directly. If neither is supplied, the script performs a raw-only audit and marks normalized comparison as not checked.

## Run both energies and decompose the difference

```bash
RUN130=results/generated/runs/YOUR_9x130_RUN_ID
RUN275=results/generated/runs/YOUR_9x275_RUN_ID
bash tools/normalization_audit/run_two_energies.sh \
  "$RUN130" "$RUN275" results/normalization-audit output
```

This legacy two-energy wrapper is **1S-only**. Replace the final `output` with `output_newID` only if both runs were made with those inputs. The wrapper continues to the second energy even when the first detects a failed check. It returns nonzero if either audit fails. Use the CLI with the same selected state for both reports to compare 2S or 3S energies; a reference from another state is rejected as a mismatched comparison.

Alternatively, run the CLI for 9x130 using `--proton-energy 130`, `--assigned-visible-xsec-fb 119.685`, and `--w-edges 20,25,30,35,40,45,50,55,60,65,70`. Then add `--reference-report results/normalization-audit/9x130/audit.json` to the 9x275 command. Reference and output reports must have different paths.

For common bins, the uniform-preparation truth ratio is decomposed exactly as:

`sigma_current/sigma_reference = assigned_rate_ratio × generated_truth_fraction_ratio × inverse_flux_ratio × inverse_B_ratio`.

The actual normalized-ROOT ratio is displayed alongside this prediction. Luminosity cancels when applied consistently. Shared W-bin boundaries are required; different final-edge conventions and selection settings are flagged.

## Checks and assumptions

- Branching convention: when the prepared ROOT records `generator_branching_fraction_ee`, compare it with the supplied inclusive-conversion branching fraction and reject a mismatch.
- File identity: optional SHA256 of original and prepared ROOTs; embedded source hash versus supplied original; embedded versus external preparation audit; input size/mtime unchanged.
- Event accounting: generated count from `h_EventCount`, tuple entry count, histogram flow bins, empty samples, weight sums, sumw2, invalid/zero/negative weights and repeated weight populations. Equal counts alone do not prove complete production bookkeeping.
- Pair integrity: every common non-weight branch is compared in row order, including schema type changes, with exact differences and tolerance diagnostics. Missing branches/rows and changes to `weight_proj/weight` are reported. Unsupported field comparisons are explicitly skipped, never silently certified.
- Weight replacement: actual normalized total and every event weight are checked against the supplied sigma × L / N. Raw weight groups are diagnostic; they do not uniquely identify generator jobs.
- Shape accounting: original weights, actual normalized weights, hypothetical uniform replacement, and hypothetical global rescaling are shown separately. Uniform replacement assumes unweighted generator events and all generated rows. Global rescaling preserves original relative weights. Neither hypothesis is selected by fitting the model.
- W extraction: use only `weight`, `true_mass`, `true_W`, `reco_mass`, `reco_W`, `sigma_epz`, `true_Q2`. `weight_proj` is not used. Missing `true_Q2` blocks the physics calculation. Both truth and reconstructed selections use the stored truth Q2 cut, matching the current extractor.
- Exact default cuts: selected state mass ±0.2 GeV at truth; selected state mass ±3 times its reco width; reconstructed `sigma_epz` within 2Ee ±0.2Ee; 0≤true_Q2≤the state-specific maximum in the table above; supplied W range. Final histogram upper edge is included. These values are configurable and printed. A sample with no events in the requested truth-mass window fails the state-population check.
- The reconstructed/truth correction is R/T and may exceed one with migration. For supported bins, R/(R/T)=T. This same-sample closure is algebraic and cannot validate the absolute rate. Unsupported bins remain undefined; they are not silently assigned zero cross section.
- Flux: reproduce the framework 500-point W trapezoid, independently integrate differential EPA in log(y) and log(Q2), and expose every bin flux. Beam hypotheses are compared with stored truth scattered-electron kinematics where available. This tests stored reconstruction formulas, not the production card's provenance.
- Saved CSV: select rows for the requested state, then validate central truth/reco yields, correction, flux, cross section, optional truth closure, plotted mask and luminosity metadata. A legacy CSV without a `state` column is accepted only for 1S. At projected L different from tuple L, yield columns scale and cross sections remain unchanged. Constant factor discrepancies suggest B, luminosity, units or sample-fraction mistakes. This does not authenticate a CSV's source file hash.

The conversion is `sigma_gamma_p[pb] = corrected_visible / (1000 × L[fb^-1] × B_ee × flux)`, using the selected state's B_ee. Use `--br-ee 0.0234` (or `--branching-fraction 0.0234`) only when matching an explicitly chosen historical 1S branching convention. The supplied cross section is **visible dielectron**, so B is divided out once at this stage. The scope is one generated state, full-sample, raw-state-normalization bin-by-bin extraction. Do not pass train/test unfolded CSVs or multi-state tuples and interpret them as this calculation. Historical truth-state rescaling requires separate matching settings; the all-generated uniform counterfactual here is explicitly different.

The optional **1S-only** real-photon benchmark uses `6.4 W^0.74 [1−((m_Upsilon+m_p)/W)^2]^2 pb` above threshold, averaged with each flux. Source: https://eic-code-browser.sdcc.bnl.gov/lxr/source/estarlight/src/photonNucleusCrossSection.cpp . This is a current-source Q2=0 diagnostic, not the campaign generator revision or a full finite-Q2 prediction. For 2S/3S, or an overridden 1S mass, it is marked `NOT_APPLICABLE`; model values and ratios are null. Flux integration and normalization checks still run. No 1S benchmark is relabelled as an excited-state prediction. Do not choose a normalization by forcing agreement with it.

To compare an independent production rate, add `--generator-visible-xsec-fb VALUE --generator-source 'exact card/revision/log reference'`. Both rates must cover the same generated phase space and include the same branching convention. The provenance string is recorded, not authenticated. The reference logs, cards and revision establish the quoted July baseline; the original production log/revision is still needed to authenticate the existing EOS sample provenance.

## Outputs and return codes

Each output directory contains `REPORT.md`, `audit.json`, `weight_summary.csv`, `weight_groups.csv`, `w_bin_audit.csv`, and `cutflow.csv`. Pair, saved-CSV and cross-energy comparison CSVs are written when those checks run. Numerical evidence and optional beam comparisons are complete in `audit.json`.

Start with `REPORT.md`, then `w_bin_audit.csv`. For review, send the two `audit.json` files, both `w_bin_audit.csv` files and the 9x275 `energy_comparison.csv`.

Exit 0 means no internal FAIL/ERROR, but WARN/NOT_CHECKED can remain; exit 1 means a failed consistency check; exit 2 means an execution/input error. External rate provenance is deliberately not certified. The script does not test unfolding response matrices or revalidate Garwood/statistical-uncertainty formulas.

SHA256 and ordered full-branch comparison read the files several times. On EOS this can take minutes. There is no sampling of rows in the pair comparison; output metadata previews and weight-group retention are bounded and marked when truncated.

## Verification

`tests/test_normalization_audit.py` builds controlled ROOT fixtures covering unequal job sizes, shape changes, valid uniform preparation, missing/changed branches, wrong denominators, weight/projection mutations, source hashes, malformed luminosity metadata, projected luminosity, CSV factors and two-energy decomposition. `tests/test_all_states.py` checks state mass/branching/Q² defaults, state CSV discovery and filtering, wrong-state input rejection, explicit Q² overrides and the absence of excited-state model predictions. Run with `python3 -m unittest discover -s tools/normalization_audit/tests -v`.

`verification/` preserves earlier test logs and read-only checks on three archived ROOT files using their historical settings. `verification/jul26_tests.log` records the current regression tests. These are not executions on current EOS inputs. No new normalized ROOT was generated as part of those archived-input checks.
