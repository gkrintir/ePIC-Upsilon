# Verification of the external 9×130 update

Software verification completed on 2026-11-01 using the supplied proton-tagged
framework snapshot. **18 tests passed.** They cover:

- Source energy, proton contract and prepared-input provenance guards; partial
  source runs accepted when preparation succeeded.
- Synthetic normalized ROOT input validation at source 10 fb⁻¹ and target 1 fb⁻¹.
- Entire external runner execution through the unavailable-case report, including
  9×130 captions, overall status propagation, input preservation, output collision
  refusal, identical-setting resume and changed-luminosity resume rejection.
- Exact legacy support reconstruction against the framework's existing response
  implementation; zero training support versus proton-cut losses, finite-value
  prefilter effects, negative reconstructed pT and outside-range feed-in.
- Required category enforcement, invalid tag schema and fresh-output safety.
- Independent numerical EPA integration in y and log(Q²), explicit state rates
  and branching fractions, additive W-bin flux and rapidity normalization.
- Rejection of an unsupported fitted truth bin, retained outside-range truth
  feed-in and reconstructed control interval, and distinction between truth pT²
  and truth |t|.

The shell quick-start passes `bash -n`. All **1,605 original files** present in
the supplied framework ZIP were compared byte-for-byte against the local
framework snapshot after testing: zero missing, zero changed. Test fixtures and
output reside in temporary directories. No framework source files, requirements,
existing results or ROOT inputs are changed by this add-on.

Numerical flux checks for Ee=9 GeV, Ep=130 GeV, W20–70 GeV:

| Q²max [GeV²] | Framework Φ | Independent integration Φ |
|---:|---:|---:|
| 1 | 0.06880106724446178 | 0.06880106724445878 |
| 0.01 | 0.04959731276197783 | 0.04959731276197786 |

These are numerical implementation checks, not validation of the production
executable, acceptance model or data sensitivity. The CERN ROOT files were not
available here. No optimized edges or physics uncertainty claim is supplied;
those must come from the new run and its validation status.

Optional local test command, from the directory containing this package and the
framework, using the analysis Python environment:

```bash
YPSILON_TEST_FRAMEWORK="$PWD/Ypsilon-Fit-Framework-ProtonTagged" \
  python3 -B -m unittest discover \
  -s Ypsilon_9x130_Reoptimization_Update/tests -v
```

No need to run these tests before the quick-start; its `validate` step checks
the actual normalized ROOT input before the audit and optimization.
