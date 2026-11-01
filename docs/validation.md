# Repository validation

Run `python3 tools/run_checks.py` from the repository root after installing
`requirements-dev.txt`. Each component runs in its own process and working
directory so historical modules with the same names cannot shadow one another.

## Checks performed during assembly

The local checks used Python 3.12. Python source was also parsed against the
Python 3.9 grammar; this syntax check is not a full Python 3.9 runtime test.

| Suite | Result | Scope |
|---|---:|---|
| Repository interface | 32 passed | Workspace isolation and source guards, actual backend planning, dependency imports, plot profiles, production planning/merge bookkeeping and source collection |
| Preserved analysis backend | 32 passed | Recoil-proton selection, unchanged truth normalization, orchestration, category contracts and reoptimization selection |
| Plotting | 81 passed | Saved numerical inputs/statuses, covariance-aware display merges, unavailable cases, mass fits/intervals and reference conventions |
| Efficiency extension | 29 passed | Final selections, topology consistency, truth denominators, chunking, coordinates, maps and count intervals |
| 9×130 extension | 18 passed | Configuration/source guards, synthetic ROOT bookkeeping, support audits and unavailable inference |

The repository checker confirms that all **122 files in the current imported
analysis baseline are byte-identical** to their supplied source snapshot. It
also checks Python syntax, clean notebook outputs, main documentation links and
absence of detector ROOT files or embedded archives in the source tree.

The mass/reconstruction, W and differential figures were rendered from the
supplied 9×275 proton-tagged numerical exports. PDF/PNG and CSV/JSON companions
were produced. W/rapidity rendering used explicit diagnostic opt-in because the
saved cases are NOT_VALIDATED. These statuses and unavailable higher-state pT²
cases remained in the audits; rendering did not promote them. The differential
figure and the README SVG were visually inspected.

The W notebook's display-rebin switch had a later hard-coded override; this was
removed. The active rebinning module was older than its supplied notebook/tests;
the archived update for saved-edge adjacent merges and proton-contract metadata
was recovered. Origins and before/after hashes are in
`plotting/source_inventory.json` and `plotting/SOURCE_NOTES.md`. No inference
formula was changed for repository assembly.

The standalone four-map configuration recipe and its exact plan command were
tested with explicit temporary paths. It did not read CERN event data.

## Not established by these checks

- No full physics workflow was rerun against the large EOS ROOT samples in this
  environment. Stored figure exports are not a substitute for that calculation.
- ROOT/podio compilation, a representative real producer job, CERN authentication,
  shared-filesystem behavior and a Condor pilot still require the actual CERN
  environment. Synthetic worker/merge tests check bookkeeping only.
- The original CERN batch scripts, exact producer revisions by energy and tested
  environment setup remain to be reconciled. The recovered 275 GeV producer is
  not accepted for new 9×130 production.
- Combined notebooks 06/07 were not executed end to end against the original
  ROOT files. Notebook 07's W–xV cache/consistency calculation requires that input.
- Existing conditional or failed statistical validation, signal-only assumptions,
  response limitations and omitted experimental uncertainties are unchanged.

Keep the original analysis runs and their numerical audits. Use fresh workspace
and run names for any subsequent source or physics change.
