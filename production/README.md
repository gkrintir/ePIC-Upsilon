# Production sources

Start with [the production guide](../docs/production.md). Existing merged ROOT files can go directly to the analysis workflow; this directory is needed when regenerating tuples.

| Content | Status |
| --- | --- |
| `macros/DVMP_Upsilon_Analysis.C`, `ElectronFinder.C` | Latest recovered uploads, preserved byte for byte |
| `macros/ePICStyle.C`, `ePICStyle.h` | Recovered from the user's GitHub repository, pinned by commit |
| `cards/` | Six official Jul26 cards, pinned by commit |
| `download/download_upsilon_fast.py` | Recovered October 21 downloader, preserved byte for byte |
| `jobs.py` | **New** explicit production planner, Condor worker and checked merge helper |
| `sources.json` | Origins, SHA256 hashes and unresolved provenance |

The original CERN Condor submit, worker/setup and merge scripts have **not** been recovered. `jobs.py` is a new implementation, not a reconstruction of those files. The preserved producer hardcodes **9×275**, uses historical per-job raw weights and is not authenticated as the exact executable used for the EOS samples. New 9×130 production is blocked until the actual per-energy producer is supplied and audited. Existing 9×130 analysis is unaffected.

No ROOT file, Condor job or remote download is created merely by reading this directory. `jobs.py plan` is read only; `--write` creates fresh job descriptions, and worker/merge require `--execute`. Submission remains an explicit `condor_submit` command.
