# Bookkeeping and existing work

The repository contains source code and declared configurations. Runtime
workspaces contain a fixed copy of that source, generated normalized inputs,
run manifests and outputs. Initialize a new workspace when the source changes;
keep old workspaces and run IDs intact.

Existing CERN directories remain independent archives:

- `/afs/cern.ch/work/g/gkrintir/private/Ypsilon-Fit-Framework`
- `/afs/cern.ch/work/g/gkrintir/private/Ypsilon-Fit-Framework-ProtonTagged`
- the corresponding EOS/SWAN result mirrors

Do not install this repository over them or copy a new producer into their run
folders. Plotting profiles can explicitly read existing runs without changing
those runs. The source receipts in old runs remain valid for their recorded code.

## Resume versus a new calculation

Use resume only when the runner supports it and input/configuration/code hashes
match. Changes to selections, luminosity, source files or response/binning
settings require a fresh run. `--keep-going` permits independent analysis tasks
to continue while affected dependent plots remain blocked. It does not convert a
failed task to success.

## Reference outcomes already reviewed

The initial repository is not a new physics validation. Relevant available
records include:

| Run | Meaning |
|---|---|
| `allstates-proton-projected-L2p5-v2-{1s,2s,3s}` | Existing 9×275 full-analysis source runs |
| `reopt-proton-9x275-L2p5-v1` | Existing tagged optimization, with per-observable validation outcomes |
| `allstates-proton-9x130-projected-L2p5-20261101T164736Z-1s` | Existing 9×130 source; its actual projection command used 1 fb⁻¹ |
| `reopt-proton-9x130-L1p0-20261101T172820Z` | 1S, 1 fb⁻¹: both slope observables have no feasible binning |

The last run has 8.562 expected selected truth-fiducial events at 1 fb⁻¹. Its pT²
candidates fail expected-count requirements; true-|t| via pT² candidates also
fail migration purity. There is no new slope or differential uncertainty to quote
from that run. The original unsupported true-|t|1–5 tail has zero surviving
training events after tagging and two held-out survivors. This is not a reason
to bypass the response-support guard.

Keep ROOT files on EOS/AFS and reference them through configuration. Commit
source changes, input definitions and small reference tables. Generated plots,
large production ROOT files, local paths, credentials and job logs are excluded
from Git. Store an exported run manifest with any figures used in the note.
