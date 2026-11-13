# Source revisions and remote hosting

The public source repository uses `main` for development. Clone it to start a new analysis.

```bash
git clone https://github.com/gkrintir/ePIC-Upsilon.git
cd ePIC-Upsilon
git status
git log -1 --oneline
```

The clone configures `origin` automatically. Sign in to GitHub when pushing changes.

Code, reference tables, source provenance and documentation are tracked. ROOT
data, runtime workspaces, generated results, local paths and environments are
ignored. The few historical replay logs under `analysis/reference_replay/` are
fixed reference evidence, not new runtime logs.

## Using the earlier ZIP

The previously supplied `ePIC-Upsilon-repository.zip` contains a standalone
snapshot and its local, prepublication Git history. It remains usable offline
when unpacked into its own new directory. It is not a checkout of the hosted
repository's history.

For GitHub work, clone the hosted repository separately and transfer any intended
local source edits into a new branch for review. Preserve the ZIP directory and
its runtime workspaces as bookkeeping. Do not force-push the archive's history
over `origin/main` or overwrite an existing framework with either source copy.

## Make a reviewed source change

```bash
git switch -c describe-your-change
# Edit code or documentation.
python3 tools/run_checks.py
git diff
git status
```

Commit only the intended files. The repository checker verifies the recovered
analysis sources against `provenance/analysis_snapshot.json`. Keep the original
`copied_files` hashes when making a reviewed source update. Add an `amended_files`
record with `path`, `baseline_sha256`, `new_sha256` and a nonempty `reason`, and
perform the relevant numerical validation. The checker compares the recorded
baseline hash with `copied_files` and the current file with `new_sha256`.

Initialize a fresh analysis workspace after a code update. Existing workspaces
retain their source snapshots, selection contracts and results. Plotting can
still read explicit old runs; it does not relabel them as new results.

Record a source revision alongside an analysis decision with `git rev-parse HEAD`.
The runtime manifests additionally record the effective configuration and input
identities. A Git commit alone does not identify the ROOT sample or certify the
physics validation.
