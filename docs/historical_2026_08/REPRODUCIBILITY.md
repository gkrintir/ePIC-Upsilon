# Reproducibility workflow

## 1. Pin input identities

Copy `config/local.example.toml` to `config/local.toml` and set the exact ROOT
files used for low and high beam energy. Do not infer sample identity from a
similar filename. Record SHA-256 values in `[checksums]`; the runner refuses a
configured mismatch.

The six validated report inputs are pinned in
`config/verified_root_inputs.toml`. Keep legacy and updated tuples in their
separate `low_legacy_root`/`high_legacy_root` and `low_root`/`high_root` slots.

The schema preflight profiles tree names and branches before extraction:

```bash
python workflow/validate_root_schema.py input.root \
  --profile differential \
  --out branch_inventory.json
```

## 2. Inspect the resolved plan

Global options must precede the subcommand:

```bash
python workflow/run_analysis.py --config config/local.toml plan report
python workflow/run_analysis.py --config config/local.toml validate report
python workflow/run_analysis.py --config config/local.toml plan report_crosscheck
```

`plan` resolves the dependency graph and prints every command without importing
ROOT. `validate` compiles promoted scripts, checks declared modules, executables,
input files, output locations, and configured input hashes.

## 3. Execute in isolated stage directories

```bash
python workflow/run_analysis.py --config config/local.toml \
  run report --keep-going --hash-inputs
```

Each task runs in its own directory below
`results/generated/runs/<run-id>/stages/`. This isolates ROOT global state and
prevents notebook-style working-directory leakage. A failed independent task
does not block unrelated branches when `--keep-going` is used.

The runner uses these statuses:

| Status | Meaning |
|---|---|
| `READY` | Prerequisites passed; dry-run command is available |
| `CACHED` | Every declared output exists and is nonempty |
| `DONE` | Process completed and declared outputs were produced |
| `BLOCKED_INPUT` | Code, file, software, value, or output-location check failed |
| `BLOCKED_DEPENDENCY` | A prerequisite stage did not produce a usable result |
| `SKIPPED_EXPERIMENTAL` | Opt-in task was intentionally omitted |
| `INVALID_OUTPUT` | Output missing or log matched a configured failure pattern |
| `FAILED` | Subprocess returned a nonzero exit code |

Mass-fit stages treat `NOT POSDEF`, `CALL LIMIT`, and `Abnormal termination` as
invalid output even if ROOT returns process status zero.

## 4. Review provenance manifests

Each executed stage records:

- resolved command and working directory;
- script SHA-256;
- configured input paths, file sizes, and timestamps;
- optional input SHA-256 when `--hash-inputs` is selected;
- output paths, sizes, timestamps, and SHA-256 values;
- start/end time, duration, log path, status, and validation reasons.

The run manifest embeds the complete resolved configuration and all stage
records. Keep it with any promoted figures or tables.

## 5. Gate report promotion

Before copying a generated artifact into a report:

1. Verify the originating stage is `DONE` or a checksum-identical `CACHED` run.
2. Confirm every requested state had positive generator weight and selected
   reconstructed events.
3. Confirm efficiencies, purity, and stability are finite and within [0, 1].
4. Confirm response-row sums match efficiencies for response unfolding.
5. Confirm all plotted points are finite and the filename phase space matches
   metadata.
6. Confirm fit minimization and covariance are valid for fitted results.
7. Confirm normalization scope, visible/inclusive convention, W range, Q2 cut,
   luminosity, random seed, iterations, and toy count in the run manifest.
8. Record the report figure mapping in `docs/report_output_map.csv` or a future
   report lock file.

## Environment notes

The supplied archive had no dependency lock. The framework declares the Python
packages inferred from imports, but PyROOT must normally come from a compatible
CERN ROOT, LCG, or SWAN environment. eSTARlight is a separate external program.
The current package can still run static, configuration, and CSV-quality tests
without either dependency. The macro-independent report cross-check requires
uproot and NumPy but not PyROOT. `requirements-reference.txt` pins the versions
used for byte-exact archived PNG reproduction.
