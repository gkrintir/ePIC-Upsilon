# Rerun the Υ analysis with corrected ECAL matches

Start from the July campaign's original `*.eicrecon.edm4eic.root` files. They contain the track projections and ECAL clusters needed to rebuild the associations affected by [EICrecon #3011](https://github.com/eic/EICrecon/issues/3011). The merged `UpsilonTree` files contain the resulting electron candidates and serve as the reference. Generation and detector simulation can be reused.

The producer option `--electron-id reco-rematch3011` replays the original matching with its stored cluster order and verifies the ordered track and cluster IDs against the saved links in each ECAL, event by event. A disagreement stops that job. It then uses the physical track-segment index on every iteration to build the corrected links. The reconstructed electron ID applies the same unique-match, E/p and isolation cuts to those links. The tuple records reconstructed electron ID as `ele_id_mode=1`; the job and merge records distinguish stored from rematched links.

## Start on lxplus and run one file

Start in the root of the checkout to which you applied the rematching patch. On the lxplus host, check the patched option and enter the July ePIC software release. Keep your checkout on AFS or on the bound EOS mount.

```bash
python3 production/jobs.py plan --help | grep reco-rematch3011
image=/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:26.07.1-stable
apptainer exec --bind /afs,/eos/cms --pwd "$PWD" "$image" eic-shell
```

The last line opens a new shell in the same directory with ROOT and podio loaded. This uses the [EIC shell](https://github.com/eic/EICrecon/blob/main/docs/tutorial/02-work-environment.md) and the [published July image](https://hub.docker.com/r/eicweb/eic_xl/tags?name=26.07.1-stable). Its [release manifest](https://github.com/eic/containers/blob/2def1e04f273ffa1c87f388f5c65e8d4b89180ec/spack-environment/packages.yaml) specifies ROOT 6.38.04, podio 1.7 and C++20. The archived October jobs used the floating `nightly` tag. The interactive example uses the named release.

Inside that shell, run one original Υ(2S), 9×130 file. The following block downloads the public copy into a fresh working directory and checks its identity. If you have the original mounted on EOS, you can instead set `input` to its absolute filename and skip the download.

```bash
(
  set -e
  root-config --version
  python3 -c 'import uproot'
  patch --version | head -1
  mkdir -p workspaces
  trial=$(mktemp -d "$PWD/workspaces/rematch-2S-9x130.XXXXXX")
  name=eSTARlight1.3.0-1.0_Upsilon2S_9x130_q2_0to0.01_hiAcc_ab_run015.0927.eicrecon.edm4eic.root
  input="$trial/$name"
  curl -fL --retry 2 "https://cernbox.cern.ch/remote.php/dav/public-files/XVLz7BDkpAHgy17/$name" -o "$input"
  printf '%s  %s\n' b5f8450deba55cc4555532f11f87be1d81713ef3cefc51774f84aa7b756a01ff "$input" | sha256sum -c -
  printf '%s\n' "$input" > "$trial/input.txt"
  batch="$trial/corrected"
  python3 production/jobs.py plan --input-list "$trial/input.txt" --work-dir "$batch" \
    --state 2S --beam 9x130 --electron-id reco-rematch3011 --femc full \
    --environment-label eic_xl:26.07.1-stable --hash-inputs --write
  python3 "$batch/snapshot/jobs.py" worker "$batch/jobs/job_000000/job.json" --execute
  echo "Corrected tuple: $batch/jobs/job_000000/output.root"
)
```

The expected worker result is `DONE` with 100 output rows. The corresponding comparison below gives 5 proton-tagged events with the stored links and 9 with the correction. The one-file command runs on the login node. When submitting a full Condor batch, configure the worker to use the same tested environment. The generated submit file does not select this image.

## Compare the two pilots

For the two-file stored-versus-corrected check, set `PILOT_1S` and `PILOT_2S` to the absolute paths of the two [original inputs](https://cernbox.cern.ch/s/XVLz7BDkpAHgy17). Use release label `eic_xl:26.07.1-stable` if you opened the shell above. In an initialized shell you can run the planner directly. If `PRODUCTION_SETUP` is set, the following example sources that script for ROOT subprocesses.

```bash
set -e
: "${PILOT_1S:?Set path to 1S pilot}" "${PILOT_2S:?Set path to 2S pilot}"
environment_args=(--environment-label "${PRODUCTION_ENV_LABEL:-eic_xl:26.07.1-stable}")
if [ -n "${PRODUCTION_SETUP:-}" ]; then environment_args+=(--setup-script "$PRODUCTION_SETUP"); fi
pilot_base="$PWD/workspaces/rematch3011-pilots"
mkdir -p "$pilot_base"
for state in 1S 2S; do
  if [ "$state" = 1S ]; then beam=9x275; input="$PILOT_1S"
  else beam=9x130; input="$PILOT_2S"; fi
  list="$pilot_base/$state-$beam.txt"
  printf '%s\n' "$input" > "$list"
  for mode in reco reco-rematch3011; do
    batch="$pilot_base/$state-$beam-$mode"
    python3 production/jobs.py plan --input-list "$list" --work-dir "$batch" \
      --beam "$beam" --state "$state" --electron-id "$mode" --femc full \
      "${environment_args[@]}" \
      --hash-inputs --write
    python3 "$batch/snapshot/jobs.py" worker "$batch/jobs/job_000000/job.json" --execute
  done
  python3 production/check_rematch_pilot.py \
    "$pilot_base/$state-$beam-reco/jobs/job_000000/job.json" \
    "$pilot_base/$state-$beam-reco-rematch3011/jobs/job_000000/job.json"
done
```

The check compares the selected **event indices** to an independent replay, at the event cuts and after tagging the recoil proton. It also compares the generated quantities and ordinary event weight for every row. With 100 generated events per file, the expected counts are:

| Pilot | Before proton tag, stored → corrected | Proton tagged, stored → corrected |
|---|---:|---:|
| Υ(1S), 9×275 | 13 → 23 | 8 → 13 |
| Υ(2S), 9×130 | 23 → 33 | 5 → 9 |

These counts use 8–12 GeV in dielectron mass, 20–100 GeV in reconstructed W and the longitudinal-balance cut. They apply only to the two identified files.

## Produce the six corrected tuples

Use one full input list for each of the three states at each beam energy. Check the file lists against the campaign inventory. Set `CORRECTED_BASE` to a **fresh absolute directory on shared storage**. For each list, set `state`, `beam` and `list`. For example, for the first production batch:

```bash
state=1S
beam=9x275
list="$PWD/workspaces/manifests/9x275-1S-inputs.txt"
batch="$PWD/workspaces/production/rematch3011-$state-$beam"
python3 production/jobs.py plan --input-list "$list" --work-dir "$batch" \
  --beam "$beam" --state "$state" --electron-id reco-rematch3011 --femc full \
  --setup-script "$PRODUCTION_SETUP" --environment-label "$PRODUCTION_ENV_LABEL" \
  --hash-inputs --write
worker="$batch/snapshot/jobs.py"
python3 "$worker" worker "$batch/jobs/job_000000/job.json" --execute
python3 "$worker" pending --batch "$batch" --write
# Submit once for a multi-file batch after checking the pilot job.
condor_submit "$batch/pending.sub"
```

When all jobs in that batch finish, merge them into the layout used by the analysis. The merge validates the generated counts. Compare the final row total with the original state/energy tuple before interpreting a yield change.

```bash
: "${CORRECTED_BASE:?Set a new shared output directory}"
if [ "$state" = 1S ]; then subdir=output_newID; else subdir=output; fi
out="$CORRECTED_BASE/Upsilon$state/$beam/hiAcc/$subdir"
mkdir -p "$out"
python3 "$worker" merge --batch "$batch" --output "$out/merged_Upsilon_final.root" --execute
```

## Recalculate the physics results

Use a new analysis configuration and workspace. Once all three corrected 9×275 tuples exist, the tagged paper setting runs as follows:

```bash
export CORRECTED_BASE
config="$PWD/workspaces/rematch3011-analysis.json"
python3 - "$config" <<'PY'
import json, os, sys
from pathlib import Path
base = Path(os.environ["CORRECTED_BASE"]).resolve(strict=True)
with Path(sys.argv[1]).open("x") as out:
    json.dump({"workspace_root": "workspaces", "input_base": str(base)}, out)
PY
run() { python3 upsilon.py --config "$config" "$@"; }
ws=tagged-rematch3011
run init --name "$ws" --selection proton-tagged
run doctor --workspace "$ws"
run analysis run --workspace "$ws" --energy 9x275 --states 1S 2S 3S \
  --lumi 2.5 --batch-id rematch275-L2p5 --hash-inputs
run reopt run --workspace "$ws" --energy 9x275 --states 1S 2S 3S \
  --observables pt2 abs_t_via_pt2 w y --lumi 2.5 \
  --source-batch rematch275-L2p5 --run-id rematch275-fit-L2p5 --hash-inputs
run status --workspace "$ws" --run-id rematch275-fit-L2p5
```

Use separate batch and run IDs for each luminosity or beam configuration, and a separate `inclusive` workspace for the proton-optional comparison. Start the old/new comparison with the same selections, paper bin boundaries and generated normalization. Recompute efficiency, migration and feed-in from the corrected tuples. The measured change in the two pilots is not an efficiency correction for the full samples. The historical `paper/replay_paper.py` checks the original tuple hashes and published numbers. Use it on the original inputs as a baseline.

The analysis command computes the standard observables. For the Appendix D trace-anomaly figure, run its separate differential-|t| fit on the corrected Υ(1S) tuples at each beam energy with 2.5 and 100 fb⁻¹.
