# Rerun the Υ analysis with corrected ECAL matches

**Starting with Υ(1S) at 9×275 on lxplus?** Follow the [Condor walkthrough](condor-rematch-1s-9x275.md) for the one-file test, worker logs, full submission and merge.

Start from the July campaign's original `*.eicrecon.edm4eic.root` files. They contain the track projections and ECAL clusters needed to rebuild the associations affected by [EICrecon #3011](https://github.com/eic/EICrecon/issues/3011). The merged `UpsilonTree` files contain the resulting electron candidates and serve as the reference. Generation and detector simulation can be reused.

The producer option `--electron-id reco-rematch3011` replays the original matching with its stored cluster order and verifies the ordered track and cluster IDs against the saved links in each ECAL, event by event. A disagreement stops that job. It then uses the physical track-segment index on every iteration to build the corrected links. The reconstructed electron ID applies the same unique-match, E/p and isolation cuts to those links. The tuple records reconstructed electron ID as `ele_id_mode=1`; the job and merge records distinguish stored from rematched links.

## Start on lxplus and run one file

From the repository root on **bare lxplus**, read the original Υ(2S), 9×130 file under `/eos/cms`. If you are already at the ePIC shell prompt, type `exit` first. Enter the mounted path to your campaign's `upsilon_data` directory at the prompt. The expected four bytes are `root`.

```bash
read -r -p "Absolute /eos/cms path to upsilon_data: " RECO_INPUT_BASE
name=eSTARlight1.3.0-1.0_Upsilon2S_9x130_q2_0to0.01_hiAcc_ab_run015.0927.eicrecon.edm4eic.root
PILOT_2S="${RECO_INPUT_BASE%/}/Upsilon2S/9x130/hiAcc/$name"
head -c 4 -- "$PILOT_2S" && echo
```

If the host reports `Permission denied`, check the host ticket with `klist -s`. If it has expired, run `kinit "$USER@CERN.CH"` on lxplus and repeat the read. If a valid ticket still gets denied, check the account with `eos root://eoscms.cern.ch whoami` and the read/traverse permissions for this group path with its owner.

Once the host prints `root`, launch the July ePIC image from the **same lxplus shell**. This block mounts the CMS EOS filesystem and passes the host's Kerberos FILE cache at exactly the path the session uses. The ticket is mounted read-only for this interactive session. [Reported lxplus fix](https://github.com/apptainer/apptainer/issues/1217#issuecomment-3690453554)

```bash
(
  set -e
  cache=$(LC_ALL=C klist | sed -n 's/^Ticket cache: //p')
  case "$cache" in
    FILE:/*) ticket=${cache#FILE:} ;;
    DIR:*|KCM:*|KEYRING:*)
      private_dir=$(mktemp -d /tmp/eos-pilot-krb.XXXXXXXX)
      ticket="$private_dir/krb5cc"
      cache="FILE:$ticket"
      trap 'kdestroy -c "$cache" 2>/dev/null || :; rmdir -- "$private_dir" 2>/dev/null || :' EXIT
      kinit -c "$cache" "$USER@CERN.CH"
      ;;
    *) echo 'Could not identify a supported ticket cache. Run klist on lxplus.' >&2; exit 1 ;;
  esac
  if ! [[ -f "$ticket" && -r "$ticket" && -O "$ticket" ]]; then
    echo 'The host ticket file is not readable or does not belong to this user.' >&2
    exit 1
  fi
  image=/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:26.07.1-stable
  unset -f which 2>/dev/null || true
  apptainer exec --bind /afs,/eos/cms --bind "$ticket:$ticket:ro" \
    --env "KRB5CCNAME=$cache" --env "PILOT_2S=$PILOT_2S" --pwd "$PWD" "$image" eic-shell
)
```

For `KCM`, `DIR` or `KEYRING`, the block requests a separate ticket on lxplus, keeps it in a private temporary directory, and destroys it when you exit the container. The [EIC shell](https://github.com/eic/EICrecon/blob/main/docs/tutorial/02-work-environment.md) starts in your checkout. The [July release](https://github.com/eic/containers/blob/2def1e04f273ffa1c87f388f5c65e8d4b89180ec/spack-environment/packages.yaml) specifies ROOT 6.38.04, podio 1.7 and C++20.

**Inside the ePIC shell**, read the same file again. Continue to the planner only when this also prints `root`. If only this read fails, check that `KRB5CCNAME` points to the mounted file and that `klist -s` succeeds inside. If you renew the ticket on the host, re-enter the container to pick up the current file.

If startup prints `Illegal option --`, check the container name at the current prompt before launching another shell. An exported lxplus `which` function can cause this message. The `unset -f which` above removes that function for the launch only. [EIC example](https://eic-code-browser.sdcc.bnl.gov/lxr/source/epic-lfhcal-tbana/tools/run-in-eic-shell.sh)

```bash
printf 'container: %s\n' "${APPTAINER_CONTAINER:-${SINGULARITY_CONTAINER:-host}}"
head -c 4 -- "$PILOT_2S" && echo
```

Run the tuple worker on that EOS path. Only the tuple output and job records go into the checkout's `workspaces` directory.

```bash
(
  set -e
  root-config --version
  python3 -c 'import uproot'
  patch --version | head -1
  mkdir -p workspaces
  trial=$(mktemp -d "$PWD/workspaces/rematch-2S-9x130.XXXXXX")
  input="${PILOT_2S:?Set the original EOS pilot path when entering the container}"
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

For the full 1S 9×275 sample, the [Condor launcher](condor-rematch-1s-9x275.md#plan-and-submit-the-full-sample) writes each per-input tuple directly to EOS and groups up to 100 inputs per Condor job. For the other samples, add `--files-per-job 100` to the planner command below if appropriate for the worker runtime.

New rematching plans request `+AccountingGroup = "group_u_CMS.CAF.COMM"` in both `submit.sub` and `pending.sub`. To check before submission, run `grep '^+AccountingGroup' "$batch/pending.sub"` after creating the pending file. A `worker --execute` run in your shell uses no Condor account.

Use one full input list for each of the three states at each beam energy. Check the file lists against the campaign inventory. Set `CORRECTED_BASE` to a **fresh absolute directory on shared storage**. For each list, set `state`, `beam` and `list`. For example, for the first production batch:

```bash
state=1S
beam=9x275
list="$PWD/workspaces/manifests/9x275-1S-inputs.txt"
batch="$PWD/workspaces/production/rematch3011-$state-$beam"
: "${CORRECTED_BASE:?Set the absolute corrected output base on EOS}"
if [ "$state" = 1S ]; then subdir=output_newID; else subdir=output; fi
out="$CORRECTED_BASE/Upsilon$state/$beam/hiAcc/$subdir"
python3 production/jobs.py plan --input-list "$list" --work-dir "$batch" \
  --tuple-output-dir "$out/intermediate" \
  --beam "$beam" --state "$state" --electron-id reco-rematch3011 --femc full \
  --setup-script "$PRODUCTION_SETUP" --environment-label "$PRODUCTION_ENV_LABEL" \
  --hash-inputs --write
worker="$batch/snapshot/jobs.py"
python3 "$worker" worker "$batch/jobs/job_000000/job.json" --execute
python3 "$worker" pending --batch "$batch" --write
# On bare lxplus, submit the untouched jobs once using the July image.
image=/cvmfs/singularity.opensciencegrid.org/eicweb/eic_xl:26.07.1-stable
condor_submit -append 'getenv = False' \
  -append "+SingularityImage = \"$image\"" "$batch/pending.sub"
```

When all jobs in that batch finish, merge them into the layout used by the analysis. The merge validates the generated counts. Compare the final row total with the original state/energy tuple before interpreting a yield change.

```bash
python3 production/jobs.py merge --batch "$batch" --output "$out/merged_Upsilon_final.root" --execute
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
