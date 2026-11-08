# CERN tuple production for Υ(2S), 9×130 GeV

On October 21–22, 2026, 76 HTCondor jobs processed 10,583 reconstructed files. Each supplied 100 events. For every input the logs record a completed tuple, reconstructed electron ID and the full fEMC. All jobs exited with status zero. The original driver printed a success message even after an individual ROOT failure, so the audit matched each input with its own tuple-write message.

The [original Condor wrapper](condor_wrapper.sh) entered the transferred analysis directory and invoked the [per-file worker](run_Upsilon_Analysis_new_condor.sh). The worker passed only the input and output paths to the macro. The jobs selected the floating `eic_xl:nightly` image. The [campaign record](campaign.json) lists the recorded settings and aggregate checks. For a new run use the [production instructions](../../docs/production.md).

## Check the archived merge locally

Run the following from the repository root, using the original archive and the distributed 2S 9×130 merged ROOT file.

```bash
python3 production/cern-2026-10/check_merge.py /path/to/upsilon-tuple-sources.zip /path/to/merged_Upsilon_final.root
```

The checker reads the submission lists and logs directly from the local ZIP. It verifies the reference ROOT checksum and sums four electron-ID counters in each successive block of 100 events. All 10,583 blocks match the logged totals in lexicographic output filename order. Only nine match submission order. There are 3,065 distinct combinations of the four counters, and 1,201 occur just once. This links the campaign to the distributed tuple while leaving the original merge command and intermediate file checksums unresolved.

To get the first input path for a pilot on CERN infrastructure, read it from your local archive:

```bash
unzip -p /path/to/upsilon-tuple-sources.zip eic/list.txt | head -n 1
```

The other recovered submission templates record possible campaigns, but only the 2S 9×130 sample has the completed worker records checked here.
