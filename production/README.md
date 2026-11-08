# Production sources

To produce new tuples, start with the [one-file 2S 9×130 pilot](../docs/production.md#one-file-pilot). To analyse the distributed merged files, start with the [paper replay](../paper/README.md).

| Component | Role |
| --- | --- |
| [9×130 macro](macros/9x130/DVMP_Upsilon_Analysis.C) and [9×275 macro](macros/DVMP_Upsilon_Analysis.C) | Beam-specific producer snapshots. Both use the same reconstructed [electron finder](macros/ElectronFinder.C). |
| [Original 2S 9×130 campaign](cern-2026-10/README.md) | Original workers and a verifier that checks the private job records against the merged sample |
| [Production planner](jobs.py) | Creates new Condor batches with frozen macros and checks event counts when merging |
| [Simulation cards](cards/) and [downloader](download/download_upsilon_fast.py) | Six official July cards and the reconstructed-input download script |
| [Source inventory](sources.json) | Origin and SHA-256 of each preserved source |

Two 100-event runs, one at each beam energy, reproduced every saved tuple branch. To compare a fresh pilot with the distributed events, follow the [one-file commands](../docs/production.md#one-file-pilot). [Results and scope by sample](../docs/missing-production-sources.md#two-reconstructed-inputs)
