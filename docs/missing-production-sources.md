# Production sources still to reconcile

The current analysis and plotting sources are recovered. The exact CERN producer
and batch environment used for the latest EOS production are not fully recovered.
This distinction matters because a source file with a familiar name is not proof
that it generated a given dataset.

| Component | Available | Remaining check |
|---|---|---|
| DVMP_Upsilon_Analysis.C | Uploaded October 18 version, preserved byte-for-byte | Hardcoded proton energy 275 and old raw-weight constants; need actual 9×130/9×275 CERN copies |
| ElectronFinder.C | Uploaded October 17 updated reconstructed-association version | Compare checksum with the actual production copy |
| ePICStyle.C/.h | Recovered from gkrintir/lAger_JPsi_MuMu_Dataset at recorded commit | Compare with CERN copy |
| Jul26 eSTARlight cards | Six cards at commit f62965a1612fc734dbab18e141f4809740d6b234 | Verify production executable/container revision and exact cards used |
| Condor submit/worker scripts | New documented planning wrapper in this repository | Original CERN submission and environment scripts still needed |
| ROOT merging | New explicit-input helper where provided | Original file list, duplicate handling and job-completeness checks still needed |
| ROOT/podio/eic-shell environment | Dependencies identified | Actual setup script or container tag/digest needed |

Do not infer that 9×130 EOS data used a 275 GeV producer from the recovered macro.
Separate local source copies may have been used. Do not change analysis rates,
selection branches or normalization to compensate for this uncertainty.

## Export the actual CERN production sources

Point `--source` at the directory that contains the ROOT macro and Condor code,
not at the whole home directory. This tool collects source/config files only;
ROOT data, generated outputs and logs are excluded. It does not modify the source.

```bash
python3 tools/collect_production_sources.py plan \
  --source /absolute/path/to/your/production-code

python3 tools/collect_production_sources.py pack \
  --source /absolute/path/to/your/production-code \
  --output /eos/user/g/gkrintir/upsilon-production-sources.zip
```

If energy variants or job wrappers live in separate directories, make one archive
per directory. Include the environment setup script and describe which copy was
used for each energy/ID scenario. The archive contains a SHA256 inventory for a
precise comparison before promoting recovered code as the production baseline.
