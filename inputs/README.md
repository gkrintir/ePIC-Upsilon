# Paper inputs

The six merged ROOT files for the two-energy analysis note are hosted on [CERNBox](https://cernbox.cern.ch/s/UHPvC04xj29cxrp). The three 9×275 files supply the paper figures. `paper_inputs.csv` lists each path, byte count and SHA-256 digest, along with the generated row count and reconstruction modes.

From the repository root, choose a destination with at least 1.5 GiB free and run:

```sh
python3 inputs/fetch_paper_inputs.py work/paper-inputs
```

The script keeps the three state directories (`Upsilon1S`, `Upsilon2S` and `Upsilon3S`) under `work/paper-inputs`. Each file downloads in chunks to a temporary file in its destination directory and is installed only after its size and checksum match the manifest. Running the command again skips verified files and restarts incomplete downloads. If an existing file fails verification, move it aside before retrying. The destination filesystem must support hard links.

To check an existing copy entirely offline, run:

```sh
python3 inputs/fetch_paper_inputs.py work/paper-inputs --check-only
```

A nonzero exit status indicates that at least one file is missing, invalid or could not be downloaded.
