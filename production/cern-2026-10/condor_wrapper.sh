#!/bin/bash
INPUT_LIST=$1
OUTPUT_DIR=$2

echo "Starting job on $(hostname) at $(date)"
mkdir -p "$OUTPUT_DIR"

# HTCondor transferred the whole folder, so we enter it
cd lAger_JPsi_MuMu_Dataset/
chmod +x run_Upsilon_Analysis_new_condor.sh

# list.txt was transferred to the top-level scratch directory, so we reference it with ../
./run_Upsilon_Analysis_new_condor.sh ../$INPUT_LIST $OUTPUT_DIR

echo "Job finished at $(date)"
