#!/bin/bash

INPUT_LIST=$1
OUTPUT_DIR=$2

# 1. Check if input is a valid file (changed from -d to -f)
if [ ! -f "$INPUT_LIST" ]; then
    echo "Error: Input file '$INPUT_LIST' does not exist."
    exit 1
fi

# 2. Check if output directory is provided and create it if missing
if [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Output directory not provided."
    echo "Usage: ./run_Upsilon_Analysis_new.sh <input_list.txt> <output_dir>"
    exit 1
fi
mkdir -p "$OUTPUT_DIR"

# 3. Read the text file line by line
while IFS= read -r filepath || [ -n "$filepath" ]; do
    # Skip empty lines or lines starting with #
    [[ -z "$filepath" || "$filepath" == \#* ]] && continue

    # Extract the base filename
    filename=$(basename "$filepath")
    
    # Sanitize the filename to prevent ROOT errors (replacing dots with underscores)
    name_no_ext="${filename%.root}"
    safe_name=$(echo "$name_no_ext" | tr '.' '_')
    
    # Define the output path
    output_file="${OUTPUT_DIR}/AnalysisOut_${safe_name}.root"

    echo "=========================================================="
    echo "Processing: $filepath"
    echo "Output: $output_file"
    
    # Run the compiled analysis macro 
    # (Adjust the macro name or arguments if your DVMP_Upsilon_Analysis.C signature differs)
    root -l -b -q "DVMP_Upsilon_Analysis.C+(\"${filepath}\", \"${output_file}\")"

done < "$INPUT_LIST"

echo "All files in $INPUT_LIST processed successfully!"
