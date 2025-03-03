#!/bin/bash

# Input parameters
input_file=$1      # Full path to the input VCF file (gzipped)
output_file=$2     # Full path to the output file (gzipped)

# Check if the input file exists
if [[ ! -f $input_file ]]; then
    echo "Error: Input file $input_file not found!"
    exit 1
fi

# Process the data and write the compressed output
{
    # Print the header for the output file (without "#" in CHROM)
    echo -e "CHROM\tPOS\tID\tREF\tALT\tES\tSE\tLP\tSS"

    # Process the data lines
    zcat "$input_file" | awk '!/^#/ {
        # Split the FORMAT column (9th column) into individual fields
        split($9, format, ":");
        # Split the sample column (10th column) into individual values
        split($10, values, ":");
        # Extract the required fields
        chrom = $1;
        pos = $2;
        id = $3;
        ref = $4;
        alt = $5;
        es = values[1];
        se = values[2];
        lp = values[3];
        ss = values[5];  # Sample size (SS) from the sample column
        # Print the extracted fields
        print chrom "\t" pos "\t" id "\t" ref "\t" alt "\t" es "\t" se "\t" lp "\t" ss;
    }'
} | gzip > "$output_file"

echo "Processed file saved as $output_file"