#!/bin/bash

# Define the output file
output_file="combined_results.txt"

# Loop through all files ending with _enrichment.results
for file in *_enrichment.results; do
    # Append the file name to the output file
    echo "File: $file" >> "$output_file"
    
    # Append the contents of the file to the output file
    cat "$file" >> "$output_file"
    
    # Add a separator between files for better readability
    echo -e "\n---\n" >> "$output_file"
done

echo "Results have been saved to $output_file"
