#!/bin/bash

# Base directory where files will be processed
base_dir="/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/data"
# Directory for extracted files
extracted_dir="$base_dir/extracted"
# Directory for final file destination
final_dir="/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/raw-sumstats/MVP"
# Path to your CSV
csv_file="/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/data/MVP_pheno_list.csv"

# Create directories if they don't exist
mkdir -p "$base_dir" "$extracted_dir" "$final_dir"

# Function to move specific files to destination based on CSV data
move_files_to_destination() {
  local extracted_path="$1"
  local destination="$2"
  local current_tar="$3"
  
  echo "Extracted Path: $extracted_path"
  echo "Destination: $destination"
  echo "Current Tar: $current_tar"
  
  # Get the files associated with the current tar file from the associative array
  file_list="${files_to_move["$current_tar"]}"
  echo "Files to Move: $file_list"
  
  # Determine the subdirectory created by extraction
  subdir=$(find "$extracted_path" -mindepth 1 -maxdepth 1 -type d)
  echo "Subdirectory: $subdir"
  
  # Move each file to the destination directory
  for file in $file_list; do
    fullfile="$subdir/$file"
    if [ -f "$fullfile" ]; then
      mv "$fullfile" "$destination/"
      echo "Moved $file to $destination/"
    else
      echo "File $file not found in $subdir"
    fi
  done
}

# Declare the associative array to store files to move
declare -A files_to_move

# Read and parse the CSV file to populate the associative array
{
  read -r header # Skip the header row
  while IFS= read -r line; do
    # Use awk to properly parse CSV with quoted fields
    tar_file=$(echo "$line" | awk -v FPAT='([^,]+)|("[^"]+")' '{print $18}' | tr -d '"')
    file_name=$(echo "$line" | awk -v FPAT='([^,]+)|("[^"]+")' '{print $19}' | tr -d '"')
    
    # Skip rows with empty tar_file or file_name
    if [[ -z "$tar_file" || -z "$file_name" ]]; then
      echo "Skipping row with empty tar_file or file_name: $tar_file, $file_name"
      continue
    fi

    # Use the tar_file as a key in the associative array and append file_name
    if [[ -n "${files_to_move["$tar_file"]}" ]]; then
      files_to_move["$tar_file"]+=" $file_name"
    else
      files_to_move["$tar_file"]="$file_name"
    fi
  done
} < "$csv_file"


# URLs for one .tar and .md5 files (for testing one pair)
urls=(
  "https://ftp-ncbi-nlm-nih-gov.dartmouth.idm.oclc.org/dbgap/studies/phs002453/phs002453.v1.p1/analyses/GIA/phs002453.MVP_R4.1000G_AGR.GIA.Labs_batch1.analysis-PI.MULTI.tar"
  "https://ftp-ncbi-nlm-nih-gov.dartmouth.idm.oclc.org/dbgap/studies/phs002453/phs002453.v1.p1/analyses/GIA/phs002453.MVP_R4.1000G_AGR.GIA.Labs_batch1.analysis-PI.MULTI.tar.md5"
)

# Iterate over the URLs in pairs (only one pair today)
for ((i=0; i<${#urls[@]}; i+=2)); do
  tar_url=${urls[i]}
  md5_url=${urls[i+1]}
  tar_filename=$(basename "$tar_url")
  md5_filename=$(basename "$md5_url")
  
  # Download the .tar and .tar.md5 files
  wget -P "$base_dir" "$tar_url"
  wget -P "$base_dir" "$md5_url"
  
  cd "$base_dir" || exit
  
  # Verify the MD5 checksum
  md5sum -c "$md5_filename" &> /dev/null
  if [ $? -eq 0 ]; then
    echo "MD5 checksum for $tar_filename passed."
    
    # Extract the tar file
    tar -xvf "$tar_filename" -C "$extracted_dir"
    
    # Use the associative array to move files
    move_files_to_destination "$extracted_dir" "$final_dir" "$tar_filename"
    
    # Clean up if needed
    rm -f "$tar_filename" "$md5_filename"
    rm -rf "$extracted_dir"/*
  else
    echo "MD5 checksum for $tar_filename failed." >&2
  fi
  
  cd ..
done