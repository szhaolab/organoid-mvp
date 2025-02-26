#!/bin/bash -l
#SBATCH --job-name=run_ldsc_array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-97%20
#SBATCH --output=/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/code/pipeline/array-outputs/run_ldsc_array_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL

################################
# Enter your code to run below #
# Usage: Run snakemake snake-ldsc-kh.smk
################################

eval "$(conda shell.bash hook)"
conda activate snakemake

# Raw sumstats directory
INPUT_DIR="/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/raw-sumstats"
FILES=("$INPUT_DIR"/*.txt.gz) # Ignore the metadata files
FILENAME="${FILES[$SLURM_ARRAY_TASK_ID-1]}" # Adjust for zero indexing of the job array
BASENAME=$(basename "$FILENAME" .txt.gz)  # Remove the .txt.gz extension
# Extract the relevant part of the filename
NAME_PART=${BASENAME#MVP_R4.1000G_AGR.}  # Remove the prefix
NAME_PART=${NAME_PART%%.EUR.GIA.dbGaP}

snakemake -p -c1 -s snake-ldsc-kh.smk --config trait="$NAME_PART"

