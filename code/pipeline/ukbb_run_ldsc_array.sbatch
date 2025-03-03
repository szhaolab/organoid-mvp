#!/bin/bash -l
#SBATCH --job-name=ukbb_run_ldsc_array
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=1-24%20
#SBATCH --output=/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/code/pipeline/array-outputs/ukbb_run_ldsc_array_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL

################################
# Enter your code to run below #
# Usage: Run snakemake snake-ldsc-kh.smk
################################

eval "$(conda shell.bash hook)"
conda activate snakemake

# Hardcoded list of filenames
FILES=(
    "ukb-d-30600_irnt.vcf.gz"
    "ukb-d-30610_irnt.vcf.gz"
    "ukb-d-30620_irnt.vcf.gz"
    "ukb-d-30630_irnt.vcf.gz"
    "ukb-d-30640_irnt.vcf.gz"
    "ukb-d-30650_irnt.vcf.gz"
    "ukb-d-30670_irnt.vcf.gz"
    "ukb-d-30680_irnt.vcf.gz"
    "ukb-d-30690_irnt.vcf.gz"
    "ukb-d-30700_irnt.vcf.gz"
    "ukb-d-30710_irnt.vcf.gz"
    "ukb-d-30720_irnt.vcf.gz"
    "ukb-d-30730_irnt.vcf.gz"
    "ukb-d-30740_irnt.vcf.gz"
    "ukb-d-30750_irnt.vcf.gz"
    "ukb-d-30760_irnt.vcf.gz"
    "ukb-d-30770_irnt.vcf.gz"
    "ukb-d-30780_irnt.vcf.gz"
    "ukb-d-30810_irnt.vcf.gz"
    "ukb-d-30830_irnt.vcf.gz"
    "ukb-d-30840_irnt.vcf.gz"
    "ukb-d-30860_irnt.vcf.gz"
    "ukb-d-30870_irnt.vcf.gz"
    "ukb-d-30880_irnt.vcf.gz"
)

# Get the filename for this job array task
FILENAME="${FILES[$SLURM_ARRAY_TASK_ID-1]}"

# Extract the numeric trait identifier (e.g., 30600, 30610, etc.)
TRAIT=$(echo "$FILENAME" | sed -n 's/^ukb-d-\([0-9]\{5\}\)_irnt\.vcf\.gz$/\1/p')

# Run snakemake with the extracted trait
snakemake -p -c1 -s ukbb-snake-ldsc-kh.smk --config trait="$TRAIT"