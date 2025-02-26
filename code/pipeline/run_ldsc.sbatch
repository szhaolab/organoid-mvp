#!/bin/bash -l
#SBATCH --job-name=run_ldsc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=run_ldsc_%j.log
#SBATCH --mail-type=BEGIN,END,FAIL

################################
# Enter your code to run below #
# Usage: Run snakemake snake-ldsc-kh.smk
################################

eval "$(conda shell.bash hook)"
conda activate snakemake

input_smk=$1

# output_folder=$1
# cd  ${output_folder} 

snakemake -p -c1 -s $input_smk --config trait="Phe_578_8"

