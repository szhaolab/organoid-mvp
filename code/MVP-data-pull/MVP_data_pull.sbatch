#!/bin/bash
#SBATCH --job-name=process_tar_files      # Name of the job
#SBATCH --output=job_output_%j.out        # Output file name (%j will be replaced with the job ID)
#SBATCH --error=job_error_%j.err          # Error file name (%j will be replaced with the job ID)
#SBATCH --time=24:00:00                   # Time limit (hh:mm:ss)
#SBATCH --mem=200G                        # Memory limit
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=END,FAIL              # Notifications for job done & fail
#SBATCH --mail-user=anne.k.hooker.gr@dartmouth.edu # Send-to email address

# Change to the directory where the bash script is located
cd /dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/code/MVP-data-pull

# Execute your bash script
bash MVP_data_pull_TEST.sh
