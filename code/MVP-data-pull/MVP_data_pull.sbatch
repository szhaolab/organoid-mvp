#!/bin/bash
#SBATCH --job-name=process_tar_files      # Name of the job
#SBATCH --output=/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/code/MVP-data-pull/job_output_%j.out        # Output file name (%j will be replaced with the job ID)
#SBATCH --error=/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/code/MVP-data-pull/job_error_%j.err          # Error file name (%j will be replaced with the job ID)
#SBATCH --time=24:00:00                   # Time limit (hh:mm:ss)
#SBATCH --mem-per-cpu=20G                        # Memory limit
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL              # Notifications for job done & fail
#SBATCH --mail-user=anne.k.hooker.gr@dartmouth.edu # Send-to email address

# Change to the directory where the bash script is located
cd /dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/code/MVP-data-pull

# Execute your bash script
bash MVP_data_pull.sh
