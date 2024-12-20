# Rotation Notes

## Coding Notes

### Dry run of snakemake on the login node
    snakemake -np -s /dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/snake-ldsc-kh.smk
-n → dry run\
-p → prints the shell commands that will be executed as part of the workflow

### Set up interactive mode and run the snakemake on a designated compute node
    srun --nodes=2 --ntasks-per-node=4 --mem-per-cpu=1GB --cpus-per-task=1 --time=03:00:00 --pty /bin/bash
    snakemake -s /dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/snake-ldsc-kh.smk


## Meetings with Siming

### 12/20/2024
- Liyang's pipeline uses the hg19 build, need to make sure LDSC agrees with this
- How does the UK biobank define the effect allele vs. reference allele (a0/a1, major/minor)?
    - This will impact the sign of the beta value and z-score
    - Need to make sure this agrees with the clean_sumstats.R script and the LDSC pipeline
- Walk through what the snakemake pipeline is doing, step by step
- Copy the clean_sumstats.R script from cancer_polyfun 
- Will need to filter out some already-irrelevant phenotypes from the UK biobank that we KNOW are irrelevant to colon epithelial cells
    - Should include all digestive traits at minimum
    - Can reasonably test a couple hundred traits
- Rmd (workflowr package) or JupyterNotebook to document results