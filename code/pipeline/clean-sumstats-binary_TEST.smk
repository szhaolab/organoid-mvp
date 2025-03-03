# A test script for the clean_sumstats rule of the snakemake protocol

import os
import glob
import pandas as pd
import re

# Below creates log directory but not used?
logdir = "log/"
if not os.path.isdir(logdir):
  os.mkdir(logdir)

pdnew = "/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/code/R/" # contain R scripts for rules
pdorig = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/" # contain required metadata 
polyfun = '/dartfs/rc/lab/S/Szhao/katieh/polyfun/' # contain polyfun software (download from github)
pdout = "/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/outputs/MVP/" # desired output folder
pdclean = "/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/code/R/" # contains the clean_sumstats R script

# Inputs (summary stats, a directory of annotations)
raw_sumstats = "/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/raw-sumstats/MVP/"
cleaned_sumstats = "/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/cleaned-sumstats/MVP/" 
bed_dir = "/dartfs/rc/lab/S/Szhao/fine-mapping/organoid-polyfun/bed_dir/" # Where the ASoC annotations are stored

# REQUIRED metadata
onekg = pdorig + '1000G_EUR_Phase3_plink/1000G.EUR.QC'
baseline = pdorig + 'Gazel_LD_baseline/' # should be referring to Gazel paper
weights = pdorig + 'weights/weights.'
ldblocks = pdorig + 'Euro_LD_Chunks.bed',

# Outputs created
results = pdout + 'results/'
annotations = pdout + 'annotations/'
munged_sumstats = pdout + 'munged_sumstats/'
#finemapping = pdout + 'finemapping/'

# Global wildcards
traits = "Phe_578_8"

chrom = list(map(str, range(1,23)))

# Cleaned sumstats columns
COLUMNS = "chrom,pos,beta,sebeta,ref,alt,SNP_ID,pval"
BINARY_COLUMNS = "chrom,pos,or,ci,ref,alt,SNP_ID,pval"

rule all:
  input:
      expand(cleaned_sumstats + '{t}_sumstats.txt.gz', t=traits)

rule clean_sumstats_binary:
    input:
        raw_sumstats + 'MVP_R4.1000G_AGR.{binary_traits}.EUR.GIA.dbGaP.txt.gz'
    output:
        cleaned_sumstats + '{binary_traits}_sumstats.txt.gz'
    params:
      columns = BINARY_COLUMNS
    shell:
        "Rscript {pdclean}clean_sumstats_binary.R {input[0]} {params.columns} {output}"