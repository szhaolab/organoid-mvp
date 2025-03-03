# A modified snakemake protocol for performing enrichment analysis and functionally-informed finemapping
# Uses PolyFun

import os
import glob
import pandas as pd
import re

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
      expand(results + '{t}_enrichment.results', t=traits),
      expand(results + '{t}_signif_loci.txt', t=traits),
      expand(cleaned_sumstats + '{t}_sumstats.txt.gz', t=traits)
      # expand(results + "{t}/{t}.{c}.snpvar_ridge.gz",t=traits, c=chrom),
      # expand(results + "{t}/{t}.{c}.snpvar_ridge_constrained.gz", t=traits, c=chrom),
      # expand(finemapping + 'processed/{t}_finemapped_susie_L1.txt.gz',t=traits)

rule clean_sumstats:
    input:
        raw_sumstats + 'MVP_R4.1000G_AGR.{trait}.EUR.GIA.dbGaP.txt.gz'
    output:
        cleaned_sumstats + '{trait}_sumstats.txt.gz'
    params:
        # Determine which R script and columns to use based on the trait name
        script = lambda wildcards: (
            f"{pdclean}clean_sumstats.R" if "INT" in wildcards.trait
            else f"{pdclean}clean_sumstats_binary.R"
        ),
        columns = lambda wildcards: (
            COLUMNS if "INT" in wildcards.trait
            else BINARY_COLUMNS
        )
    shell:
        """
        Rscript {params.script} {input[0]} {params.columns} {output}
        """

rule munge_sumstats:
    input:
       cleaned_sumstats + '{trait}_sumstats.txt.gz'
    output:
       munged_sumstats +  '{trait}_munged_sumstats.parquet'
    params:
       n=60000
    shell:
       "python {polyfun}munge_polyfun_sumstats.py --sumstats {input} --n {params.n} --out {output}"

rule get_significant_loci:
    input:
      cleaned_sumstats + '{trait}_sumstats.txt.gz',
      ldblocks
    output:
      results + '{trait}_signif_loci.txt'
    shell:
      "Rscript {pdnew}get_signif_loci_KH.R {input[0]} {input[1]} {output}"

rule create_annotations:
    input:
      onekg + '.{chrom}.bim',
      bed_dir,
      baseline + 'baseline_MAF_LD.{chrom}.annot.gz'
    output:
      annotations + '{trait}/{trait}.{chrom}.annot.gz',
      annotations + '{trait}/{trait}.{chrom}.l2.M'
    shell:
      "Rscript {pdnew}create_annotations_KH.R {input[0]} {input[1]} {input[2]} {output[0]} {output[1]}"

# Compute LD scores for the AsoC annotation from a reference panel
rule compute_ld_scores:
    input:
        annotations + '{trait}/{trait}.{chrom}.annot.gz'
    output:
        annotations + '{trait}/{trait}.{chrom}.l2.ldscore.parquet'
    params:
        prefix = onekg + '.{chrom}',
    shell:
        """
        python {polyfun}compute_ldscores.py \
               --bfile {params.prefix} \
               --annot {input} \
               --out {output}
        """

rule enrichment_ldsc:
    input:
      munged_sumstats + '{trait}_munged_sumstats.parquet',
      expand(annotations + '{{trait}}/{{trait}}.{c}.l2.ldscore.parquet', c=chrom)
    output:
      results + '{trait}_enrichment.results'
    params:
      annot_prefix = annotations + '{trait}/{trait}.',
      out_prefix = results + '{trait}_enrichment'
    shell:
      """
      python {polyfun}ldsc.py --h2 {input[0]} \
             --ref-ld-chr {params.annot_prefix} \
             --w-ld-chr {weights} \
             --not-M-5-50 \
             --out {params.out_prefix} \
             --overlap-annot \
      && sed -i 's/_0//g' {output}
      """

# The not-M-5-50 argument specifies estimating functional enrichment of the heritability causally explained by all SNPs (not just common MAF < 5%)





# rule polyfun_ldsc:
#     input:
#        munged_sumstats +  '{trait}_munged_sumstats.parquet',
#        expand(annotations + '{{trait}}/{{trait}}.{c}.l2.ldscore.parquet', c=chrom)
#     output:
#        expand(results + "{{trait}}/{{trait}}.{c}.snpvar_ridge.gz", c=chrom),
#        expand(results + "{{trait}}/{{trait}}.{c}.snpvar_ridge_constrained.gz", c=chrom)
#     params:
#        ld_prefix = annotations + "{trait}/{trait}.",
#        out_prefix = results + "{trait}/{trait}"
#     shell:
#         """
#         python {polyfun}polyfun.py \
#                --compute-h2-L2 \
#                --no-partitions \
#                --output-prefix {params.out_prefix} \
#                --sumstats {input[0]} \
#                --ref-ld-chr {params.ld_prefix} \
#                --w-ld-chr {weights} \
#                --allow-missing
#          """

# rule run_finemapping:
#     input:
#       expand(results + "{{trait}}/{{trait}}.{c}.snpvar_ridge_constrained.gz", c=chrom),
#       results + '{trait}_signif_loci.txt'
#     output:
#       finemapping + 'processed/{trait}_finemapped_susie_L1.txt.gz'
#     params:
#       trait = "{trait}"
#     shell:
#       "/dartfs/rc/lab/S/Szhao/liyang/organoid_stimulus_proj/ATAC/Combined_samples_202405/filtered/bam2vcf/UT/ASoC_samp18/polyfun/finemap-ly.sh {params.trait}"
