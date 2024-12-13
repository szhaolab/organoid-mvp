# A modified snakemake protocol for performing enrichment analysis and functionally-informed finemapping
# Uses PolyFun

import os
import glob
import pandas as pd

# Below creates log directory but not used?
logdir = "log/"
if not os.path.isdir(logdir):
  os.mkdir(logdir)

pdnew = "/dartfs/rc/lab/S/Szhao/liyang/enrichment/ARID1A_proj/" # contain R scripts for rules
pdorig = "/dartfs/rc/lab/S/Szhao/fine-mapping/cancer-polyfun/" # contain required metadata 
polyfun = '/dartfs/rc/lab/S/Szhao/katieh/polyfun/' # contain polyfun software (download from github)
pdout = "/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/outputs/" # replace with desired output folder

# Inputs (summary stats, a directory of annotations)
cleaned_sumstats = "/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/cleaned_sumstats/" 
bed_dir = "/dartfs/rc/lab/S/Szhao/katieh/organoid-rotation/bed_dir/"

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
traits = "IBD3"
chrom = list(map(str, range(1,23)))


rule all:
  input:
      expand(results + '{t}_enrichment.results', t=traits),
      expand(results + '{t}_signif_loci.txt', t=traits),
      # expand(results + "{t}/{t}.{c}.snpvar_ridge.gz",t=traits, c=chrom),
      # expand(results + "{t}/{t}.{c}.snpvar_ridge_constrained.gz", t=traits, c=chrom),
      # expand(finemapping + 'processed/{t}_finemapped_susie_L1.txt.gz',t=traits)

rule munge_sumstats:
    input:
       cleaned_sumstats + '{traits}_sumstats.txt.gz'
    output:
       munged_sumstats +  '{traits}_munged_sumstats.parquet'
    params:
       n=60000
    shell:
       "python {polyfun}munge_polyfun_sumstats.py --sumstats {input} --n {params.n} --out {output}"

rule get_significant_loci:
    input:
      cleaned_sumstats + '{traits}_sumstats.txt.gz',
      ldblocks
    output:
      results + '{traits}_signif_loci.txt'
    shell:
      "Rscript {pdnew}get_signif_loci_ly.R {input[0]} {input[1]} {output}"

rule create_annotations:
    input:
      onekg + '.{chrom}.bim',
      bed_dir,
      baseline + 'baseline_MAF_LD.{chrom}.annot.gz'
    output:
      annotations + '{traits}/{traits}.{chrom}.annot.gz',
      annotations + '{traits}/{traits}.{chrom}.l2.M'
    shell:
      "Rscript {pdnew}create_annotations_LY.R {input[0]} {input[1]} {input[2]} {output[0]} {output[1]}"

rule compute_ld_scores:
    input:
        annotations + '{traits}/{traits}.{chrom}.annot.gz'
    output:
        annotations + '{traits}/{traits}.{chrom}.l2.ldscore.parquet'
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
      munged_sumstats + '{traits}_munged_sumstats.parquet',
      expand(annotations + '{{traits}}/{{traits}}.{c}.l2.ldscore.parquet', c=chrom)
    output:
      results + '{traits}_enrichment.results'
    params:
      annot_prefix = annotations + '{traits}/{traits}.',
      out_prefix = results + '{traits}_enrichment'
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

# rule polyfun_ldsc:
#     input:
#        munged_sumstats +  '{traits}_munged_sumstats.parquet',
#        expand(annotations + '{{traits}}/{{traits}}.{c}.l2.ldscore.parquet', c=chrom)
#     output:
#        expand(results + "{{traits}}/{{traits}}.{c}.snpvar_ridge.gz", c=chrom),
#        expand(results + "{{traits}}/{{traits}}.{c}.snpvar_ridge_constrained.gz", c=chrom)
#     params:
#        ld_prefix = annotations + "{traits}/{traits}.",
#        out_prefix = results + "{traits}/{traits}"
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
#       expand(results + "{{traits}}/{{traits}}.{c}.snpvar_ridge_constrained.gz", c=chrom),
#       results + '{traits}_signif_loci.txt'
#     output:
#       finemapping + 'processed/{traits}_finemapped_susie_L1.txt.gz'
#     params:
#       trait = "{traits}"
#     shell:
#       "/dartfs/rc/lab/S/Szhao/liyang/organoid_stimulus_proj/ATAC/Combined_samples_202405/filtered/bam2vcf/UT/ASoC_samp18/polyfun/finemap-ly.sh {params.trait}"
