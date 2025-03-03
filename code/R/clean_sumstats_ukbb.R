library(tidyverse)

# Cleans summary statistics and normalizes them for downstream analysis
clean_sumstats <- function(sumstats, cols.to.keep){
  
  stopifnot(!is.null(sumstats))
  stopifnot(length(cols.to.keep) == 8)
  
  # Map input columns to expected column names
  chr <- cols.to.keep[1]  # CHROM
  pos <- cols.to.keep[2]  # POS
  es <- cols.to.keep[3]   # ES (effect size)
  se <- cols.to.keep[4]   # SE (standard error)
  a0 <- cols.to.keep[5]   # REF (reference allele)
  a1 <- cols.to.keep[6]   # ALT (alternate allele)
  rs <- cols.to.keep[7]   # ID (SNP ID)
  lp <- cols.to.keep[8]   # LP (-log10 p-value)
  
  # Extract relevant columns and rename them
  clean.sumstats <- sumstats[ ,c(chr, pos, es, se, a0, a1, rs, lp)]
  colnames(clean.sumstats) <- c('chr','pos','es','se','a0','a1','snp', 'lp')  # Lowercase es and lp
  
  # drop XY chromosomes
  clean.sumstats <- clean.sumstats[!(clean.sumstats$chr %in% c("X","Y")), ]
  
  # make chromosomes integers
  clean.sumstats$chr <- as.integer(clean.sumstats$chr)
  
  # Convert lp (-log10(pval)) to pval
  clean.sumstats$pval <- 10^(-clean.sumstats$lp)
  
  # Compute Z-scores using qnorm(pval/2) and match the sign of es
  zscore <- qnorm(clean.sumstats$pval / 2, lower.tail = FALSE)
  zscore <- ifelse(clean.sumstats$es < 0, -zscore, zscore)  # Match sign of es
  clean.sumstats['zscore'] <- zscore
  clean.sumstats <- clean.sumstats[!is.na(zscore),]
  
  # Keep SNPs only, no indels
  nucs <- c('A','C','T','G')
  bol <- (clean.sumstats$a0 %in% nucs) & (clean.sumstats$a1 %in% nucs)
  clean.sumstats <- clean.sumstats[bol,]

  # sort by chromosome and position
  clean.sumstats <- clean.sumstats[order(clean.sumstats$chr, clean.sumstats$pos), ]
  
  # drop duplicate SNPs
  chrpos <- paste0(clean.sumstats$chr, '_', clean.sumstats$pos)
  clean.sumstats <- clean.sumstats[!duplicated(chrpos), ]
    
  return(clean.sumstats)
}

args <- commandArgs(trailingOnly=T)
sumstats <- vroom::vroom(args[1], col_names = T)
cols.to.keep <- unlist(strsplit(as.character(args[2]), split=','))
outName <- args[3]

cleaned_sumstats <- clean_sumstats(sumstats, cols.to.keep)
vroom::vroom_write(cleaned_sumstats, outName)