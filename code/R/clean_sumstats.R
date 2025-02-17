
library(tidyverse)

# Cleans summary statistics and normalizes them for downstream analysis
clean_sumstats <- function(sumstats, cols.to.keep){
  
  stopifnot(!is.null(sumstats))
  stopifnot(length(cols.to.keep) == 8)
  
  chr <- cols.to.keep[1]
  pos <- cols.to.keep[2]
  beta <- cols.to.keep[3]
  se <- cols.to.keep[4]
  a0 <- cols.to.keep[5]
  a1 <- cols.to.keep[6]
  rs <- cols.to.keep[7]
  pval <- cols.to.keep[8]
  
  # keep SNPs in 1kg
  #sumstats <- inner_join(sumstats, snps.to.keep, by=rs)
  # Extract relevant columns
  clean.sumstats <- sumstats[ ,c(chr, pos, beta, se, a0, a1, rs, pval)]
  colnames(clean.sumstats) <- c('chr','pos','beta','se','a0','a1','snp', 'pval')
  
  # drop XY chromosomes
  clean.sumstats <- clean.sumstats[!(clean.sumstats$chr %in% c("X","Y")), ]
  
  # make chromosomes integers
  clean.sumstats$chr <- as.integer(clean.sumstats$chr)
  
  # Compute Zscores
  zscore <- clean.sumstats$beta/clean.sumstats$se
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


