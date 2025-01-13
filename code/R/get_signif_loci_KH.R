
# Copied and annotated from /fine-mapping/organoid-polyfun/get_signif_loci_ly.R
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))

# Creates a GRanges object
make_ranges <- function(seqname, start, end){
  return(GRanges(seqnames = seqname, ranges = IRanges(start = start, end = end)))
}

assign.locus.snp <- function(cleaned.sumstats, ldBed){
  
  cleaned.sumstats <- cleaned.sumstats %>% drop_na() # Removing rows with missing values
  cleaned.sumstats <- cleaned.sumstats[cleaned.sumstats$pval > 0, ] # Filtering out SNPs with non positive p-vals
  
  ldRanges <- make_ranges(ldBed$X1, ldBed$X2, ldBed$X3) # Creates GRanges for LD blocks based on ldBed coordinates
  ldRanges <- plyranges::mutate(ldRanges, locus=ldBed$X4) # Adds locus column (LD block identifier) to ldRanges
  
  snpRanges <- GRanges(seqnames = cleaned.sumstats$chr, 
                       ranges   = IRanges(start = cleaned.sumstats$pos, 
                                          end   = cleaned.sumstats$pos,
                                          names = cleaned.sumstats$snp)) # Creates genomic ranges for the SNPs
  
  snpRanges <- plyranges::mutate(snpRanges, snp=names(snpRanges)) # Adds SNP identifiers
  
  snp.ld.overlap <- plyranges::join_overlap_inner(snpRanges, ldRanges) # Find overlap and return SNPs that fall within the LD block regions
  snp.ld.block <- as_tibble(snp.ld.overlap@elementMetadata)
  snp.ld.block <- snp.ld.block[!duplicated(snp.ld.block$snp), ] # some SNPs are in multiple ld-blocks due to edge of ld blocks
  cleaned.annot.sumstats <- inner_join(cleaned.sumstats, snp.ld.block, 'snp') # Merges clean sum stats with the LD block info
  
  return(cleaned.annot.sumstats)
}


args <- commandArgs(trailingOnly = T)

sumstats <- vroom::vroom(args[1], col_names = T)
ldblocks <- vroom::vroom(args[2], col_names = F)
outName <- args[3]

cleaned.sumstats <- assign.locus.snp(sumstats, ldblocks) 
signif_loci <- cleaned.sumstats %>% group_by(locus) %>% summarise(n_signif=sum(pval<=5e-8)) %>% filter(n_signif > 0) # Group by locus, count number of significant SNPs in the locus, and retain only loci with at least one significant SNP 

# If more than 150 significant loci, sort them by number of significant SNPs and keep only the top 150
if(nrow(signif_loci)>150){
  signif_loci <- signif_loci[order(signif_loci$n_signif, decreasing = T), ]
  top.loci <- signif_loci$locus[1:150]
} else{
  top.loci <- signif_loci$locus
}

sig_loci_df <- ldblocks[ldblocks$X4 %in% top.loci, ]

vroom::vroom_write(sig_loci_df, outName, col_names = F)

