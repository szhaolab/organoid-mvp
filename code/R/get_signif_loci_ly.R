
suppressMessages(library(tidyverse))
suppressMessages(library(GenomicRanges))

make_ranges <- function(seqname, start, end){
  return(GRanges(seqnames = seqname, ranges = IRanges(start = start, end = end)))
}

assign.locus.snp <- function(cleaned.sumstats, ldBed){
  
  cleaned.sumstats <- cleaned.sumstats %>% drop_na()
  cleaned.sumstats <- cleaned.sumstats[cleaned.sumstats$pval > 0, ]
  
  ldRanges <- make_ranges(ldBed$X1, ldBed$X2, ldBed$X3)
  ldRanges <- plyranges::mutate(ldRanges, locus=ldBed$X4)
  
  snpRanges <- GRanges(seqnames = cleaned.sumstats$chr, 
                       ranges   = IRanges(start = cleaned.sumstats$pos, 
                                          end   = cleaned.sumstats$pos,
                                          names = cleaned.sumstats$snp))
  
  snpRanges <- plyranges::mutate(snpRanges, snp=names(snpRanges))
  
  snp.ld.overlap <- plyranges::join_overlap_inner(snpRanges, ldRanges)
  snp.ld.block <- as_tibble(snp.ld.overlap@elementMetadata)
  snp.ld.block <- snp.ld.block[!duplicated(snp.ld.block$snp), ] # some SNPs are in multiple ld-blocks due to edge of ld blocks
  cleaned.annot.sumstats <- inner_join(cleaned.sumstats, snp.ld.block, 'snp')
  
  return(cleaned.annot.sumstats)
}


args <- commandArgs(trailingOnly = T)

sumstats <- vroom::vroom(args[1], col_names = T)
ldblocks <- vroom::vroom(args[2], col_names = F)
outName <- args[3]

cleaned.sumstats <- assign.locus.snp(sumstats, ldblocks) 
signif_loci <- cleaned.sumstats %>% group_by(locus) %>% summarise(n_signif=sum(pval<=5e-8)) %>% filter(n_signif > 0)

if(nrow(signif_loci)>150){
  signif_loci <- signif_loci[order(signif_loci$n_signif, decreasing = T), ]
  top.loci <- signif_loci$locus[1:150]
} else{
  top.loci <- signif_loci$locus
}

sig_loci_df <- ldblocks[ldblocks$X4 %in% top.loci, ]

vroom::vroom_write(sig_loci_df, outName, col_names = F)

