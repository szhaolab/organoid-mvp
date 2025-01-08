suppressMessages(library(tidyverse))
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))

args <- commandArgs(trailingOnly = T)
bim <- args[1]
bed_dir <- args[2]
baseline <- args[3]
outName <- args[4]
MfileName <- args[5]

# Each annotation gets assigned SNPs based on overlap
annotator <- function(bim_file, annotations=NULL){
  
  if(is.null(annotations)){
    bim_file$weights <- 1
    return(bim_file)
  }
  
  snpRanges <- GRanges(seqnames = bim_file$CHR, 
                       ranges = IRanges(start = bim_file$BP, end = bim_file$BP))
  snpRanges <- plyranges::mutate(snpRanges, SNP=bim_file$SNP)
  
  for(f in annotations){
    
    name <- basename(f)
    name <- strsplit(name, split='[.]')[[1]][1]
    curr <- import(f, format='bed')
    subdf <- subsetByOverlaps(snpRanges, curr)
    snpsIn <- unique(subdf$SNP)
    bim_file <- mutate(bim_file, !!name := ifelse(SNP %in% snpsIn,1,0))
  }
  return(bim_file)
}

annot <- vroom::vroom(bim, delim = '\t', col_names = F)
annot <- annot[ ,c(2,1,4,5,6)] 
colnames(annot) <- c("SNP", "CHR","BP","A1","A2")

if(dir.exists(bed_dir)){
  bed_files <- list.files(bed_dir, pattern = '*.bed', full.names = T)
} else{
  bed_files <- NULL
}
# if(file.exists(bed_dir)){
#   bed_files = bed_dir
# } else{
#   bed_files = NULL
# }

base_annot <- annotator(annot, bed_files)

baseline_annots <- vroom::vroom(baseline)
base_annot <- inner_join(base_annot, baseline_annots, by=c('SNP','CHR','BP','A1','A2'))
vroom::vroom_write(base_annot, outName)

M_file <- unname(colSums(base_annot[, c(6:ncol(base_annot))])) 
writeLines(paste0(M_file, collapse = ' '), con = MfileName)


