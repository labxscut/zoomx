#!/usr/bin/env Rscript
version="REPLACE_WITH_COMMIT_OR_VERSION"
options(warn=-1)

#source("https://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
library(optparse)
library(GenomicRanges)
library(rtracklayer)

option_list <- list(
)

parser <- OptionParser(
  usage="%prog [options] gridFile.bed fragFile.bed overlap.mtx.gz \n producing grid-frag overlapping matrix file", 
  option_list=option_list)

cat("-Info: vesion:",version,"\n")
cat("-Info: invoking command:",commandArgs(),"\n")

args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, 
                 positional_arguments = TRUE)
if(length(cmd$args)==3){ 
  gridFile=cmd$args[1]
  fragFile=cmd$args[2]
  outputFile=cmd$args[3]
} else {
  print_help(parser);
  gridFile="~/ws/hg/hg19/human_g1k_v37.ucsc.fasta.bed.w100Ks50K"
  fragFile="P04979_CR1.wg10.bam.fragBar.bed.a50000"
  outputFile="P04979_CR1.wg10.bam.grid100K.overlap.mtx.gz"
}

cat("gridFile:", gridFile, "\n")
cat("fragFile:", fragFile, "\n")
cat("outputFile:", outputFile, "\n")

grid <- import(con=gridFile,format="bed")
cat("gridSize:", length(grid), "\n")
fragdata <- read.table(fragFile,header=F)
cat("fragdataSize:", nrow(fragdata), "\n")
colnames(fragdata) <- c('chr','start','end','bar','fragsize','fraglength','tagfragfreq')
frags <- with(fragdata, GRanges(chr, IRanges(start+1, end)))
overlap <- as.matrix(findOverlaps(grid,frags))  #findOverlaps(query,subject) -> [gidx, fidx]
cat("overlapSize:", nrow(overlap), "\n")

gz <- gzfile(outputFile, "w")
write.table(overlap, file=gz, quote=F, row.names=F, sep="\t")
close(gz)
cat("gridOverlapDONE\n")
