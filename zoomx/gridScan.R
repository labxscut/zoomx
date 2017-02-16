#!/usr/bin/env Rscript
version="REPLACE_WITH_COMMIT_OR_VERSION"
options(warn=-1)

#source("https://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
library(optparse)
library(GenomicRanges)
library(rtracklayer)
library(rhdf5)
library(hash)
library(Matrix)
library(plyr)

option_list <- list(
  #make_option(c("-o", "--output"), default="gridScan", help="output prefix, [default: %default] \n"),
  #make_option(c("-h", "--hdf5File"), default="none", help="hdf5 file (output from parse10x.py), [default: %default] \n"),
  #make_option(c("-f", "--fragFile"), default="none", help="fragment file (BED formatted: chr, start, end, barcode, fragsize, barfragcount), [default: %default] \n"),
  #make_option(c("-g1", "--gridFile1"), default="none", help="grid file region 1 (BED formatted: chr, start, end), [default: %default] \n"),
  #make_option(c("-g2", "--gridFile2"), default="none", help="grid file region 2 (BED formatted: chr, start, end), [default: %default] \n"),
  make_option(c("-mtf", "--maxTagFraction"), default="none", help="max tag fraction to be considerrred, [default: %default] \n"),
  make_option(c("-efd", "--expectFalseDiscover"), default=1000, help="expected false discoveries, [default: %default] \n"),
  make_option(c("-G", "--genomeSize"), default=3000000000, help="bp size of the whole genome, [default: %default] \n"),
  make_option(c("-R", "--regionSize"), default=100000, help="bp size of the whole genome, [default: %default] \n") 
)

parser <- OptionParser(
  usage="%prog [options] gridFile1.bed gridFile2.bed fragFile hdf5File outputFile \n Screening potential grid pairs by Poisson scan", 
  option_list=option_list)

cat("-Info: vesion:",version,"\n")
cat("-Info: invoking command:",commandArgs(),"\n")

args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, 
                 positional_arguments = TRUE)
if(length(cmd$args)==5){ 
  gridFile1=cmd$args[1]
  gridFile2=cmd$args[2]
  fragFile=cmd$args[3]
  hdf5File=cmd$args[4]
  outputFile=cmd$args[5]
  maxTagFraction=cmd$options$maxTagFraction
  expectFalseDiscover=cmd$options$expectFalseDiscover
  genomeSize=cmd$options$genomeSize
  regionSize=cmd$options$regionSize
} else {
  print_help(parser);
  setwd("~/setup/genomics/10x")
  gridFile1="~/ws/hg/hg19/human_g1k_v37.ucsc.fasta.bed.w100Ks50K"
  gridFile2="~/ws/hg/hg19/human_g1k_v37.ucsc.fasta.bed.w100Ks50K"
  fragFile="P04979_CR1.wg10.bam.fragBar.bed.a50000"
  hdf5File="P04979_CR1.wg10.bam.a50000.h5"
  outputFile="P04979_CR1.wg10.bam.grid100K.txt"
  genomeSize=3000000000   #size of genome 3G
  regionSize=100000       #size of the region 100K
  expectFalseDiscover=1000 #tolerated false discoveries
}

cat("gridFile1:", gridFile1, "\n")
cat("gridFile2:", gridFile2, "\n")
cat("fragFile:", fragFile, "\n")
cat("hdf5File:", hdf5File, "\n")
cat("outputFile:", outputFile, "\n")
cat("genomeSize:", genomeSize, "\n")
cat("regionSize:", regionSize, "\n")
cat("expectFalseDiscover:", expectFalseDiscover, "\n")

grid1 <- import(con=gridFile1,format="bed")
cat("grid1Size:", length(grid1), "\n")
grid2 <- import(con=gridFile2,format="bed")
cat("grid2Size:", length(grid2), "\n")
fragdata <- read.table(fragFile,header=F)
cat("fragdataSize:", nrow(fragdata), "\n")
colnames(fragdata) <- c('chr','start','end','bar','fragsize','fragcount')
frags <- with(fragdata, GRanges(chr, IRanges(start+1, end)))
overlap1 <- as.matrix(findOverlaps(grid1,frags))  #findOverlaps(query,subject) -> [gidx, fidx]
cat("overlap1Size:", nrow(overlap1), "\n")
overlap2 <- as.matrix(findOverlaps(grid2,frags))  #findOverlaps(query,subject) -> [gidx, fidx]
cat("overlap2Size:", nrow(overlap2), "\n")

#requirement hdf5 to contain ALL possible Tb's
#fragFile to caontain ALL Fragments

#DONE: modify parse10x.py to generate necessary statistics: N_f, N_b, mean_f, u0=N_f*mean(L_f)/G, N_f^2f^2sum(Tb)^2
regionFraction=regionSize/genomeSize              #fraction of the size of a region compared to a Gennome
cat("regionFraction:", regionFraction, "\n")
suppressMessages({tagFractions=as.data.frame(h5read(hdf5File,"T",bit64conversion='double'))})
cat("tagFractions:\n")
summary(tagFractions[4])
suppressMessages({fragNumber=as.numeric(as.data.frame(h5read(hdf5File,"N_f",bit64conversion='double'))[4])})
cat("fragNumber:", fragNumber, "\n")
suppressMessages({averageFragLength=as.numeric(as.data.frame(h5read(hdf5File,"L_f",bit64conversion='double'))[4])})        
#size of a fragment in bases
cat("averageFragLength:", averageFragLength, "\n")
tagFractionSquareSum=sum(tagFractions[4]*tagFractions[4])
cat("tagFractionSquareSum:", tagFractionSquareSum, "\n")
expectShareFrag=(fragNumber**2)*(regionFraction**2)*tagFractionSquareSum      #u0, expected shared Fragment by random chance
cat("expectShareFrag:", expectShareFrag, "\n")
expectJunctionFrag=fragNumber*averageFragLength/genomeSize                    #u, expected shared Fragment by junction formation
cat("expectJunctionFrag:", expectJunctionFrag, "\n")

#DONE: create a hash table Tb->Tb_index; NOTE, we need to cutoff excessive Tb here
cat("creating hash table for tag and fraction...")
if(cmd$options$maxTagFraction == "none") maxTagFraction=max(quantile(tagFractions[,4],0.99),0.001) else maxTagFraction=as.numeric(maxTagFraction)
#for(i in seq_len(nrow(tagFractions))) { print(tagFractions[which(tagFractions$axis0==toString(fragdata$bar[overlaps[i,2]])),]) }
system.time({
  tidx2tag=rep(0,nrow(tagFractions))
  tag2tidx=rep(0,nrow(tagFractions))
  tidx=1
  for(i in seq_len(nrow(tagFractions))){
    if(tagFractions[i,4]<=maxTagFraction){
      tidx2tag[tidx]=i
      tag2tidx[i]=tidx
      tidx=tidx+1
    } #reducing working barcodes by filtering out over abundant ones; tag is the original index, tidx is index in computation
  }
})
tidx2tag=tidx2tag[tidx2tag!=0]
cat("done\n")

#DEBUG: testing limit of system
#  library(plyr)
#  library(Matrix)
#  X_i=rep(1:61928,630)
#  X_j=sample(1:554715,length(X_i),replace=TRUE)
#  X_data=data.frame(i=X_i,j=X_j)
#  X_data=X_data[complete.cases(X_data),]
#  X_data=count(X_data, vars=c("i","j"))
#  X=sparseMatrix(i=X_data$i, j=X_data$j, x=X_data$freq)
#  Y_i=X_j
#  Y_j=X_i
#  Y_data=data.frame(i=Y_i,j=Y_j)
#  Y_data=Y_data[complete.cases(Y_data),]
#  Y_data=count(Y_data, vars=c("i","j"))
#  Y_data$freq=(Y_data$freq>=1)
#  Y=sparseMatrix(i=Y_data$i, j=Y_data$j, x=Y_data$freq)

cat("creating matrix X...")
#DONE: X(i,j) is overlaps of Region1(i) and Fragment with Tb(j)
system.time({
  X_i=overlap1[,1]
  X_j=sapply(overlap1[,2],FUN=function(x) { tag2tidx[fragdata$bar[x]] })
  cat("finshed X_j\n")
  X_data=data.frame(i=X_i,j=X_j)
  X_data=X_data[complete.cases(X_data),]
  X_data=count(X_data, vars=c("i","j"))
  X=sparseMatrix(i=X_data$i, j=X_data$j, x=X_data$freq)
  #for(i in seq_len(nrow(overlap1))){
  #  gidx=overlap1[i,1]
  #  fidx=overlap1[i,2]
  #  tidx=tag2tidx[[toString(fragdata$bar[fidx])]] #cann't find tag in tag2tidx
  #  if(!is.null(tidx)) { 
  #    X[gidx,tidx]=X[gidx,tidx]+1
  #  }
  #}
  #X <- Matrix(0, nrow = length(grid1), ncol = length(tag2tidx), sparse = TRUE)
})
print(dim(X))
print(nnzero(X))
cat("done\n")

cat("creating matrix Y...")
#DONE: Y(j,k) is indicator if overlap of Region2(k) and Fragment with Tb(j)
system.time({
  Y_j=overlap2[,1]
  Y_i=sapply(overlap2[,2],FUN=function(x) { tag2tidx[fragdata$bar[x]] })
  cat("finshed Y_j\n")
  Y_data=data.frame(i=Y_i,j=Y_j)
  Y_data=Y_data[complete.cases(Y_data),]
  Y_data=count(Y_data, vars=c("i","j"))
  Y_data$freq=(Y_data$freq>=1)
  Y=sparseMatrix(i=Y_data$i, j=Y_data$j, x=Y_data$freq)
#Y <- Matrix(0, nrow = length(tag2tidx), ncol = length(grid2), sparse = TRUE)
#for(i in seq_len(nrow(overlap2))){
#  gidx=overlap2[i,1]
#  fidx=overlap2[i,2]
#  tidx=tag2tidx[[toString(fragdata$bar[fidx])]] #cann't find tag in tag2tidx
#  if(!is.null(tidx)) { 
#    Y[tidx,gidx]=1
#  }
#}
})
print(dim(Y))
print(nnzero(Y))
cat("done\n")

#DONE create a sparse matrix M = X * t(Y)
cat("multiplying matrix X and Y...")
M <- X %*% Y
cat("done\n")

#compoute tcut from FDR control: affordable computations N_c = M^2 * (1-ppois(u0,t))
tcut = qpois(1-expectFalseDiscover/(nrow(M)**2),as.numeric(expectShareFrag))
cat("tcut:", tcut, "\n")

#Select region pair (i,j) that have Y(i,j)>t
sel=which(M > tcut, arr.ind=T)
res=as.data.frame(cbind(sel,M[sel]))
names(res)=c("grid1","grid2","m12")
res=res[order(-res$m12),]
cat("writing results...\n")
save.image(file=paste(outputFile,"RData",sep="."))
write.table(res, file=outputFile,quote=F,row.names=F,sep="\t")
cat("done\n")
cat("GRIDSCANDONE\n")
