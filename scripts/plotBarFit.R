#!/usr/bin/env Rscript

####Externalities####
ptm=proc.time() #time keeper

set_colClass <- function(d, colClasses) {
  colClasses <- rep(colClasses, len=length(d))
  d[] <- lapply(seq_along(d), function(i) switch(colClasses[i],
                                                 numeric=as.numeric(d[[i]]),
                                                 character=as.character(d[[i]]),
                                                 Date=as.Date(d[[i]], origin='1970-01-01'),
                                                 POSIXct=as.POSIXct(d[[i]], origin='1970-01-01'),
                                                 factor=as.factor(d[[i]]),
                                                 as(d[[i]], colClasses[i]) ))
  return(d)
}

myrequire = function(pkg, repo="CRAN", version=NULL, lib=Sys.getenv("R_LIBS_USER"), ...){
  if(!is.null(version)) { #always install
    install_version(pkg, version = version, repos=repo)
  }
  if(repo!="CRAN" & repo!="Bioc") { #always install
    install_git(repo)
  }
  tryCatch(suppressMessages(library(pkg,character.only=T)), error=function(e) {
    print(e)
    if(repo=="Bioc"){
      source("http://bioconductor.org/biocLite.R")
      biocLite(pkg,...)
    } else if (repo=="CRAN") {
      install.packages(pkg,repo="http://cran.us.r-project.org",lib=lib,...)
    }
  })
  tryCatch(suppressMessages(library(pkg,character.only=T)), error=function(e) {
    print(e)
    stop(pkg," was not installed and cannot install on the fly!\n")
  })
}

for(p in c("optparse","devtools","fitdistrplus","ggplot2")) myrequire(p)
for(p in c("GenomicRanges","rtracklayer")) myrequire(p, repo="Bioc")


option_list2 <- list(
  make_option(c("-p", "--outPrefix"), default="none",
              help="outPrefix for separate plots separated by ';', none=mcvFiles [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)

parser <- OptionParser(usage = "%prog [options] moleFile1;moleFile2;... gridFile1;gridFile2;...",
                       option_list=option_list2)
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)
if(length(cmd$args)!=2 & !cmd$options$debug){ print_help(parser); quit() }
cat("-Info: invoking command:",commandArgs(),"\n")

moleFiles=strsplit(cmd$args[1],";")[[1]]
gridFiles=strsplit(cmd$args[2],";")[[1]]
outPrefix=strsplit(cmd$options$outPrefix,";")[[1]]
debug=cmd$options$debug
noQuiet=cmd$options$noQuiet
if(outPrefix=="none") outPrefix=moleFiles

if(debug){
  setwd("~/ws/10x6")
  moleFiles = c("NA12878_from10X.wg12.bam.multiFragBar.bed.a0")
  gridFiles = c("~/ws/hg/hg19/human_g1k_v37.ucsc.fasta.bed.w10Ks10K")
  outPrefix = moleFiles
}

moleTableHeader=c("chr1","st1","ed1","barID","molSize","molLen","barMolNum","molID","molHP","molPS","molBX")
refGridHeader=c("chr1","st1","ed1")

readRefGrid = function(i, gridFiles) {
  tmp=read.table(gridFiles[i],sep="\t",stringsAsFactors = FALSE)
  colnames(tmp) = refGridHeader
  grRefGrid = GRanges(seqnames=tmp$chr1
                      , ranges=IRanges(start=tmp$st1, end=tmp$ed1))
  grRefGrid
} #grRefGrid = readRefGrid(1, gridFiles)

readMoleTable = function(i, moleFiles) {
  #tmp=import.bed(moleFiles[i]) #to use import.bed, has to follow UCSC bed, bed15, bedgraph or bedpe formats
  tmp=read.table(moleFiles[i],sep="\t",stringsAsFactors = FALSE)
  colnames(tmp) = moleTableHeader
  grMoleTable = GRanges(seqnames=tmp$chr1
                      , ranges=IRanges(start=tmp$st1, end=tmp$st1+1)
                      , bx=tmp$molBX)
  grMoleTable
} #grMoleTable=readMoleTable(1, moleFiles)

runBarFit = function(i, moleFiles, gridFiles) {
  message("runBarFit",i,"-th moleFile")
  grRefGrid = readRefGrid(i, gridFiles)
  grMoleTable=readMoleTable(i, moleFiles)
  ovlaps = findOverlaps(grMoleTable,grRefGrid)
  ovlapGridBX = cbind(subjectHits(ovlaps), grMoleTable[queryHits(ovlaps),]$bx)
  names(ovlapGridBX) = c("subject","bx")
  uniqueBX = unique(grMoleTable$bx) 
  stats = lapply(uniqueBX, testBarFit, ovlapGridBX=ovlapGridBX, grRefGrid=grRefGrid)
  stats
}#allTests = runBarFit(1, moleFiles, gridFiles) 

testBarFit = function(bx, ovlapGridBX, grRefGrid) {
  message(bx)
  tmp = matrix(ovlapGridBX[ovlapGridBX[,2] == bx,],ncol=2)
  #print(tmp)
  res=list(pvalue=NA, lambda=NA)
  if(nrow(tmp)>0) {
    bx_hist = as.data.frame(table(tmp[,1]))
    bx_data = rep(0,length(grRefGrid))
    bx_data[as.numeric(as.character(bx_hist[,1]))]=bx_hist[,2]
    #test if Poisson Distributed
    pois.fit = fitdist(bx_data,"pois",method="mme")
    pois.stat = gofstat(pois.fit, chisqbreaks=c(0,sort(unique(bx_hist[,2])))) #goodness-of-fit test
    res=list(pvalue=pois.stat$chisqpvalue,lambda=pois.fit$estimate)
  }
  res
  #chisq.fit = chisq.test(x=c((length(grRefGrid)-nrow(bx_hist)), nrow(bx_hist)), p=dpois(0:1,lambda=bx_mean)) 
}#bx=uniqueBX[1]; testBarFit(bx, ovlapGridBX)

#need to do a density plot

allTests = lapply(seq_along(moleFiles),runBarFit,moleFiles,gridFiles)
for(i in seq_along(moleFiles)) {
  pvalues = sapply(allTests[[i]],"[[",1)
  lambdas = sapply(allTests[[i]],"[[",2)
  save(pvalues, lambdas, file=paste(outPrefix[i],".barfit.RData",sep=""))
  plt = ggplot()
  plt = plt+geom_histogram(aes(pvalues),breaks=seq(0,1,0.01))+scale_x_continuous(limits=c(0,1.01))
  ggsave(file=paste(outPrefix[i],".barfit.tiff",sep=""))
}

quit()

#m=c(40,34,32,30,29,27,26,50,34,27,27,26,26)
#r=c(18,5,2,0,1,0,4,26,7,0,2,11,2)
#t.test(m,r,paired=T,alternative="greater")
