#!/usr/bin/env Rscript
# Assumes you've already run coverageBed -hist, and grep'd '^all'. E.g. something like:
# find *.bam | parallel 'bedtools -abam {} -b capture.bed -hist | grep ^all > {}.all.txt'
# find *.bam | parallel 'bedtools -abam {} -b capture.bed -hist | grep ^all > {}.all.txt'

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

myrequire = function(pkg, repo="CRAN", lib=Sys.getenv("R_LIBS_USER"), ...){
  cat("requiring package", pkg, "\n")
  tryCatch(suppressMessages(library(pkg,character.only=T)), error=function(e) {
    print(e)
    if(repo!="CRAN"){
      source("http://bioconductor.org/biocLite.R")
      biocLite(pkg,...)
    } else {
      install.packages(pkg,repo="http://cran.us.r-project.org",lib=lib,...)
    }
  })
  tryCatch(suppressMessages(library(pkg,character.only=T)), error=function(e) {
    print(e)
    stop(pkg," was not installed and cannot install on the fly!\n")
  })
}

for(p in c("RColorBrewer","optparse")) myrequire(p)
#for(p in c("GenomicRanges", "ggbio", "GenomeInfoDb")) myrequire(p, repo="bioC")

option_list2 <- list(
  make_option(c("-o", "--outPrefix1"), default="none",
              help="outPrefix for a merged plot without normalization, none=not plotting [default %default]"),
  make_option(c("-p", "--outPrefix2"), default="none",
              help="outPrefix for separate plots without normalization, none=not plotting [default %default]"),
  make_option(c("-r", "--outPrefix3"), default="none",
              help="outPrefix for separate plots with normalization, none=not plotting [default %default]"),
  make_option(c("-l", "--labels"), default="none",
              help="none=$histSet1File1,$histSet1File2;$histSet2File1,$histSet2File2;... [default %default]"),
  make_option(c("-i", "--insertSize"), type="integer", default=340,
              help="mean insert size [default %default]"),
  make_option(c("-d", "--readPairSize"), type="integer", default=280,
              help="read pair size, rl1 + rl2 - barcode [default %default]"),
  make_option(c("-m", "--maxCoverage"), type="integer", default=400,
              help="maximum cvoerage [default %default]"),
  make_option(c("-f", "--format"), default="tiff",
              help="output format, [default %default]"),
  make_option(c("-s", "--plotScale"), default=10,
              help="output scaling, [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)

parser <- OptionParser(usage = "%prog [options] $histSet1File1,$histSet1File2;$histSet2File1,$histSet2File2;...",
                       option_list=option_list2)
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)
if(length(cmd$args)!=1 & !cmd$options$debug){ print_help(parser); quit() }
histFiles=strsplit(cmd$args[1],",")[[1]]
labels=histFiles
outPrefix1=cmd$options$outPrefix1                                #plotMerge
outPrefix2=strsplit(cmd$options$outPrefix2,";")                  #plotSeparate
outPrefix3=strsplit(cmd$options$outPrefix3,";")                  #plotSepNorm
format=cmd$options$format
plotScale=cmd$options$plotScale
debug=cmd$options$debug
noQuiet=cmd$options$noQuiet
readPairSize=cmd$options$readPairSize
insertSize=cmd$options$insertSize
maxCoverage=cmd$options$maxCoverage
colScheme="Dark2"
if(cmd$options$labels!="none") labels=strsplit(unlist(strsplit(cmd$options$labels,";")),",")
if(outPrefix1=="none") plotMerge=FALSE
if(outPrefix2=="none") plotSeparate=FALSE
if(outPrefix3=="none") plotSepNorm=FALSE

# Get a list of the bedtools output histFiles you'd like to read in
if(debug){
  setwd("~/ws/10x6")
  histFiles = list(c("NA12878_WGS_phased_possorted_bam.bam.hist.all.txt","NA12878_WGS_phased_possorted_bam.bam.multiFragBar.bed.a0.hist.all.txt"),c("P00526_2386_nt.10X_v2.bam.hist.all.txt","P00526_2386_nt.10X_v2.bam.multiFragBar.bed.a0.hist.all.txt"),c("P00526_2721_mt.st.10X_v2.bam.hist.all.txt","P00526_2721_mt.st.10X_v2.bam.multiFragBar.bed.a0.hist.all.txt"),c("P00526_2725_mt.10X_v2.bam.hist.all.txt","P00526_2725_mt.10X_v2.bam.multiFragBar.bed.a0.hist.all.txt"))
  labels = list(c("NA12878_PEfrag","NA12878_molecule"),c("Norm2386_PEfrag","Norm2386_molecule"),c("MetR2721_PEfrag","MetR2721_molecule"),c("MetL2725_PEfrag","MetL2725_molecule"))
  outPrefix1 = "zoomx.cumcov.plot"
  plotMerge=TRUE
  outPrefix2 = c("NA12878_WGS_phased_possorted_bam.cumcov.plot","P00526_2386_nt.10X_v2.bam.cumcov.plot","P00526_2721_mt.st.10X_v2.bam.cumcov.plot","P00526_2725_mt.10X_v2.bam.cumcov.plot")
  plotSeparate=TRUE
  outPrefix3 = c("NA12878_WGS_phased_possorted_bam.cumncov.plot","P00526_2386_nt.10X_v2.bam.cumncov.plot","P00526_2721_mt.st.10X_v2.bam.cumncov.plot","P00526_2725_mt.10X_v2.bam.cumncov.plot")
  plotSepNorm=TRUE
  plotScale=10
  readPairSize=280
  insertSize=340
  maxCoverage=400
  format="tiff"
}
insertFactor=insertSize/readPairSize

# define the plotCumCov function
# histFiles=mergeHistFiles
# labels=mergeLabels
plotCumCov = function(histFiles, labels, pf, insertFactor, format, addNorm=T, plotScale=10) { 
  
  # Create lists to hold coverage and cumulative coverage for each alignment,
  # and read the data into these lists.
  cov <- list()
  cov_cumul <- list()
  for (i in 1:length(histFiles)) {
    cov[[i]] <- read.table(histFiles[i])
    names(cov[[i]]) = c("region","coverage","bases","totalBases","density")
    cov[[i]][["cumulative"]] = 1-cumsum(cov[[i]][,5])
    if(i%%2==1) {
      x = cov[[i]][["coverage"]]*insertFactor
      y = cov[[i]][["cumulative"]]
      sp <- smooth.spline(y[1:(2*maxCoverage+1)]~x[1:(2*maxCoverage+1)], all.knots=FALSE)
      xmax = floor(length(cov[[i]][["cumulative"]])*insertFactor)+1 #fit range from 1:xmax
      tmpCumCov = pmax(predict(sp, 0:xmax)$y,0)
      tmpCov = list()
      tmpCov[["region"]] = rep("all",2*maxCoverage+1)
      tmpCov[["coverage"]] = 0:(2*maxCoverage)
      tmpCov[["bases"]] = NA
      tmpCov[["totalBases"]] = NA
      tmpCov[["density"]] = -diff(c(1, tmpCumCov))[1:(2*maxCoverage+1)]
      tmpCov[["cumulative"]] = tmpCumCov[1:(2*maxCoverage+1)]
      cov[[i]] = as.data.frame(tmpCov)
    }
    cov_cumul[[i]] = cov[[i]][["cumulative"]]
  }
  #max(which(cov[[1]]$cumulative>0.5))
  #max(which(tmpCov$cumulative>0.5))
  #max(which(cov[[1]]$cumulative>0.5))
  #max(which(tmpCov$cumulative>0.8))
  #max(which(cov[[1]]$cumulative>0.8))
  #max(which(cov[[2]]$cumulative>0.8))
  #sqrt(27)=5.19 ~=5
  #sqrt(176)=13.3 ~=13
  #sqrt(33)=5.74 ~=6
  #-diff(cov[[1]][["cumulative"]][c(27-5+1,27+5+1)])
  #-diff(tmpCov[["cumulative"]][c(33-6+1,33+6+1)])
  #-diff(cov[[2]][["cumulative"]][c(176-13+1,176+13+1)])
  
  # Pick some colors
  # Ugly:
  # cols <- 1:length(cov)
  # Prettier:
  # ?colorRampPalette
  # display.brewer.all()
  cols <- brewer.pal(length(cov), colScheme)
  
  # Save the graph to a file
  if(format=="tiff") {
    tiff(paste(pf, "tiff", sep="."), w=ceiling(5.6*plotScale*72/2)*2, h=ceiling(5.6*plotScale*72/2)*2, units="px", 
         res=72*plotScale, type="Xlib", compression="lzw")
    message("tiffFile=",paste(pf, "tiff", sep="."))
  } else {
    pdf(paste(pf, "pdf", sep="."), w=5.6, h=5.6, units="in")
    message("pdfFile=",paste(pf, "pdf", sep="."))
  }
  
  # Create plot area, but do not plot anything. Add gridlines and axis labels.
  plot(cov[[1]][2:(maxCoverage+1), 2], cov_cumul[[1]][1:maxCoverage], type='n', xlab="Coverage", 
       ylab="Fraction of genome at coverage", ylim=c(0,1.0), main="Cumulative Genome Coverage", cex=plotScale)
  
  abline(v = 50, col = "gray")
  abline(v = 100, col = "gray")
  abline(h = 0.50, col = "gray")
  abline(h = 0.80, col = "gray")
  #axis(1, at=c(20,50,80), labels=c(20,50,80))
  axis(1, at=c(50), labels=c(50))
  axis(2, at=c(0.80), labels=c(0.80))
  axis(2, at=c(0.50), labels=c(0.50))
  
  # Actually plot the data for each of the alignments (stored in the lists).
  tmp_cov=cov; tmp_cov_cumul=cov_cumul; tmp_labels=labels; tmp_lty=rep(1,length(cov))
  if(addNorm) { #add normaliation step to cov[[ai]], cov_cumul[[ai]], labels[ai] based on cov[[bi]]
    ai=1; bi=2; 
    # align ai to bi
    bmid = which(cov_cumul[[bi]]<0.5)[1]   #29
    amid = which(cov_cumul[[ai]]<0.5)[1]    #178
    al_step = amid/bmid                    #alignment step
    alab = paste(labels[ai],"(normalized)",sep="")
    tmp_cov[[ai]]=cov[[ai]]
    tmp_cov_cumul[[ai]]=rep(NA,maxCoverage)
    for(s in 1:maxCoverage) {
      sj = which( (1:length(cov_cumul[[ai]]))/al_step > s )[1]
      if( is.na(sj) ) sj = length(cov_cumul[[ai]])
      tmp_cov_cumul[[ai]][s] = cov_cumul[[ai]][sj]
    }
    #x = 1:401; y=tmp_cov_cumul[[ai]]; lo <- loess(y~x); tmp_cov_cumul[[ai]]=predict(lo)
    x = 1:maxCoverage; y=tmp_cov_cumul[[ai]]; 
    sp <- smooth.spline(y~x, all.knots=FALSE); tmp_cov_cumul[[ai]]=predict(sp)$y
    tmp_labels[ai] = alab
    tmp_lty[ai] = 2
  } 
  for (i in 1:length(tmp_cov)) 
    points(tmp_cov[[i]][2:(maxCoverage+1), 2], tmp_cov_cumul[[i]][1:maxCoverage], type='l', lwd=3, lty=tmp_lty[i], col=cols[i])
  
  # Add a legend using the nice sample labeles rather than the full filenames.
  legend("topright", legend=tmp_labels, col=cols, lty=tmp_lty, lwd=3, bty="n") #,box.lty=NA,box.col="white",bg="white")
  
  dev.off()
}

# plotMerge
if(plotMerge) {
  mergeHistFiles = unlist(histFiles)
  mergeLabels = unlist(labels)
  plotCumCov( mergeHistFiles, mergeLabels, outPrefix1, insertFactor=insertFactor, format=format, addNorm=F, plotScale=plotScale )
}

if(plotSeparate) {
  lapply( seq(length(histFiles)), function(x) { plotCumCov( histFiles[[x]], labels[[x]], outPrefix2[x], insertFactor=insertFactor, format=format, addNorm=F, plotScale=plotScale ) })
}

if(plotSepNorm) {
  lapply( seq(length(histFiles)), function(x) { plotCumCov( histFiles[[x]], labels[[x]], outPrefix3[x], insertFactor=insertFactor, format=format, addNorm=T, plotScale=plotScale ) })
}

quit()

# obsolete codes: proof of idea
pf = "all.cumcov"
  giabpf = "giab.cumcov"
  giabnormpf = "giab.norm.cumcov"
  giabhistFiles = histFiles[grepl("NA12878",histFiles)]
  giablabels = labels[grepl("NA12878",labels)]
  normpf = "norm.cumcov"
  normnormpf = "norm.norm.cumcov"
  normhistFiles = histFiles[grepl("2386",histFiles)]
  normlabels = labels[grepl("2386",labels)]
  metRpf = "metR.cumcov"
  metRnormpf = "metR.norm.cumcov"
  metRhistFiles = histFiles[grepl("2721",histFiles)]
  metRlabels = labels[grepl("2721",labels)]
  metLpf = "metL.cumcov"
  metLnormpf = "metL.norm.cumcov"
  metLhistFiles = histFiles[grepl("2725",histFiles)]
  metLlabels = labels[grepl("2725",labels)]

plotCumCov( giabhistFiles, giablabels, giabpf, addNorm=F )

plotCumCov( normhistFiles, normlabels, normpf, addNorm=F )

plotCumCov( metRhistFiles, metRlabels, metRpf, addNorm=F )

plotCumCov( metLhistFiles, metLlabels, metLpf, addNorm=F )

plotCumCov( giabhistFiles, giablabels, giabnormpf, addNorm=T )

plotCumCov( normhistFiles, normlabels, normnormpf, addNorm=T )

plotCumCov( metRhistFiles, metRlabels, metRnormpf, addNorm=T )

plotCumCov( metLhistFiles, metLlabels, metLnormpf, addNorm=T )


