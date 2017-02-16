#!/usr/bin/env Rscript
# Assumes you've already run coverageBed -hist, and grep'd '^all'. E.g. something like:
# find *.bam | parallel 'bedtools -abam {} -b capture.bed -hist | grep ^all > {}.all.txt'

# Get a list of the bedtools output files you'd like to read in
setwd("~/ws/10x6")
files = list.files(pattern="all.txt$")
labels = gsub("_WGS_phased_possorted_bam.bam.hist.all.txt","_read",files)
labels = gsub("_WGS_phased_possorted_bam.bam.multiFragBar.bed.a0.hist.all.txt","_molecule",labels)
labels = gsub(".10X_v2.bam.hist.all.txt","_read",labels)
labels = gsub(".10X_v2.bam.multiFragBar.bed.a0.hist.all.txt","_molecule",labels)
pf = "all.cumcov"
giabpf = "giab.cumcov"
giabnormpf = "giab.norm.cumcov"
giabfiles = files[grepl("NA12878",files)]
giablabels = labels[grepl("NA12878",labels)]
normpf = "norm.cumcov"
normnormpf = "norm.norm.cumcov"
normfiles = files[grepl("2386",files)]
normlabels = labels[grepl("2386",labels)]
metRpf = "metR.cumcov"
metRnormpf = "metR.norm.cumcov"
metRfiles = files[grepl("2721",files)]
metRlabels = labels[grepl("2721",labels)]
metLpf = "metL.cumcov"
metLnormpf = "metL.norm.cumcov"
metLfiles = files[grepl("2725",files)]
metLlabels = labels[grepl("2725",labels)]

plotCumCov = function(files, labels, pf, addNorm=T) { 

  # Create lists to hold coverage and cumulative coverage for each alignment,
  # and read the data into these lists.
  cov <- list()
  cov_cumul <- list()
  for (i in 1:length(files)) {
      cov[[i]] <- read.table(files[i])
      cov_cumul[[i]] <- 1-cumsum(cov[[i]][,5])
  }

  # Pick some colors
  # Ugly:
  # cols <- 1:length(cov)
  # Prettier:
  # ?colorRampPalette
  # display.brewer.all()
  library(RColorBrewer)
  cols <- brewer.pal(length(cov)+1, "Set1")

  # Save the graph to a file
  png(paste(pf, "png", sep="."), h=750, w=750, pointsize=20)

  # Create plot area, but do not plot anything. Add gridlines and axis labels.
  plot(cov[[1]][2:401, 2], cov_cumul[[1]][1:400], type='n', xlab="Coverage", 
       ylab="Fraction of genome \u2265 coverage", ylim=c(0,1.0), main="Genome Coverage (Cumulative)")
  abline(v = 20, col = "gray60")
  abline(v = 50, col = "gray60")
  abline(v = 80, col = "gray60")
  abline(v = 100, col = "gray60")
  abline(h = 0.50, col = "gray60")
  abline(h = 0.90, col = "gray60")
  axis(1, at=c(20,50,80), labels=c(20,50,80))
  axis(2, at=c(0.90), labels=c(0.90))
  axis(2, at=c(0.50), labels=c(0.50))

  # Actually plot the data for each of the alignments (stored in the lists).
  tmp_cov=cov; tmp_cov_cumul=cov_cumul; tmp_labels=labels; tmp_lty=rep(1,length(cov))
  if(addNorm) { #add normaliation step to cov[[ai]], cov_cumul[[ai]], labels[ai] based on cov[[bi]]
    ai=1; bi=2; 
    # align ai to bi
    bmid = which(cov_cumul[[bi]]<0.5)[1]   #29
    amid = which(cov_cumul[[ai]]<0.5)[1]    #178
    al_step = amid/bmid                    #alignment step
    alab = paste(labels[ai],"normalized",sep="_")
    tmp_cov[[ai]]=cov[[ai]]
    tmp_cov_cumul[[ai]]=rep(NA,400)
    for(s in 1:400) {
      sj = which( (1:length(cov_cumul[[ai]]))/al_step > s )[1]
      if( is.na(sj) ) sj = length(cov_cumul[[ai]])
      tmp_cov_cumul[[ai]][s] = cov_cumul[[ai]][sj]
    }
    #x = 1:401; y=tmp_cov_cumul[[ai]]; lo <- loess(y~x); tmp_cov_cumul[[ai]]=predict(lo)
    x = 1:400; y=tmp_cov_cumul[[ai]]; sp <- smooth.spline(y~x, all.knots=FALSE); tmp_cov_cumul[[ai]]=predict(sp)$y
    tmp_labels[ai] = alab
    tmp_lty[ai] = 2
  } 
  for (i in 1:length(tmp_cov)) points(tmp_cov[[i]][2:401, 2], tmp_cov_cumul[[i]][1:400], type='l', lwd=3, lty=tmp_lty[i], col=cols[i])

  # Add a legend using the nice sample labeles rather than the full filenames.
  legend("topright", legend=tmp_labels, col=cols, lty=tmp_lty, lwd=3)
  
  dev.off()

}

plotCumCov( files, labels, pf, addNorm=F )

plotCumCov( giabfiles, giablabels, giabpf, addNorm=F )

plotCumCov( normfiles, normlabels, normpf, addNorm=F )

plotCumCov( metRfiles, metRlabels, metRpf, addNorm=F )

plotCumCov( metLfiles, metLlabels, metLpf, addNorm=F )

plotCumCov( giabfiles, giablabels, giabnormpf, addNorm=T )

plotCumCov( normfiles, normlabels, normnormpf, addNorm=T )

plotCumCov( metRfiles, metRlabels, metRnormpf, addNorm=T )

plotCumCov( metLfiles, metLlabels, metLnormpf, addNorm=T )


