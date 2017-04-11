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

for(p in c("RColorBrewer","optparse","ggplot2","MASS","rjson")) myrequire(p)

option_list2 <- list(
  make_option(c("-p", "--outPrefix"), default="none",
              help="outPrefix for separate plots separated by ';', none=mcvFiles [default %default]"),
  make_option(c("-v", "--overlapFactor"), type="integer", default=1,
              help="a factor of 2 means windows are half overlapping [default %default]"),
  make_option(c("-d", "--referenceGridsize"), type="integer", default=10000,
              help="gridSize used ad reference [default %default]"),
  make_option(c("-n", "--genomeSize"), type="numeric", default=3200000000,
              help="genomeSize for plot [default %default]"),
  make_option(c("-f", "--format"), default="tiff",
              help="output format, [default %default]"),
  make_option(c("-s", "--plotScale"), default=1,
              help="output scaling, [default %default]"),
  make_option(c("-l", "--label"), default="none",
              help="output label, [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)

parser <- OptionParser(usage = "%prog [options] mcvFile1;mcvFile2;... gridSize1;gridSize2;... cfgFile",
                       option_list=option_list2)
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)
if(length(cmd$args)!=3 & !cmd$options$debug){ print_help(parser); quit() }

mcvFiles=strsplit(cmd$args[1],";")[[1]]
gridSize=as.integer(strsplit(cmd$args[2],";")[[1]])
cfgFile=cmd$args[3]
outPrefix=strsplit(cmd$options$outPrefix,";")[[1]]
format=cmd$options$format
genomeSize=cmd$options$genomeSize
overlapFactor=cmd$options$overlapFactor
referenceGridsize=cmd$options$referenceGridsize
plotScale=cmd$options$plotScale
label=cmd$options$label
debug=cmd$options$debug
noQuiet=cmd$options$noQuiet
if(outPrefix=="none") outPrefix=mcvFiles
if(label=="none") label=""

if(debug){
  setwd("~/ws/10x6")
  mcvFiles = c("NA12878_WGS_phased_possorted_bam.bam.w100Ks50K.mcv","NA12878_WGS_phased_possorted_bam.bam.w100Ms50M.mcv",
               "NA12878_WGS_phased_possorted_bam.bam.w10Ks5K.mcv","NA12878_WGS_phased_possorted_bam.bam.w10Ms5M.mcv",
               "NA12878_WGS_phased_possorted_bam.bam.w1Ms500K.mcv")
  outPrefix = mcvFiles
  cfgFile = "NA12878_WGS_phased_possorted_bam.bam.gridScan.cfg.json"
  gridSize = as.integer(c("100000","100000000","10000","10000000","1000000"))
  referenceGridsize = 10000
  overlapFactor = 1
  genomeSize = 3200000000
  plotScale = 10
  format = "tiff"
  label = "NA12878_"
}
colScheme="Dark2"
cols <- brewer.pal(8, colScheme)[3:4]
colnames = paste(label,c("Theory","Empirical"),sep="")
names(cols) = colnames

mainTextSize=14
legendTextSize=16
mainLineSize=1.5
plotMargin = unit(c(0.28,0.28,0.28,0.28), "in")
cfgOptions = fromJSON(readLines(cfgFile))
averageFragLength = cfgOptions$mLf
#fragNumber = cfgOptions$aNf
fragNumber = cfgOptions$aNf
coveragePerBlock = round(fragNumber / (genomeSize/gridSize) )
names(coveragePerBlock) = as.character(gridSize)
tickFactor = overlapFactor / (gridSize/referenceGridsize)  # times tick to correct for x labels

data = lapply(mcvFiles, "read.table")
cov = lapply(seq_along(data), FUN=function(x){ data[[x]][,4][data[[x]][,4]>10] } )

plotDenCov = function(i, cov, outPrefix, coveragePerBlock, cols, tickFactor, plotScale, format){
  prefix=outPrefix[[i]];
  colnames=names(cols)
  Coverage = cov[[i]]; cpb=coveragePerBlock[[i]]; tf=tickFactor[[i]]; gs=names(coveragePerBlock)[[i]]
  xlim=round(c(0.25*cpb,1.75*cpb)/10)*10  #round to 10-based
  xticks=seq(xlim[1],xlim[2],by=50/tf) 
  xticklabs=paste(as.integer(seq(xlim[1]*tf,xlim[2]*tf,by=50)))
  ylim=c(0, 1.25*max(dpois(seq(1,1.5*cpb),lambda=cpb)))
  cat("theory: Poisson lamba=", cpb, "\n")
  cat("xlim:", xlim, "\n")
  cat("ylim:", ylim, "\n")
  p0 = qplot(Coverage, geom = 'blank') +   
    geom_histogram(aes(y = ..density..), alpha = 0.4, fill="blue", color="gray", binwidth=4/tf) +
    geom_line(aes(y = ..density.., colour=colnames[2]), stat="density", size=mainLineSize) + #Empricial
    geom_line(aes(x = seq(xlim[1], xlim[2]), 
                  y = dpois(seq(xlim[1], xlim[2]), lambda=cpb), 
                  colour=colnames[1]), #Theory
              size=mainLineSize) +
    theme_classic() +
    scale_colour_manual(name = '', values = cols ) +
    scale_x_continuous('Count', breaks=xticks, label=xticklabs, limit=xlim) +
    ylim(ylim) + ggtitle(paste("Molecule Count per", as.integer(gs)/1000,"KB")) + ylab("Density") + 
    theme(legend.position = c(0.85,0.85), text=element_text(size=legendTextSize)) + 
    theme(plot.margin=plotMargin) + 
    theme(plot.title = element_text(hjust=0.5, size=as.integer(mainTextSize), face="bold")) #+
  message("plotFile=",paste(prefix, format, sep="."))
  ggsave(filename=paste(prefix, format, sep="."),plot=p0,w=5.6, h=5.6,units="in",dpi=72*plotScale,
         type="Xlib", compression="lzw")
}

sapply(seq(length(mcvFiles)), "plotDenCov", cov, outPrefix, coveragePerBlock, cols=cols, 
       tickFactor=tickFactor, plotScale=plotScale, format=format)

quit()

#obsolete proof-of-idea codes
sapply(1, "plotDenCov", cov, outPrefix, coveragePerBlock, cols=cols, plotScale=plotScale)

testDenCov = function(i, cov, files, coveragePerBlock) {
  cv = cov[[i]]; prefix=files[[i]]; cpb=coveragePerBlock[[i]]
  my.mle<-fitdistr(cv, densfun="poisson")
  BIC(my.mle)
  cv_select = seq(round(0.5*median(cv)),round(1.5*median(cv))) #exclude outliers
  tcv = cv[ cv %in% cv_select ]
  tb = as.data.frame(table(tcv))
  tb$tcv = as.integer(as.character(tb$tcv))
  tb$Freq = as.integer(as.character(tb$Freq))
  tcv_max = max(tb$tcv)
  tcv_min = min(tb$tcv)
  for(j in seq(tcv_min, tcv_max)) if(! j %in% tb$tcv) tb=rbind(tb, data.frame(tcv=j,Freq=0))
  tb$Expect = dpois(tb$tcv, median(tcv)) # or mean of (tcv)
  tb$P = tb$Expect/(sum(tb$Expect))
  chisq.test(tb$Freq, p=tb$P, simulate.p.value=TRUE)
  cat("testing", i, "-th figure", paste(prefix, "pdf", sep="."), "\n")
  #ks.test(cv+rnorm(length(cv))/1000, "ppois", mean(cv))
  cat("theory: Poisson lamba=", cpb, "\n")
  cat("xlim:", xlim, "\n")
  cat("ylim:", ylim, "\n")
}


# plotDenCov = function(i, cov, files, coveragePerBlock){
#   cols <- brewer.pal(3, "Set1")
#   cv = cov[[i]]; prefix=files[[i]]; cpb=coveragePerBlock[[i]]
#   cat("plotting",paste(prefix, "png", sep="."),"\n")
#   png(paste(prefix, "png", sep="."), h=1000, w=1000, pointsize=20)
#   hp = hist(cv, breaks=100, ylim=c(0,0.01), xlim=c(0,2*cpb), prob=T, col=cols[1])            # prob=TRUE for probabilities not counts
#   print(hp)
#   poix = seq(min(hp$breaks),max(hp$breaks)) #hp$breaks
#   poiy = dpois(poix,lambda=cpb) #cpbp coveragePerBlock
#   #lines(density(cv), col=cols[2])             # add a density estimate with defaults
#   lines(density(cv, adjust=2), lty="dotted", ylim=hp$ylim*2, col=cols[2])   # add another "smoother" density
#   #points(poix, poiy, col=cols[3])
#   lines(poix, poiy, lty="dotted", col=cols[3])   # add another "smoother" density
#   # Add a legend using the nice sample labeles rather than the full filenames.
#   #legend("topright", legend=labels, col=cols, lty=1, lwd=4)
#   dev.off()
# }
# 
# lapply(seq_along(cov), "plotDenCov", cov, files, coveragePerBlock)

#ggplot2
# create some data to work with

x = rnorm(1000);

# overlay histogram, empirical density and normal density
p0 = qplot(x, geom = 'blank') +   
  geom_line(aes(y = ..density.., colour = 'Empirical'), stat = 'density') +  
  stat_function(fun = dnorm, aes(colour = 'Normal')) +                       
  geom_histogram(aes(y = ..density..), alpha = 0.4) +                        
  scale_colour_manual(name = 'Density', values = c('red', 'blue')) + 
  theme(legend.position = c(0.85, 0.85))

print(p0)
