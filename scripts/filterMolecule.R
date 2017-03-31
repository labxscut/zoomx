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

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

findAneuploid = function(bv, y) {  #find windows in boolean vector bv with at leat y consecutive true  
  pbv=c(0,bv,0)
  edges=diff(pbv)
  rising=which(edges==1)
  falling=which(edges==-1)
  spanWidth = falling - rising
  wideEnough = spanWidth >= y
  startPos = rising[wideEnough]
  endPos = falling[wideEnough]-1
  as.vector(unlist(sapply(seq_along(startPos),FUN=function(x){ startPos[x]:endPos[x] })))
}

for(p in c("optparse", "IRanges")) myrequire(p)
for(p in c("GenomicRanges")) myrequire(p, repo="bioC")

####Options####
option_list2 <- list(
  make_option(c("-o", "--outPrefix"), type="integer", default="none",
              help="none=$molFile.clean [default %default]"),
  make_option(c("-x", "--xMadExclude"), type="integer", default=1,
              help="if coverage is x mad away from median, exclude regional molecules [default %default]"),
  make_option(c("-y", "--yMedianInclude"), type="integer", default=6,
              help="if coverage is y*median+z*mad or less from median and span more than w windows, include aneuploidity regions [default %default]"),
  make_option(c("-z", "--zMadInclude"), type="integer", default=6,
              help="if coverage is y*median+z*mad or less from median and span more than w windows, include aneuploidity regions [default %default]"),
  make_option(c("-w", "--windowsInclude"), type="integer", default=35,
              help="if coverage is y*median+z*mad or less from median and span more than w windows, include aneuploidity regions [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)

parser <- OptionParser(usage = "%prog [options] covFile moleFile", option_list=option_list2)
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)

outPrefix=cmd$options$outPrefix
xMadExclude=cmd$options$xMadExclude
yMedianInclude=cmd$options$yMedianInclude
zMadInclude=cmd$options$zMadInclude
windowsInclude=cmd$options$windowsInclude
dubug=cmd$options$debug
noQuiet=cmd$options$noQuiet

if(length(cmd$args)!=2 & !debug){ print_help(parser); quit() }
covFile=cmd$args[1]
moleFile=cmd$args[2]

####Debug####
if(debug){ 
  #options(warn=2)
  outPrefix = "none"
  xMadExclude = 7     
  yMedianInclude = 6
  zMadInclude = 6
  windowsInclude = 35
  covFile="P00526_2721_mt.st.10X_v2.bam.cov100"
  moleFile="P00526_2721_mt.st.10X_v2.bam.multiFragBar.bed.a0"
} else { 
  options(warn=-1)
}

if(outPrefix=="none") {
  cleanFile=paste(moleFile,"clean",sep=".")
} else {
  cleanFile=paste(outPrefix,"clean",sep=".")
}

cov=read.table(covFile)  #takes a while
mol=read.table(moleFile)  #takes a while
names(cov)=c("CHR","START","END","COV")
stats=list()
stats[["mean0"]]=mean(cov$COV)
stats[["mean1"]]=mean(cov$COV[covCOV!=0])
stats[["sd0"]]= sd(cov$COV)
stats[["sd1"]]= sd(cov$COV[ cov$COV!=0])
stats[["mode1"]]= Mode(cov$COV[cov$COV!=0])
stats[["mode0"]]= Mode(cov$COV)
stats[["mad0"]]= mad(cov$COV[cov$COV!=0])
stats[["mad1"]]= mad(cov$COV[cov$COV!=0])
stats[["median1"]]= median(cov$COV[ cov$COV!=0])
stats[["median0"]]= median(cov$COV)
statsTab=as.data.frame(stats)
excludeThresh=statsTab$median1+statsTab$mad1*xMadExclude
includeThresh=statsTab$median1*yMedianInclude+statsTab$mad1*zMadInclude

cov[["EX"]]=(cov$COV>excludeThresh)
cov[["IN"]] = FALSE
cov[["IN"]][findAneuploid(cov$COV<includeThresh,windowsInclude)] = TRUE
cov[["EX2"]] = cov[["EX"]] & (!cov[["IN"]])
cat("f=",moleFile,"windowsInclude=",windowsInclude,"allExcluded=",sum(cov[["EX2"]]),"\n")

regExclude = GRanges(seqnames = cov[["CHR"]][cov[["EX2"]]],
        ranges = IRanges(start = cov[["START"]][cov[["EX2"]]], end=cov[["START"]][cov[["EX2"]]]))
posPad=100
molStart = GRanges(seqnames = mol[,1],
                   ranges = IRanges(start = mol[,2]-posPad, end = mol[,2]+posPad))
molEnd = GRanges(seqnames = mol[,1],
                   ranges = IRanges(start = mol[,2]-posPad, end = mol[,2]+posPad))

overlapStart=!is.na(as.matrix(GenomicRanges::findOverlaps(molStart, regExclude, ignore.strand=TRUE, select="first"))[,1])  #
overlapEnd=!is.na(as.matrix(GenomicRanges::findOverlaps(molStart, regExclude, ignore.strand=TRUE, select="first"))[,1])    #
overlapEither = (overlapEnd | overlapStart)
cat("overlapEither=",sum(overlapEither),"totalMolecules=",length(overlapEither),"\n")

write.table(file=cleanFile,mol[!overlapEither,],sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
cat("covFilterDone in",proc.time()[3]-ptm[3],"seconds! \n")







quit()

#legacy proof-of-concept scripts
moleFiles=c()
files=c("MetL2725.wg12.bam.cov100", "MetR2721.wg12.bam.cov100", "Norm2386.wg12.bam.cov100")
files=c("P00526_2386_nt.10X_v2.bam.cov100","P00526_2721_mt.st.10X_v2.bam.cov100","P00526_2725_mt.10X_v2.bam.cov100")
moleFiles=c("P00526_2386_nt.10X_v2.bam.multiFragBar.bed.a0","P00526_2721_mt.st.10X_v2.bam.multiFragBar.bed.a0","P00526_2725_mt.10X_v2.bam.multiFragBar.bed.a0")
cov=lapply(files, read.table)
names(cov)=files
names(cov[[1]])=c("CHR","START","END","COV")
names(cov[[2]])=c("CHR","START","END","COV")
names(cov[[3]])=c("CHR","START","END","COV")
stats=list()
stats[["mean0"]]=sapply(cov, FUN=function(x) { mean(x$COV) })
stats[["mean1"]]=sapply(cov, FUN=function(x) { mean(x$COV[x$COV!=0]) })
stats[["sd0"]]=sapply(cov, FUN=function(x) { sd(x$COV) })
stats[["sd1"]]=sapply(cov, FUN=function(x) { sd(x$COV[x$COV!=0]) })
stats[["mode1"]]=sapply(cov, FUN=function(x) { Mode(x$COV[x$COV!=0]) })
stats[["mode0"]]=sapply(cov, FUN=function(x) { Mode(x$COV) })
stats[["mad0"]]=sapply(cov, FUN=function(x) { mad(x$COV[x$COV!=0]) })
stats[["mad1"]]=sapply(cov, FUN=function(x) { mad(x$COV[x$COV!=0]) })
stats[["median1"]]=sapply(cov, FUN=function(x) { median(x$COV[x$COV!=0]) })
stats[["median0"]]=sapply(cov, FUN=function(x) { median(x$COV) })
statsTab=as.data.frame(stats)
thresh67=statsTab$median1*6+statsTab$mad1*7
thresh66=statsTab$median1*6+statsTab$mad1*6
thresh65=statsTab$median1*6+statsTab$mad1*5
cut66=sapply(seq(3), FUN=function(x){ sum(cov[[x]]$COV>thresh66[x]) })
thresh15=statsTab$median1+statsTab$mad1*5
thresh14=statsTab$median1+statsTab$mad1*4
thresh16=statsTab$median1+statsTab$mad1*6
cut15=sapply(seq(3), FUN=function(x){ sum(cov[[x]]$COV>thresh15[x]) })
for(i in seq(3)) { cov[[i]][["EX"]]=(cov[[i]]$COV>thresh15[i]) }
  
v2=c(1,2,3,2,0,3,1,2,0)
ex2=setdiff(seq_along(v2),findAneuploid(v2<3,2)) #
findAneuploid(cov[[1]]$EX & cov[[1]]$COV<thresh66[1],10)

length(findAneuploid(cov[[1]]$EX & cov[[1]]$COV<thresh66[1],10))
length(findAneuploid(cov[[2]]$EX & cov[[2]]$COV<thresh66[2],10))
length(findAneuploid(cov[[3]]$EX & cov[[3]]$COV<thresh66[3],10))

length(findAneuploid(cov[[1]]$EX & cov[[1]]$COV<thresh66[1],20))
length(findAneuploid(cov[[2]]$EX & cov[[2]]$COV<thresh66[2],20))
length(findAneuploid(cov[[3]]$EX & cov[[3]]$COV<thresh66[3],20))

length(findAneuploid(cov[[1]]$EX & cov[[1]]$COV<thresh66[1],25))
length(findAneuploid(cov[[2]]$EX & cov[[2]]$COV<thresh66[2],25))
length(findAneuploid(cov[[3]]$EX & cov[[3]]$COV<thresh66[3],25))


for(t in c(30,35,40)) 
for(i in seq(3)) 
  {
    cov[[i]][["EX"]]=(cov[[i]]$COV>thresh15[i])
    cov[[i]][["IN"]] = FALSE
    cov[[i]][["IN"]][findAneuploid(cov[[i]]$COV<thresh66[i],t)] = TRUE
    cov[[i]][["EX2"]] = cov[[i]][["EX"]] & (!cov[[i]][["IN"]])
    cat("f=",files[i],"t=",t,"ex2=",sum(cov[[i]][["EX2"]]),"ex2chr10=",sum(cov[[i]][["EX2"]]&cov[[i]][["CHR"]]=="10"),"\n")
}

for(t in c(30,35,40)) 
  for(i in seq(3)) 
  {
    cov[[i]][["EX"]]=(cov[[i]]$COV>thresh14[i])
    cov[[i]][["IN"]] = FALSE
    cov[[i]][["IN"]][findAneuploid(cov[[i]]$COV<thresh66[i],t)] = TRUE
    cov[[i]][["EX2"]] = cov[[i]][["EX"]] & (!cov[[i]][["IN"]])
    cat("f=",files[i],"t=",t,"ex2=",sum(cov[[i]][["EX2"]]),"ex2chr10=",sum(cov[[i]][["EX2"]]&cov[[i]][["CHR"]]=="10"),"\n")
  }

for(t in c(20,25,30,35,40)) 
  for(i in seq(3)) 
  {
    cov[[i]][["EX"]]=(cov[[i]]$COV>thresh16[i])
    cov[[i]][["IN"]] = FALSE
    cov[[i]][["IN"]][findAneuploid(cov[[i]]$COV<thresh66[i],t)] = TRUE
    cov[[i]][["EX2"]] = cov[[i]][["EX"]] & (!cov[[i]][["IN"]])
    cat("f=",files[i],"t=",t,"ex2=",sum(cov[[i]][["EX2"]]),"ex2chr10=",sum(cov[[i]][["EX2"]]&cov[[i]][["CHR"]]=="10"),"\n")
  }

for(t in c(30,35,40)) 
  for(i in seq(3)) 
  {
    cov[[i]][["IN"]] = FALSE
    cov[[i]][["IN"]][findAneuploid(cov[[i]]$COV<thresh67[i],t)] = TRUE
    cov[[i]][["EX2"]] = cov[[i]][["EX"]] & (!cov[[i]][["IN"]])
    cat("f=",files[i],"t=",t,"ex2=",sum(cov[[i]][["EX2"]]),"ex2chr10=",sum(cov[[i]][["EX2"]]&cov[[i]][["CHR"]]=="10"),"\n")
  }

for(t in c(30,35,40)) 
  for(i in seq(3)) 
  {
    cov[[i]][["IN"]] = FALSE
    cov[[i]][["IN"]][findAneuploid(cov[[i]]$COV<thresh65[i],t)] = TRUE
    cov[[i]][["EX2"]] = cov[[i]][["EX"]] & (!cov[[i]][["IN"]])
    cat("f=",files[i],"t=",t,"ex2=",sum(cov[[i]][["EX2"]]),"ex2chr10=",sum(cov[[i]][["EX2"]]&cov[[i]][["CHR"]]=="10"),"\n")
  }

#t=35, thresh66, thresh16
#IRanges
#read in mole file
#filter out if (exclude)overlap(start) or (exclude)overlap(end)
#save new mole file

a=(findAneuploid(cov[[2]]$COV<thresh66[2],30))
cov[[2]][["IN"]]=FALSE
cov[[2]][["IN"]][a]=TRUE
cov[[2]][["EX2"]]=cov[[2]][["EX"]] & (!cov[[2]][["IN"]])
cov[[2]][cov[[2]]$CHR=="10" & cov[[2]]$START>93300000 & cov[[2]]$END <123900000 & cov[[2]]$EX2, ]
cov[[2]][cov[[2]]$EX2,]
aggregate(. ~ CHR, cov[[2]][cov[[2]]$EX2,])
