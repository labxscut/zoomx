#!/usr/bin/env Rscript

#add SV phasing to script use MolHP, MolPS

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

#for(p in c("devtools", "optparse", "rjson", "ggplot2", "RColorBrewer", "cowplot", "igraph", "raster","httr","plyr","stringr","brms")) require(p,character.only=T)
for(p in c("devtools", "optparse", "rjson", "ggplot2", "RColorBrewer", "cowplot", "igraph")) require(p,character.only=T)
#for(p in c("XML")) myrequire(p,version="3.96-0.2",repo="http://cran.us.r-project.org")
for(p in c("GenomicRanges", "ggbio", "Rsamtools", "Biostrings","BSgenome.Hsapiens.UCSC.hg19","BSgenome.Hsapiens.UCSC.hg38")) require(p,character.only=T)
#for(p in c("GenomicRanges", "ggbio", "Rsamtools", "Biostrings","BSgenome.Hsapiens.UCSC.hg19","BSgenome.Hsapiens.UCSC.hg38")) myrequire(p,repo="Bioc")
#for(p in c("primerTree")) myrequire(p,repo="https://github.com/chaelir/primerTree.git")

#genome plot packages: ggbio, Sushi, genoplotR, circlize 

####Options####
option_list2 <- list(
  make_option(c("-f", "--flushSD"), type="integer", default=5000,
              help="maximum of flushing sd [default %default]"),
  make_option(c("-g", "--gapSize"), type="integer", default=100000,
              help="2 times the boundary padding bp in molecule plot [default %default]"),
  make_option(c("-s", "--scaleSize"), type="integer", default=10000,
              help="boundary rounded to bp [default %default]"),
  make_option(c("-b", "--breakCI"), type="numeric", default=0.95,
              help="limit of CI to estimate break point [default %default]"),
  make_option(c("-r", "--reference"), default="hg37",
              help="use hg19 for UCSC build [default %default]"),
  make_option(c("-l", "--breakSlopSize"), type="integer", default=1000,
              help="2 times the boundary padding bp in read pair plot [default %default]"),
  make_option(c("-p", "--padding"), default="-",
              help="signature pattern for breakpoint [default %default]"),
  make_option(c("-m", "--includeDupMolecule"), action="store_true", default=FALSE,
              help="including duplicated molecules in the plots/tables [default %default]"),
  make_option(c("-u", "--isUCSC"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-e", "--minAlleleFraction"), type="numeric", default=0.05,
              help="Hypothesized minimum alleleFraction i.e. H0>minAlleleFraction for Bayesian Factor to detect [default %default]"),
  make_option(c("-w", "--withControl"), action="store_true", default=FALSE,
              help="with a control sample, [default %default]"),
  make_option(c("-y", "--withoutContig"), action="store_false", default=TRUE,
              help="generate rdp plots and contigs, [default %default]"),
  make_option(c("-o", "--plotPrefix"), default="none",
              help="output prefix [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)

parser <- OptionParser(usage = "%prog [options] workDir statsJSON", option_list=option_list2)
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)
if(length(cmd$args)!=2 & !cmd$options$debug){ print_help(parser); quit() }
cat("-Info: invoking command:",commandArgs(),"\n")
breakCI=cmd$options$breakCI
verbose=cmd$options$noQuiet
debug=cmd$options$debug
isUCSC=cmd$options$isUCSC
gapSize=cmd$options$gapSize
scaleSize=cmd$options$scaleSize
flushSD=cmd$options$flushSD
reference = cmd$options$reference
breakSlopSize = cmd$options$breakSlopSize
padding = cmd$options$padding
plotPrefix = cmd$options$plotPrefix
includeDupMolecule = cmd$options$includeDupMolecule
withControl = cmd$options$withControl
withoutContig = cmd$options$withoutContig
minAlleleFraction = cmd$options$minAlleleFraction
workDir=cmd$args[1]
statsFile=cmd$args[2]

####Debug####
if(debug){ 
  breakCI=0.95
  gapSize = 100000     #boundary padding is 10K/2 
  scaleSize = 10000    #boundary rounded to 10K
  flushSD = 5000
  isUCSC=FALSE         #input has no chr in chromosome name
  workDir="~/setup/zoomx/test/plots/MetB7175.wg12.chr8-98M-100M.bam.dist.E10_Norm7176.wg12.chr8-98M-100M.bam.dist.E10.somatic.Plot_1"
  statsFile="~/setup/zoomx/test/data/MetB7175.wg12.chr8-98M-100M.bam.gridScan.cfg.json"
  minAlleleFraction=0.05
  withControl=TRUE
  #isUCSC=TRUE
  #workDir="~/setup/zoomx/test/plots/NA12878.wg12.chr22-22M-24M.bam.dist.E45.refined.Plot_1"
  #statsFile="~/setup/zoomx/test/data/NA12878.wg12.chr22-22M-24M.bam.gridScan.cfg.json"
  #withControl=FALSE
  includeDupMolecule=FALSE
  plotPrefix="none"
  reference="hg37"
  padding = "-"
  breakSlopSize = 1000
} else { 
  options(warn=-1)
}

#####Main#####

####plotMolecule-LoadInput####
options(stringsAsFactors = FALSE)
statsTab=rjson::fromJSON(paste(readLines(statsFile), collapse=""))
setwd(workDir)
regCols=c("chr1","start1","end1","chr2","start2","end2","name","barcnt","FR","RF","FF","RR","MI") # in refined
bedpeCols=c("chr1","st1","ed1","chr2","st2","ed2","hp")
bedpeClass=c("character","numeric","numeric","character","numeric","numeric","numeric")
molCols=c("chr1","st1","ed1","barID","molSize","molLen","barMolNum","molID","molHP","molPS","molBX")
# correct fields in BAM file { "molID" : "MI:i:", "molPS" : "PS:Z:", "molHP" : "HP:i:", "molBX" : "BX:Z:" }
molColsOld=c("chr1","st1","ed1","barID","molSize","molLen","barMolNum","molID","molHP","molPS")  #for backward compatibility

#read outs
if(withControl){
  suffixes = c("case", "control")
} else {
  suffixes = c("case")
}
message("-Info: suffixes=",paste(suffixes, collapse=","))

strandDf=c("1","-1",".")
strandGR=c("+","-","*")
supportType=c("molSOV","molNSOV","rdpFR","rdpRF","rdpFF","rdpRR")
supportGR=c("molSOV","molNSOV","rdpFR","rdpRF","rdpFF","rdpRR")
strand1Info=c(".",".","1","-1","1","-1")
strand2Info=c(".",".","-1","1","1","-1")
supportColors=RColorBrewer::brewer.pal(3,"Set1")[1:length(supportGR)]
names(supportColors)=supportGR
strandColors=RColorBrewer::brewer.pal(3,"Set2")
names(strandColors)=strandGR

refined=tryCatch(read.table("refined.bedpe",sep="\t",col.names=regCols),error=function(e) {NULL})
infoTab=refined
chr1 = unique(infoTab$chr1)[1]
chr2 = unique(infoTab$chr2)[1]
if(plotPrefix=="none") {
  plotPrefix=paste(refined[1:6],collapse="_")
} else {
  plotPrefix=plotPrefix
}

dataReadIn = function(i, suffixes) {
  message("-Info: processing", suffixes[i])
  s = suffixes[i]
  refinedReg1=tryCatch(read.table(paste("refined.reg1.mol",s,sep="."),col.names=molCols),error=function(e) {t=read.table(paste("refined.reg1.mol",s,sep="."),col.names=molColsOld);t[["molBx"]]=NA;t}, finally={NULL})
  refinedReg2=tryCatch(read.table(paste("refined.reg2.mol",s,sep="."),col.names=molCols),error=function(e) {t=read.table(paste("refined.reg2.mol",s,sep="."),col.names=molColsOld);t[["molBx"]]=NA;t}, finally={NULL})
  totalMol = nrow(refinedReg1)+nrow(refinedReg2)
  sharedBarID = unique(refinedReg1$barID[refinedReg1$barID %in% refinedReg2$barID]) #shared barcodes
  featureMol=length(sharedBarID)
  mergedMolList = lapply(sharedBarID, function(b) { #find molecules sharing the same barcode between reg1 and reg2
    sel1=which(refinedReg1$barID==b)
    sel2=which(refinedReg2$barID==b)
    mdat=cbind(refinedReg1[sel1[1],c("chr1","st1","ed1")], refinedReg2[sel2[1],c("chr1","st1","ed1")])
    colnames(mdat)=bedpeCols[1:6]
    mdat$id=paste("molID",b,sep=":")
    mdat$type="molecule"
    mdat$barcode=refinedReg1[sel1[1],c("molBX")]
    mdat$mi1=refinedReg1[sel1[1],c("molID")]
    mdat$mi2=refinedReg2[sel2[1],c("molID")]
    mdat$hp1=refinedReg1[sel1[1],c("molHP")]
    mdat$hp2=refinedReg1[sel2[1],c("molHP")]
    mdat$ps1=refinedReg1[sel1[1],c("molPS")]
    mdat$ps2=refinedReg1[sel2[1],c("molPS")]
    mdat
  })
  mergedMol = do.call("rbind", mergedMolList)                         #all molecules, i.e. BX
  mergedMI = mergedMol[mergedMol$mi1 != mergedMol$mi2, ]              #keep non-self-overlapping molecules
  mergedMI = as.data.frame(mergedMI)
  if(nrow(mergedMI)>0) mergedMI$type = "molNSOV" 
  message("-Info: mergedMI:", nrow(mergedMI))
  mergedBI = mergedMol[mergedMol$mi1 == mergedMol$mi2, ]              #keep self-overlapping molecules
  mergedBI = as.data.frame(mergedBI)
  if(nrow(mergedBI)>0) mergedBI$type = "molSOV"
  message("-Info: mergedBI:", nrow(mergedBI))
  
  # add "chr" for non UCSC format, required step to use BSgenome
  message("-Info: refined:", paste(refined))
  if(!isUCSC & nrow(refined)>0) {refined$chr1=paste("chr",refined$chr1,sep="")}
  if(!isUCSC & nrow(refined)>0) {refined$chr2=paste("chr",refined$chr2,sep="")}
  if(!isUCSC & nrow(refinedReg1)>0) {refinedReg1$chr1=paste("chr",refinedReg1$chr1,sep="")}
  if(!isUCSC & nrow(refinedReg2)>0) {refinedReg2$chr1=paste("chr",refinedReg2$chr1,sep="")}
  
  mergedRF=NULL; mergedRR=NULL; mergedFR=NULL; mergedFF=NULL; supportRDP=0
  
  #get rdp information for case file only
  if(s=="case") {
    refinedFRid=rjson::fromJSON(refined$FR)
    refinedRFid=rjson::fromJSON(refined$RF)
    refinedFFid=rjson::fromJSON(refined$FF)
    refinedRRid=rjson::fromJSON(refined$RR)
    supportRDP=length(refinedFRid)+length(refinedRFid)+length(refinedFFid)+length(refinedRRid)

    refinedFR=tryCatch(read.table(paste("refined.FR.bedpe",s,sep="."),col.names=bedpeCols, colClasses=bedpeClass),error=function(e) {NULL})
    refinedRF=tryCatch(read.table(paste("refined.RF.bedpe",s,sep="."),col.names=bedpeCols, colClasses=bedpeClass),error=function(e) {NULL})
    refinedFF=tryCatch(read.table(paste("refined.FF.bedpe",s,sep="."),col.names=bedpeCols, colClasses=bedpeClass),error=function(e) {NULL})
    refinedRR=tryCatch(read.table(paste("refined.RR.bedpe",s,sep="."),col.names=bedpeCols, colClasses=bedpeClass),error=function(e) {NULL})
    refinedMI=tryCatch(read.table(paste("refined.FR.bedpe",s,sep="."),col.names=bedpeCols, colClasses=bedpeClass),error=function(e) {NULL})

    #following assume (chr1==chr2) => st1 <= st2 in rdp output which is enforced by parse10x
    mergedFR = do.call("rbind", 
                       apply(refinedFR[,1:6], 1, FUN = function(r, chr1) { 
                         if( r[1]!=chr1 | (r[4]==chr1 & as.numeric(r[2])>as.numeric(r[5])) ) 
                         { t=r[1:3]; r[1:3]=r[4:6]; r[4:6]=t; }
                         as.data.frame(t(r)) 
                      }, chr1) )
    mergedFR = as.data.frame(mergedFR)
    if(nrow(mergedFR)>0) { #set_colClass(mergedFR, c("character","integer","integer","character","integer","integer"));
      mergedFR$id = paste("rdpFR",refinedFRid,sep=":"); mergedFR$type="rdpFR"; mergedFR$barcode=refinedFR$hp;
      mergedFR$mi1 = 0; mergedFR$mi2 = 0; mergedFR$hp1 = 0; mergedFR$hp2 = 0;}
    message("-Info: mergedFR:", nrow(mergedFR))
    
    mergedRF = do.call("rbind", 
                       apply(refinedRF[,1:6], 1, FUN = function(r, chr1) { 
                         if( r[1]!=chr1 | (r[4]==chr1 & as.numeric(r[2])>as.numeric(r[5])) ) 
                         { t=r[1:3]; r[1:3]=r[4:6]; r[4:6]=t; }
                         as.data.frame(t(r)) 
                       }, chr1) )
    mergedRF = as.data.frame(mergedRF)
    if(nrow(mergedRF)>0) { #set_colClass(mergedFR, c("character","integer","integer","character","integer","integer"));
      mergedRF$id = paste("rdpRF",refinedRFid,sep=":"); mergedRF$type="rdpRF"; mergedRF$barcode=refinedRF$hp;
      mergedRF$mi1 = 0; mergedRF$mi2 = 0; mergedRF$hp1 = 0; mergedRF$hp2 = 0;}
    message("-Info: mergedRF:", nrow(mergedRF))
    
    mergedFF = do.call("rbind", 
                       apply(refinedFF[,1:6], 1, FUN = function(r, chr1) { 
                         if( r[1]!=chr1 | (r[4]==chr1 & as.numeric(r[2])>as.numeric(r[5])) ) { t=r[1:3]; r[1:3]=r[4:6]; r[4:6]=t; }
                                                                                      as.data.frame(t(r)) }, chr1))
    mergedFF = as.data.frame(mergedFF)
    if(nrow(mergedFF)>0) { #set_colClass(mergedFF, c("character","integer","integer","character","integer","integer")); 
                             mergedFF$id = paste("rdpFF",refinedFFid,sep=":"); mergedFF$type="rdpFF"; mergedFF$barcode=refinedFF$hp;
                             mergedFF$mi1 = 0; mergedFF$mi2 = 0; mergedFF$hp1 = 0; mergedFF$hp2 = 0;} 
    message("-Info: mergedFF:", nrow(mergedFF))
    
    mergedRR = do.call("rbind", 
                       apply(refinedRR[,1:6], 1, FUN = function(r, chr1) { 
                         if( r[1]!=chr1 | (r[4]==chr1 & as.numeric(r[2])>as.numeric(r[5])) ) { t=r[1:3]; r[1:3]=r[4:6]; r[4:6]=t; } 
                       as.data.frame(t(r)) }, chr1))
    mergedRR = as.data.frame(mergedRR)
    if(nrow(mergedRR)>0) { #set_colClass(mergedRR, c("character","integer","integer","character","integer","integer")); 
                             mergedRR$id = paste("rdpRR",refinedRRid,sep=":"); mergedRR$type="rdpRR"; mergedRR$barcode=refinedRR$hp;
                             mergedRR$mi1 = 0; mergedRR$mi2 = 0; mergedRR$hp1 = 0; mergedRR$hp2 = 0;}
    message("-Info: mergedRR:", nrow(mergedRR))
  
    if(!isUCSC & nrow(mergedFR)>0) {mergedFR$chr1=paste("chr",mergedFR$chr1,sep="")
                                     mergedFR$chr2=paste("chr",mergedFR$chr2,sep="")}
    if(!isUCSC & nrow(mergedRF)>0) {mergedRF$chr1=paste("chr",mergedRF$chr1,sep="")
                                      mergedRF$chr2=paste("chr",mergedRF$chr2,sep="")}
    if(!isUCSC & nrow(mergedFF)>0) {mergedFF$chr1=paste("chr",mergedFF$chr1,sep="")
                                     mergedFF$chr2=paste("chr",mergedFF$chr2,sep="")}
    if(!isUCSC & nrow(mergedRR)>0) {mergedRR$chr1=paste("chr",mergedRR$chr1,sep="")
                                     mergedRR$chr2=paste("chr",mergedRR$chr2,sep="")}
    if(!isUCSC & nrow(mergedMI)>0) {mergedMI$chr1=paste("chr",mergedMI$chr1,sep="")
                                       mergedMI$chr2=paste("chr",mergedMI$chr2,sep="")}
    if(!isUCSC & nrow(mergedBI)>0) {mergedBI$chr1=paste("chr",mergedBI$chr1,sep="")
                                      mergedBI$chr2=paste("chr",mergedBI$chr2,sep="")}

  } #end of getting all rdps
   
  mergedAll=rbind(mergedFR,mergedRF,mergedFF,mergedRR,mergedBI,mergedMI)
  
  if(nrow(mergedAll)>0) { #plotMolecule-
    mergedAll$strand1=strand1Info[match(mergedAll$type, supportType)]
    mergedAll$strand2=strand2Info[match(mergedAll$type, supportType)]
    mergedAll$supportType=mergedAll$type
    GRcols=c("chrom1",      "start1", "end1",   "chrom2",   "start2", "end2",   "name",     "score",    "barcode",  "mi1",      "mi2",      "hp1",      "hp2",      "ps1", "ps2", "strand1",  "strand2",  "supportType")
    GRclasses=c("character","numeric","numeric","character","numeric","numeric","character","character","character","character","character","character","character","character", "character", "character","character","character")
    colnames(mergedAll) = GRcols
    mergedAll$score = 1
    mergedAll = set_colClass(mergedAll, GRclasses)
  }

  print(mergedAll)

  if(!includeDupMolecule) {  
    sel = (mergedAll$supportType != "molSOV") #ecluding self-overlapping molecules
  }
  mergedAll=mergedAll[sel,]
  message("-Info: mergedAll:", nrow(mergedAll))

  rawReg1St = min(refined$start1,refinedReg1$st1)
  rawReg1Ed = max(refined$end1,refinedReg1$ed1)
  rawReg2St = min(refined$start2,refinedReg2$st1)
  rawReg2Ed = max(refined$end2,refinedReg2$ed1)
  scaledReg1St = floor(rawReg1St/scaleSize)*scaleSize - gapSize/2
  scaledReg1Ed = ceiling(rawReg1Ed/scaleSize)*scaleSize + gapSize/2
  scaledReg2St = floor(rawReg2St/scaleSize)*scaleSize - gapSize/2
  scaledReg2Ed = ceiling(rawReg2Ed/scaleSize)*scaleSize + gapSize/2
  boundaryReg12 = c(scaledReg1St,scaledReg1Ed,scaledReg2St,scaledReg2Ed)
  shiftReg1St = 1    #Mark Reg1 st point 1 as (chr1: scaledReg1St)
  shiftReg1Ed = scaledReg1Ed - scaledReg1St #Reg1 ed point scaledReg1Ed-scaledReg2St as (chr1: scaledReg1Ed) 
  shiftReg2St = shiftReg1Ed + 1   #Mark Reg2 st point shiftReg1Ed+1 as (chr2: scaledReg2St)
  shiftReg2Ed = shiftReg1Ed + (scaledReg2Ed - scaledReg2St) #Mark as (chr2: scaledReg2Ed)
  shiftReg12 = c(shiftReg1St,shiftReg1Ed,shiftReg2St,shiftReg2Ed)
  breakpoint1 = NA
  breakpoint2 = NA
  leftBreakReg1 = NA
  leftBreakReg2 = NA
  breakReg1CILow = NA
  breakReg1CIHigh = NA
  breakReg2CILow = NA
  breakReg2CIHigh = NA
  ### The LOGIC here is temporary cause obvisouly even if the junction is the same hap
  ### the hap doesn't expect to be correctly identified  by long-ranger correclty
  ### let's break it down to haplotypes and phased segemnts at each breakpoint
  hapTab1 = table(mergedAll$hp1[mergedAll$hp1!=-1])
  hapTab2 = table(mergedAll$hp2[mergedAll$hp2!=-1]) #only hp=1 and hp=2 matters
  cat("-Info: haplotype (hp) table of breakpoint1"); print(hapTab1)
  cat("-Info: haplotype (hp) table of breakpoint2"); print(hapTab2)
  breakHapTab1 = rjson::toJSON(hapTab1)
  breakHapTab2 = rjson::toJSON(hapTab2)
  breakHapMax =  rep(0, 2)
  for(i in 1:2) {
    hapTab = list(hapTab1, hapTab2)[[i]]
    if(!is.null(hapTab)) 
      breakHapMax[i]=names(sort(hapTab, decreasing=TRUE))[1]
    cat("-Info: breakpoint", i, " haplotype"); print(breakHapMax[i])
  }

  breakPsMax =  rep(0, 2)
  psTab1 = table(mergedAll$ps1[mergedAll$ps1!=-1])
  psTab2 = table(mergedAll$ps2[mergedAll$ps2!=-1])
  cat("-Info: phased segment (ps) of breakpoint1"); print(psTab1)
  cat("-Info: phased segment (ps) of breakpoint2"); print(psTab2)
  breakPsTab1 = rjson::toJSON(psTab1)
  breakPsTab2 = rjson::toJSON(psTab2)
  breakPsMax =  rep(0, 2)
  for(i in 1:2) {
    psTab = c(psTab1, psTab2)[i]
    if(!is.null(psTab)) 
      breakPsMax[i]=names(sort(psTab, decreasing=TRUE))[1]
    cat("-Info: breakpoint ", i, " phased segment"); print(breakPsMax[i])
  }

  if(nrow(mergedAll)>0) {
    convertedAll = mergedAll  #molecule only
    convertedAll$chrom1 = "chrA"
    convertedAll$start1 = convertedAll$start1 - scaledReg1St + 1
    convertedAll$end1 = convertedAll$end1 - scaledReg1St + 150
    convertedAll$chrom2 = "chrA"
    convertedAll$start2 = convertedAll$start2 - scaledReg2St + shiftReg2St + 1
    convertedAll$end2 = convertedAll$end2 - scaledReg2St + shiftReg2St +150
    convertedAll$supportType = supportType[match(convertedAll$supportType,supportType)]
    convertedAll$strand1 = strandGR[match(convertedAll$strand1,strandDf)]
    convertedAll$strand2 = strandGR[match(convertedAll$strand2,strandDf)]
    flushSDReg1 = c(sd(convertedAll$start1), sd(convertedAll$end1), flushSD+1) #avoid NA in min
    flushSDReg2 = c(sd(convertedAll$start2), sd(convertedAll$end2), flushSD+1) #avoid NA in min
    leftBreakReg1 = if(min(flushSDReg1)<flushSD) TRUE else NA  #determine the orientation of breakpoint
    leftBreakReg2 = if(min(flushSDReg2)<flushSD) TRUE else NA
    if(!is.na(leftBreakReg1)) leftBreakReg1 = flushSDReg1[1] < flushSDReg1[2]
    if(!is.na(leftBreakReg2)) leftBreakReg2 = flushSDReg2[1] < flushSDReg2[2]
    # leftBreak: TRUE iff left edge of Reg1 is break, FALSE if right edge, NA if break is middle
    # breakpoint inferred using all molecules and read pairs available
    breakReg1 = round(if(is.na(leftBreakReg1)) { median(c(convertedAll$start1,convertedAll$end1)) } else if (leftBreakReg1) { quantile(convertedAll$start1,1-breakCI) } else { quantile(convertedAll$end1,breakCI) })
    breakReg1CILow = min(convertedAll$start1)
    breakReg1CIHigh = max(convertedAll$end1)
    breakReg2 = round(if(is.na(leftBreakReg2)) { median(c(convertedAll$start2,convertedAll$end2)) } else if (leftBreakReg2) { quantile(convertedAll$start2,1-breakCI) } else { quantile(convertedAll$end2,breakCI) })
    breakReg2CILow = min(convertedAll$start2)
    breakReg2CIHigh = max(convertedAll$end2)

    message("-Info: gr...")
    #print(convertedAll[convertedAll$end1>convertedAll$start2,c(1,2,3,4,5,6)])
    #print(mergedAll[convertedAll$end1>convertedAll$start2,c(1,2,3,4,5,6)])
    gr <- GenomicRanges::GRanges(seqnames = c(convertedAll$chrom1, convertedAll$chrom2), 
                ranges = IRanges::IRanges(start = c(convertedAll$start1, convertedAll$start2), 
                                 end = c(convertedAll$end1, convertedAll$end2)),
                #strand = factor(c(as.character(convertedAll$strand1),as.character(convertedAll$strand2)), levels=strandGR),
                strand = c(convertedAll$strand1,convertedAll$strand2),
                value = rep(0,2*nrow(convertedAll)), 
                score = c(convertedAll$score, convertedAll$score), 
                supportType = c(convertedAll$supportType, convertedAll$supportType), 
                pair = c(convertedAll$name, convertedAll$name))
    message("-Info: grComp...")
    #print(convertedAll)
    grComp <- GenomicRanges::GRanges(seqnames = convertedAll$chrom1,
                    ranges = IRanges::IRanges(start = convertedAll$end1, end=pmax(convertedAll$end1,convertedAll$start2)+1),
                    strand = rep(NA,nrow(convertedAll)),
                    value = rep(0,nrow(convertedAll)),
                    score = rep(0,nrow(convertedAll)),
                    supportType = convertedAll$supportType,
                    pair = convertedAll$name)

    tickReg1 = c(shiftReg12[1]+gapSize/2, shiftReg12[2]-gapSize/2) 
    tickReg2 = c(shiftReg12[3]+gapSize/2, shiftReg12[4]-gapSize/2)
    tickLab1 = c(boundaryReg12[1]+gapSize/2, boundaryReg12[2]-gapSize/2) 
    tickLab2 = c(boundaryReg12[3]+gapSize/2, boundaryReg12[4]-gapSize/2)
    tickLab1 = paste(chr1,tickLab1,sep=":")
    tickLab2 = paste(chr2,tickLab2,sep=":")
    
    breakpoint1 = boundaryReg12[1]+unname(breakReg1[1])
    breakpoint2 = boundaryReg12[3]+unname(breakReg2[1])-shiftReg12[2]
    
    message("-Info: plotting...")
    message("-Info: plotFile=",paste(workDir, paste(plotPrefix,"break",s,"pdf",sep="."),sep="/"))
    pdf(file=paste(plotPrefix,"break",s,"pdf",sep="."), width=16, height=12)
    #print(gr)
    #save(gr, file='/home/xlib/setup/gr.debug.rdata')
    p = ggplot2::ggplot()+ggbio::geom_arrowrect(gr, ggplot2::aes(color=supportType, fill=strand, group="pair"))
    #print(grComp)
    p = p + ggbio::geom_segment(grComp, ggplot2::aes(color=supportType, group="pair"))
    p = p + ggplot2::scale_color_manual(values=supportColors, na.value="black")
    p = p + ggplot2::scale_fill_manual(values=strandColors, na.value="black")
    p = p + ggplot2::scale_x_continuous("Genomic Region Pair", breaks=tickReg1, labels=tickLab1,
                                        limits=c(shiftReg12[1],shiftReg12[4]),
                                        sec.axis=ggplot2::sec_axis(~., breaks=tickReg2, labels=tickLab2))
    p = p + ggplot2::geom_vline(xintercept=breakReg1, linetype="dotted") + ggplot2::geom_vline(xintercept=breakReg2, linetype="dotted")
    p = p + ggplot2::geom_text(ggplot2::aes(x=breakReg1-gapSize/10, label=paste(breakpoint1), y=100), colour="blue", angle=90)
    p = p + ggplot2::geom_text(ggplot2::aes(x=breakReg2+gapSize/10, label=paste(breakpoint2), y=100), colour="blue", angle=90)
    #p = p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0, hjust = 1, color="red"))
    print(p)
    dev.off()
  }
  
  ####WriteOutput####

  message("-Info: exporting infoTab...")
  #simply genotyping by short readpair, requiring two read pairs at least
  if(infoTab$chr1!=infoTab$chr2) {
    tmpSVType = "TRA"
  } else if(sum(grepl("rdp",mergedAll$supportType))>0) {
   tmpSVType = "" #default deletion
   if(sum(mergedAll$supportType %in% c("rdpRR","rdpFF"))>1) tmpSVType = paste("INV",tmpSVType)
   if(sum(mergedAll$supportType %in% c("rdpRF"))>=1) tmpSVType = paste("DUP",tmpSVType)
   if(sum(mergedAll$supportType %in% c("rdpFR"))>=1) tmpSVType = paste("DEL",tmpSVType)
   if(tmpSVType == "") tmpSVType="UNK"
  } else {
    tmpSVType = "UNK"
  }
  
  infoTabNew = cbind(infoTab[,seq(7)], as.data.frame( list( breakpoint1=breakpoint1,
                                                            breakpoint2=breakpoint2,
                                                            isLeftBreak1=leftBreakReg1,
                                                            isLeftBreak2=leftBreakReg2,
                                                            SVHap=paste(breakHapMax,sep="<>"),
                                                            SVPs=paste(breakPsMax,sep="<>"),,
                                                            SVType=tmpSVType,
                                                            AlleleFraction=sum(mergedAll$supportType=="molNSOV")/statsTab$Cem,
                                                            ReadPair=supportRDP,
                                                            MolCnt=sum(mergedAll$supportType=="molNSOV"),
                                                            totalMol=totalMol,
                                                            featureMol=featureMol,
                                                            Sample=s)) )

  mergedAllFile = paste(plotPrefix,"support",s,"txt",sep=".")
  write.table(mergedAll,file=mergedAllFile,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE) #mergeAll
  
  jsonFile=base::file(paste(plotPrefix,"break",s,"json",sep="."))
  writeLines(rjson::toJSON(list(reg1=chr1,
                                break1=infoTabNew$breakpoint1,
                                break1CILow=boundaryReg12[1]+unname(breakReg1CILow[1]),
                                break1CIHigh=boundaryReg12[1]+unname(breakReg1CIHigh[1]),
                                isLeftBreak1=leftBreakReg1,
                                isLeftBreak2=leftBreakReg2,
                                reg2=chr2,
                                break2=infoTabNew$breakpoint2, 
                                break2CILow=boundaryReg12[3]+unname(breakReg2CILow[1])-shiftReg12[2],
                                break2CIHigh=boundaryReg12[3]+unname(breakReg2CIHigh[1])-shiftReg12[2],
                                sharedBarcode=mergedMol$barcode,
                                MolCnt=sum(mergedAll$supportType=="molNSOV"),
                                SVType=tmpSVType,
                                SVHap=paste(breakHapMax,sep="<>"),
				SVHapTab=list(hapTab1=hapTab1,hapTab2=hapTab2),
                                SVPs=paste(breakPsMax,sep="<>"),,
                                SVPsTab=list(hapTab1=hapTab1,hapTab2=hapTab2),
                                AlleleFraction=sum(mergedAll$supportType=="molNSOV")/statsTab$Cem,
                                ReadPair=supportRDP,
                                Sample=s
  )), con=jsonFile)
  close(jsonFile)
 
  return(list(mergedAll=mergedAll, infoTab=infoTabNew)) 
}

allData = lapply(seq_along(suffixes), dataReadIn, suffixes)
names(allData)=suffixes

#chr1  start1	end1	chr2	start2	end2	molcnt	breakpoint1	breakpoint2	SVType	AlleleFraction	ReadPair Sample
infoTabNew = Reduce(rbind, lapply(allData,'[[',"infoTab"))
infoFile = paste(plotPrefix,"info","txt",sep=".")
infoTabOutput=infoTabNew
if(length(suffixes)==2) { 
  #do the fisher exact test here, second row is control
  zC=infoTabNew[2,]$featureMol
  kC=infoTabNew[2,]$totalMol-2*infoTabNew[2,]$featureMol #nC-zC
  zT=infoTabNew[1,]$featureMol
  kT=infoTabNew[1,]$totalMol-1*infoTabNew[1,]$featureMol #nC-zC
  t=fisher.test(rbind(c(zC,zT),c(kC,kT)), alternative="less")
  infoTabOutput$fisherP=t$p.value
  #cleaning up for case/control output
  infoTabOutput[2,]$breakpoint1=NA
  infoTabOutput[2,]$breakpoint2=NA
  infoTabOutput[2,]$isLeftBreak1=NA
  infoTabOutput[2,]$isLeftBreak2=NA
  infoTabOutput[2,]$SVType=NA
  infoTabOutput[2,]$fisherP=NA
}

message("-Info: infoFile=",workDir,"/",infoFile); 
write.table(infoTabOutput,file=infoFile,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)
message("-Info: infoTabOutput:")
print(infoTabOutput)
if(withoutContig) { cat("PLOTMOLECULEDONE in",proc.time()[3]-ptm[3],"seconds! \n"); quit() }

####rdpBreakpoint-PlotRdp####
#case only
callInfo=refined
colScheme="Dark2"
colors=RColorBrewer::brewer.pal(length(supportType), colScheme)
names(colors)=supportType
info=allData[["case"]][["mergedAll"]]
infoTab=allData[["case"]][["infoTab"]]
info=info[!grepl("mol",info$supportType),] #only considers rdp

if(reference=="hg38") {
  refSeq = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
} else { #fall back to hg19/hg37
  refSeq = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
}
strandDf=c("1","-1",".")
strandGR=c("+","-","*")
rdpInfo=info
strandSymbols=list("+"=1,"*"=".","-"=-1)
rdpInfo$strand1=strandGR[match(rdpInfo$strand1,strandDf)]
rdpInfo$strand2=strandGR[match(rdpInfo$strand2,strandDf)]
rd1GR=GenomicRanges::GRanges(seqnames=rdpInfo$chrom1, 
                             ranges=IRanges::IRanges(start=as.numeric(rdpInfo$start1),end=as.numeric(rdpInfo$end1)),
                             name=rdpInfo$name,
                             score=as.numeric(rdpInfo$score),
                             strand=rdpInfo$strand1,
                             type=rdpInfo$supportType)

rd2GR=GenomicRanges::GRanges(seqnames=rdpInfo$chrom2, 
                             ranges=IRanges::IRanges(start=as.numeric(rdpInfo$start2),end=as.numeric(rdpInfo$end2)),
                             name=rdpInfo$name,
                             score=as.numeric(rdpInfo$score),
                             strand=rdpInfo$strand2,
                             type=rdpInfo$supportType)

message("-Info: rd1GR"); print(rd1GR)
message("-Info: rd2GR"); print(rd2GR)
if(length(rd1GR)==0 | length(rd2GR)==0) {
  message("-Info: no supporting read pairs and/or molecules, quitting")
  quit()
}

message("-Info: plotting rdps...")
g1 = ggplot2::ggplot() + 
  ggbio::geom_arrowrect(data=rd1GR,ggplot2::aes(y=name,group="name",fill=type,group.selfish=TRUE),xlab=as.character(GenomicRanges::seqnames(rd1GR)[1])) + 
  ggplot2::theme(legend.position="bottom", legend.box="vertical", legend.direction = "horizontal") + #box-v, item-h
  ggplot2::geom_vline(xintercept=infoTab$breakpoint1,linetype="dotted") +
  ggplot2::geom_text(ggplot2::aes(x=infoTab$breakpoint1, label=paste(infoTab$breakpoint1), y=1), colour="blue", angle=90)
  
g2 = ggplot2::ggplot() + 
  ggbio::geom_arrowrect(data=rd2GR,ggplot2::aes(y=name,group="name",fill=type,group.selfish=TRUE),xlab=as.character(GenomicRanges::seqnames(rd2GR)[1])) +
  ggplot2::geom_vline(xintercept=infoTab$breakpoint2,linetype="dotted") +
  ggplot2::geom_text(ggplot2::aes(x=infoTab$breakpoint2, label=paste(infoTab$breakpoint2), y=1), colour="blue", angle=90)

#ggplot2::ggsave("tmp1.pdf", g1)
#ggplot2::ggsave("tmp2.pdf", g2)

legendHeight=0.3 #control legend area height
plots=list(g1,g2)
plotGrobs1=lapply(plots,ggbio:::Grob)
message("-Info: gr1...")
legend <- cowplot::get_legend(plotGrobs1[[1]])
plotGrobs2=lapply(plots,function(p) { ggbio:::Grob(p+ggplot2::theme(legend.position="none")) })
message("-Info: gr2...")
if(length(plotGrobs2)>1) { labels="AUTO" } else { labels=NULL }
prow <- cowplot::plot_grid( plotlist = plotGrobs2, labels = labels, nrow = 1 )
message("-Info: grid...")
lrow <- cowplot::ggdraw() + cowplot::draw_grob(legend) #this only draws the legend, all aes has to be predefined.
pall = cowplot::plot_grid( prow, lrow, ncol=1, rel_heights=c(1,legendHeight))
message("-Info: saving support plots to:", paste(plotPrefix,"rdp","case","pdf",sep=".") )
cowplot::save_plot(paste(plotPrefix,"rdp","case","pdf",sep="."), pall, base_height = 10, base_aspect_ratio = 2.1)

####rdpBreakpoint-InferContig####

#convert support groups to contigs

#estimate breakpoint based on each supportType group (same supportType group)
#break1IsLeft=rep(NA,length(supportType))
#break2IsLeft=rep(NA,length(supportType))

# contruct a graph of |rdps|*|rdps| vertices
# c1 & c2  => c
# T  &  F  => F
# T  &  T  => T

#we first group by type
#  in each group by type, 
#  debug not finished
direction1=ifelse(is.na(infoTab$leftBreakReg1), "U", ifelse(infoTab$leftBreakReg1, "L", "R"))
direction2=ifelse(is.na(infoTab$leftBreakReg2), "U", ifelse(infoTab$leftBreakReg2, "L", "R"))
contigs=list()
for(i in seq_along(supportType)) {
  # i=1; i=2; i=3; i=4
  s = supportType[i]
  if( ! grepl("mol", s) ) { #some rdp support
    if(s=="rdpFF") {
      strand1 = "+"; strand2 = "+"
    } else if(s=="rdpRR") {
      strand1 = "-"; strand2 = "-"
    } else if(s=="rdpFR") {
      strand1 = "+"; strand2 = "-"
    } else if(s=="rdpRF") {
      strand1 = "-"; strand2 = "+"
    } else {
      next
      #strand1 = c("-","+"); strand2 = c("+","-")
    }
    rdps = info[info$supportType==s,]
    if(nrow(rdps)==0) next #no such rdp support
    rdpRanges1 = IRanges(start=rdps$start1, end=rdps$start1+1)
    rdpRangePairs1 = expand.grid(rdpRanges1,rdpRanges1) #
    rdpRangePairOverlaps1 = matrix(poverlaps(rdpRangePairs1[,1], rdpRangePairs1[,2], maxgap=breakSlopSize/2), nrow=length(rdpRanges1) )
    rdpRanges2 = IRanges(start=rdps$start2, end=rdps$start2+1)
    rdpRangePairs2 = expand.grid(rdpRanges2,rdpRanges2) #
    rdpRangePairOverlaps2 = matrix(poverlaps(rdpRangePairs2[,1], rdpRangePairs2[,2], maxgap=breakSlopSize/2), nrow=length(rdpRanges2) )
    rdpRangePairOverlaps = rdpRangePairOverlaps1 & rdpRangePairOverlaps2
    colnames(rdpRangePairOverlaps)=seq(1,dim(rdpRangePairOverlaps)[2]) #make sure cols and rows has the same set of names
    rownames(rdpRangePairOverlaps)=seq(1,dim(rdpRangePairOverlaps)[1])
    rdpRangePairOverlapsGraph = igraph::graph_from_adjacency_matrix(rdpRangePairOverlaps) #|rdps| by |rdps| graph
    rdpRangePairOverlapsGraphComponents=igraph::components(rdpRangePairOverlapsGraph)   #will assume 2*|rdps| nodes
    for(j in seq.int(rdpRangePairOverlapsGraphComponents$no)) { 
      sel=which(rdpRangePairOverlapsGraphComponents$membership==j) #only need first half
      selrdps=rdps[sel,]
      #if(nrow(selrdps)<2) next #require at least two rdps 
      #if less than three rdps, we cannot determine direction
      #if more, we determine direction by find small difference size
      #break1IsLeft is NA if neither Left or Right, is False if Right
      # this part is not reliable, always assume not knowing breakpoint direction
      # if(nrow(selrdps)>2) { #rdp >=3, determine direction by read mapping density
      #   start1Diff=diff(sort(selrdps$start1)) #sorted
      #   if(start1Diff[1]<tail(start1Diff,1)) direction1 = "L"
      #   if(start1Diff[1]>tail(start1Diff,1)) direction1 = "R"
      #   start2Diff=diff(sort(selrdps$start2)) #sorted
      #   if(start2Diff[1]<tail(start2Diff,1)) direction2 = "L"
      #   if(start2Diff[1]>tail(start2Diff,1)) direction2 = "R"
      # }
      sel1 = paste(strand1,direction1) %in% c("+ R", "- L", "+ U", "- U") #possible read1 configuration
      sel2 = paste(strand2,direction2) %in% c("+ R", "- L", "+ U", "- U") #possible read2 configuration
      sel = sel1 & sel2                                                   #only select possible ones
      #strand1 = strand1[sel]; strand2 = strand2[sel];                     
      direction1 = ifelse(strand1=="+","R","L"); direction2 = ifelse(strand2=="+","R","L"); 
      #once know the strand convert U to L/R based on +/-
      # 1R2R, 1L2R, 1R2L, 1L2L
      start1 = (direction1=="R")*max(selrdps$start1) + (direction1=="L")*min(selrdps$start1); end1=start1
      start2 = (direction2=="R")*max(selrdps$start2) + (direction2=="L")*min(selrdps$start2); end2=start2
      for(k in seq_along(start1)) { # handles the two possibility of unknown type read pairs
        cname=paste("contig","type",i,"group",j,"possibility",k,sep="-")
        c=list(chrom1=selrdps$chrom1[1],start1=start1[k],end1=end1[k],
               chrom2=selrdps$chrom2[1],start2=start2[k],end2=end2[k],
               name=cname, score=nrow(selrdps),
               strand1=strand1,strand2=strand2,direction1=direction1[k],direction2=direction2[k],
               type=s, support=nrow(selrdps), cluster=paste("[", paste(unlist(selrdps$name),collapse=","), "]", sep=""))
        contigs[[length(contigs)+1]]=c
        if(is.null(names(contigs))) {
          names(contigs)=cname #first element
        } else {
          names(contigs)[length(contigs)]=cname #continuing elements
        }
      }
    }
  } else { #some molecule support, no clustering needed
    next #let's don't worry about it for now, just use breakpoint as determined from plotMolecule analysis
    #currently use the derived information of plotMolecule
    # if(is.na(leftBreakReg1) | is.na(leftBreakReg2)) next #anything non flushed, we gave up inferring contigure
    # 
    # moles = info[info$supportType==s,]
    # 
    # moleLeftFlush1 = sd(moles$start1) < flushSD #first seg left flush <500
    # moleRightFlush1 = sd(moles$end1) < flushSD #first seg right flush <500
    # moleLeftFlush2 = sd(moles$start2) < flushSD #second seg left flush <500
    # moleRightFlush2 = sd(moles$end2) < flushSD #second seg right flush <500
    # 
    # j=1 #only one group - all molecules overlaps
    # k=1 #only one possibility - use the best one
    # 
    # #deal 16 combinations of flushing
    # 
    # if(moleLeftFlush1 & moleRightFlush1 & !moleLeftFlush2 & !moleRightFlush2 ) {
    #   # a short translocated fragment, derive both start and end contigs
    #   # case1: mol2L - mol1L - mol1R - mol2L : 
    #   #   junction 1.1: 
    #   #   junction 1.2:
    #   # case2: r(mol2R) - mol1L - mol1R - r(mol2L)
    #   #   junction 2.1:
    #   #   junction 2.2:
    #   c=list(chrom1=moles$chrom1[1],start1=min(moles$start2),end1=end1[k],
    #          chrom2=moles$chrom2[1],start2=start2[k],end2=end2[k],
    #          name=paste("contig","type",i,"group",j,"possibility",k,sep="-"), score=nrow(moles),
    #          strand1="*",strand2="*",direction1=direction1[k],direction2=direction2[k],
    #          type=s, support=nrow(selrdps))
    # } else if(!moleLeftFlush1 & !moleRightFlush1 & moleLeftFlush2 & moleRightFlush2 ) {
    #   # a short translocated fragment, derive both start and end contigs
    #   
    # } else if(moleLeftFlush1 & !moleRightFlush1 & moleLeftFlush2 & !moleRightFlush2) { 
    #   # derive one contig
    #   
    # } else if(moleLeftFlush1 & !moleRightFlush1 & !moleLeftFlush2 & moleRightFlush2) { 
    #   # derive one contig
    #   
    # } else if(!moleLeftFlush1 & moleRightFlush1 & moleLeftFlush2 & !moleRightFlush2) {
    #   # derive one contig
    # } else if(!moleLeftFlush1 & moleRightFlush1 & !moleLeftFlush2 & moleRightFlush2) {
    #   # derive one contig
    # } else {
    #   message("Warning: unsupported molecule configuration")
    #   next
    # }
  }
}

cat("-Info: writing contig table...\n")
contigsDF = Reduce(rbind, lapply(contigs,as.data.frame)) #save contigs
contigsDF = cbind(contigsDF, callInfo)
write.table(contigsDF, file=paste(plotPrefix,"rdp","case","ctg",sep="."), sep="\t", quote=F, row.names=F, col.names=F)

cat("-Info: getting contig seqs...\n")
contigSeqs = x <- vector(mode = "list", length = length(contigs))
#tmpSlopSize = min(breakSlopSize,maxSlopSize)
tmpSlopSize = breakSlopSize

formatCoordinate=function(chrom1,start1,end1,chrom2,start2,end2,type,support){
  paste(chrom1,start1,end1,chrom2,start2,end2,type,support,sep="|")
}

formContig=function(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,tag,type,support){
  contigSeqInfo = formatCoordinate(seg1chr, seg1start, seg1end, seg2chr, seg2start, seg2end,type,support)
  contigName = paste(tag,contigSeqInfo,sep="")
  contigSeg1 = refSeq[[seg1chr]][seg1start:seg1end]
  if(seg1rc) contigSeg1 = Biostrings::reverseComplement(contigSeg1)
  contigSeg2 = refSeq[[seg2chr]][seg2start:seg2end]
  if(seg2rc) contigSeg2 = Biostrings::reverseComplement(contigSeg2)
  setNames(list(Biostrings::xscat(contigSeg1, Biostrings::DNAString(x=padding), contigSeg2)), contigName)
}

retrieveContig = function(i, contigs, breakSlopSize, padding, refSeq) {
  contig = contigs[[i]]
  contigSeqs = NULL
  #contigInfo = NULL
  
  #only following 16 patterns are expected
  #(+,R) x (+,R) ; (+,R) x (-,R) ; (+,R) x (+,L) ; (+,R) x (-,L) 
  #(-,R) x (+,R) ; (-,R) x (-,R) ; (-,R) x (+,L) ; (-,R) x (-,L)
  #(+,L) x (+,R) ; (+,L) x (-,R) ; (+,L) x (+,L) ; (+,L) x (-,L) 
  #(-,L) x (+,R) ; (-,L) x (-,R) ; (-,L) x (+,L) ; (-,L) x (-,L)
  
  #since the breakpoint is in the middle of a read pair. i.e. in the middle of a fragment
  #only cases listed below are possible
  #this is not considering palindormes, in which case + and - are unresolvable

  #let's say the two segment reference sequences are
  #     seg1        seg2
  # 5-AC[ACGTGT]GA-AT[ATCACA]GC-3
  # 3-TG[TGCACA]CT-TA[TAGTGT]CG-5
  # [READ]
  # <junction edge>
  
  #(+,R1) x (+,R2) implies 5-AC[ACGTGT]<GA>-<GC>[TGTGAT]AT-3
  #(+,R2) x (+,R1) implies 5-ATATCACAGC-AGTGTGCACA-3
  #(+,R1) x (-,L2) implies 5-AC[ACGTGT]<GA>-<AT>[ATCACA]GC-3
  #(-,L1) x (+,R2) implies 5-ATATCACAGC-ACACGTGTGA-3
  #(-,L1) x (-,L2) implies 5-TCACACGTGT-ATATCACAGC-3
  #(-,L2) x (-,L1) implies 5-GCTGTGATAT-ACACGTGTGA-3
  
  #contig:
  #  [seg1] - read => x <= read - [seg2]; 
  #reads are part of segments
  #read order is important, let's say the order is (r1,r2) and (r2,r1).
  #we assign the order to (-)/(-), (-)/(+) and (+)/(+), 
  #then we have (1-,2-),(2-,1-),(1-,2+),(1+,2-),(1+,2+),(2+,1+) six cases
  #NOTE: for L-R pairs, mapping orders are obvious; for L-L, R-R pairs, 1 and 2 are added to dnote their mapping order
  
  #generate (+,R) x (-,L) contig for compatible modes 
  # [seg1] | R+ => x <= L- | [seg2]
  if(paste(contig$strand1,contig$direction1) %in% c("+ R","+ U","* R","* U")  
     && paste(contig$strand2,contig$direction2) %in% c("- L","- U","* L","* U") ) {
    seg1chr=contig$chrom1; seg1start=contig$start1-breakSlopSize+1; seg1end=contig$end1; seg1rc=FALSE
    seg2chr=contig$chrom2; seg2start=contig$start2; seg2end=contig$end2+breakSlopSize-1; seg2rc=FALSE
    contigSeqs = Biostrings::append(contigSeqs, formContig(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,paste("contig-",i,"(R+,L-)|",sep=""),contig$type,contig$support))
    #contigInfo = append(contigInfo, list(chrom1=contig$chrom1, start1=contig$start1, end1=contig$end1, 
  }
  #generate (-,L) x (+,R) contig for compatible modes
  # [seg1] | R+ => x <= L- | [seg2]
  if(paste(contig$strand1,contig$direction1) %in% c("- L","- U","* L","* U") 
     && paste(contig$strand2,contig$direction2) %in% c("+ R","+ U","* R","* U") ) {
    seg1chr=contig$chrom2; seg1start=contig$start2-breakSlopSize+1; seg1end=contig$end2; seg1rc=FALSE
    seg2chr=contig$chrom1; seg2start=contig$start1; seg2end=contig$end1+breakSlopSize-1; seg2rc=FALSE
    contigSeqs = Biostrings::append(contigSeqs, formContig(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,paste("contig-",i,"(L-,R+)|",sep=""),contig$type,contig$support))
  }
  #generate (+,R1) x (+,R2) contig for compatible modes
  # [seg1] | R1+ => x <= rc(R2+) | [rc(seg2)]
  if(paste(contig$strand1,contig$direction1) %in% c("+ R","+ U","* R","* U")  
     && paste(contig$strand2,contig$direction2) %in% c("+ R","+ R","* R","* U") ) {
    seg1chr=contig$chrom1; seg1start=contig$start1-breakSlopSize+1; seg1end=contig$end1; seg1rc=FALSE
    seg2chr=contig$chrom2; seg2start=contig$start2; seg2end=contig$end2+breakSlopSize-1; seg2rc=TRUE
    contigSeqs = Biostrings::append(contigSeqs, formContig(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,paste("contig-",i,"(R1+,R2+)|",sep=""),contig$type,contig$support))
  }
  #generate (+,R2) x (+,R1) contig for compatible modes
  # [seg1] | R2+ => x <= rc(R1+) | [rc(seg2)]
  if(paste(contig$strand1,contig$direction1) %in% c("+ R","+ U","* R","* U")  
     && paste(contig$strand2,contig$direction2) %in% c("+ R","+ R","* R","* U") ) {
    seg1chr=contig$chrom2; seg1start=contig$start2-breakSlopSize+1; seg1end=contig$end2; seg1rc=FALSE
    seg2chr=contig$chrom1; seg2start=contig$start1; seg2end=contig$end1+breakSlopSize-1; seg2rc=TRUE
    contigSeqs = Biostrings::append(contigSeqs, formContig(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,paste("contig-",i,"(R2+,R1+)|",sep=""),contig$type,contig$support))
  }
  #generate (-,L1) x (-,L2) contig for compatible modes
  # [rc(seg1)] | rc(L1-) => x <= L2- | [seg2]
  if(paste(contig$strand1,contig$direction1) %in% c("- L","- U","* L","* U")  
     && paste(contig$strand2,contig$direction2) %in% c("- L","- U","* L","* U") ) {
    seg1chr=contig$chrom1; seg1start=contig$start1-breakSlopSize+1; seg1end=contig$end1; seg1rc=TRUE
    seg2chr=contig$chrom2; seg2start=contig$start2; seg2end=contig$end2+breakSlopSize-1; seg2rc=FALSE
    contigSeqs = Biostrings::append(contigSeqs, formContig(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,paste("contig-",i,"(L1-,L2-)|",sep=""),contig$type,contig$support))
  }
  #generate (-,L2) x (-,L1) contig for compatible modes
  # [rc(seg1)] | rc(L2-) => x <= L1- | [seg2]
  if(paste(contig$strand1,contig$direction1) %in% c("- L","- U","* L","* U")  
     && paste(contig$strand2,contig$direction2) %in% c("- L","- U","* L","* U") ) {
    seg1chr=contig$chrom2; seg1start=contig$start2-breakSlopSize+1; seg1end=contig$end2; seg1rc=TRUE
    seg2chr=contig$chrom1; seg2start=contig$start1; seg2end=contig$end1+breakSlopSize-1; seg2rc=FALSE
    contigSeqs = Biostrings::append(contigSeqs, formContig(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,paste("contig-",i,"(L2-,L1-)|",sep=""),contig$type,contig$support))
  }
  
  contigSeqs
}

#contigSeqs[[1]] = retrieveContig(1, contigs, breakSlopSize, padding, refSeq)

#need to work on this
contigSeqs = Biostrings::DNAStringSet(unlist(lapply(seq_along(contigs),"retrieveContig", contigs, breakSlopSize, padding, refSeq)))

#output to fasta
message("-Info: write contigs to:", paste(plotPrefix,"rdp","case","fas",sep="."))
Biostrings::writeXStringSet(contigSeqs, file=paste(plotPrefix,"rdp","case","fas",sep="."), format="fasta")

#quit()
#old proof-of-idea scripts
#IMbr1=min(info[info$supportType=="rdpFR",]$start1)
#IMbr2=min(info[info$supportType=="rdpFR",]$start2)
#RRbr1=min(info[info$supportType=="rdpRR",]$start1)
#RRbr2=min(info[info$supportType=="rdpRR",]$start2)
#FFbr1=max(info[info$supportType=="rdpFF",]$start1)
#FFbr2=max(info[info$supportType=="rdpFF",]$start2)

#only trust the start location of mapping
#hassupportType=sapply(supportType, function(x) { any(grepl(x,info$supportType)) } )

cat("PLOTMOLECULEDONE in",proc.time()[3]-ptm[3],"seconds! \n")

#determine breakpoint automatically, we will just use the start of mapping

# a manual procedure example of rdpRR, as in contig2
# diff(info[info$supportType=="rdpRR",]$start1) #find outlier in differ, 1 in this case, exlude it
# diff(info[info$supportType=="rdpRR",]$start1[2:7]) #find the flushing ends, left in this case
# min(info[info$supportType=="rdpRR",]$start1[2:7]) #find the leftmost of flushing ends, 116376679
# posSort=sort(info[info$supportType=="rdpRR",]$start2, index.return=T) #if ends not sorted, sort first 
# excludeIndex=posSort$ix[which(diff(posSort$x)>500)+1] #find outlier (>500) by differ, 1 in this case
# Here we can use the results from start1 ends
# diff(info[info$supportType=="rdpRR",]$start2[-excludeIndex]) #find the flushing ends, left in this case
# min(info[info$supportType=="rdpRR",]$start2[2:7]) #find the leftmost of flushing ends, 131457138
# posSort=sort(info[info$supportType=="rdpFR",]$start1, index.return=T) #if ends not sorted, sort first 
# excludeIndex=posSort$ix[which(diff(posSort$x)>500)+1] #nothing to be excluded
# diff(sort(info[info$supportType=="rdpFR",]$start1))  #find the flushing ends, left in this case
# min(sort(info[info$supportType=="rdpFR",]$start1))

# determine from molecule
# find flushing ends
#   sd(info[info$supportType=="molecule",]$start1) #first seg left flush <500  
#   sd(info[info$supportType=="molecule",]$end1) #first seg right flush <500
#   sd(info[info$supportType=="molecule",]$start2) #second seg left flush <500  
#   sd(info[info$supportType=="molecule",]$end2) #second seg right flush <500
# if one end bloth flush, determine short segment by taking [min(start),max(end)]
#   moleStart2 = min(info[info$supportType=="molecule",]$start2)
#   moleEnd2 = max(info[info$supportType=="molecule",]$end2)
#   #deduce the two breakpoints in long segment by taking the median (median is not a bad guess)
#   moleStart1 = round(median(info[info$supportType=="molecule",]$start1+info[info$supportType=="molecule",]$end1)/2)-1500 #116376483 #KOD-hifi 1-5kb; Taq 500bp; Phusion 20kb
#   moleEnd1 = round(median(info[info$supportType=="molecule",]$start1+info[info$supportType=="molecule",]$end1)/2)+1500
#   in this case, a maxSlopSize is bounded by the size of short segment 
#   maxSlopSize = max(info[info$supportType=="molecule",]$end2-info[info$supportType=="molecule",]$start2)
#all 4 contig combinations: 
#              moleStart1 (R) x moleStart2 (L)
#              moleStart1 (R) x moleEnd2 (R)
#              moleStart2 (L) x moleEnd1 (L)
#              moleEnd2 (R)   x moleEnd1 (L)

# chr2: 116376677 | chr2: 116368067 x chr9: 131457051 | chr 9: 131457138 - chr 9: 131458625 x chr2: 116376628

# chr2: 116376628 x chr9: 131457051 # [contig1: FFbr1 + right x FFbr2 + right] 
# chr2: 116368067 x chr9: 131457138 # [contig2: RRbr1 - left x RRbr2 - right] 
# chr2: 116376677 x chr2: 131457051 # [contig3: IMbr1 * left x IMbr2 * unknown]
#
# chr2: 116368067 x chr9: 131457138 - chr9: 131458625 x chr2: 116376628 [RRbr1 - RRbr2]

#NOTE breakSlopSize needs to be capped;
#direction is where the breakpoint is located relative to the given the coordinate
#L is left, R is right, U is unknown
# contig1 = list(chrom1 = "chr2", 
#                start1 = 116376628, 
#                end1 = 116376628,
#                strand1 = "+",
#                direction1 = "R",
#                chrom2 = "chr9",
#                start2 = 131458625,
#                end2 = 131458625, 
#                strand2 = "+",
#                direction2 = "R")
# contig2 = list(chrom1 = "chr2", 
#                start1 = 116368067,
#                end1 = 116368067,
#                strand1 = "-",
#                direction1 = "L",
#                chrom2 = "chr9",
#                start2 = 131458503,
#                end2 = 131458503,
#                strand2 = "-",
#                direction2 = "L"
#                )
# contig3 = list(chrom1 = "chr2", 
#                start1 = 116376677,
#                end1 = 116376677,
#                strand1 = "*",
#                direction1 = "L",
#                chrom2 = "chr9",
#                start2 = 131457051,
#                end2 = 131457051,
#                strand2 = "*",
#                direction2 = "U"
# )
#contigs = list(contig1,contig2,contig3)

#print(contigs)


#quick calculation of these numbers
#caseFeatureMol=31 #cat refined.reg1.mol.case refined.reg2.mol.case | cut -f11 | sort | uniq -c | grep 2 | wc -l
#caseTotalMol=37+44-31=50 #`cat refined.reg1.mol.case | wc -l` + `cat refined.reg2.mol.case | wc -l` - caseFeatureMol
#controlFeatureMol=0 #cat refined.reg1.mol.control refined.reg2.mol.control | cut -f11 | sort | uniq -c | grep 2 | wc -l
#controlTotalMol=12+7-0=19 #`cat refined.reg1.mol.control | wc -l` + `cat refined.reg2.mol.control | wc -l` - controlFeatureMol
#calcPosteriorForProportion(caseFeatureMol, caseTotalMol, controlTotalMol-controlFeatureMol+1,controlFeatureMol+1) = calcPosteriorForProportion(31, 50, 20, 1)
#d <- data.frame(s = 31, k = 50)
# pd <- tibble(
#   x = seq(0, 1, by = .01),
#   Prior = dbeta(x, alpha1, beta1)
# )
# ggplot(pd, aes(x, Prior)) +
#   geom_line() +
#   coord_cartesian(xlim = 0:1, ylim = c(0, 6), expand = 0) +
#   labs(y = "Density", x = bquote(theta))
# 
# pd$Posterior <- dbeta(pd$x, 31+alpha1, (50-31)+beta1)
# pdw <- gather(pd, key=Type, value=density, Prior:Posterior)
#   
# ggplot(pdw, aes(x, density, col=Type)) +
#   geom_line() +
#   annotate("point", x=c(.5, .5), y = c(pdw$density[pdw$x==.5])) +
#   annotate("label", x=c(.5, .5), 
#   y = pdw$density[pdw$x==.5], 
#   label = round(pdw$density[pdw$x==.5], 3),
#   vjust=-.5)
# 
# filter(pd, x >=objectiveFrequency) %>% 
#   mutate(BF01 = Posterior/Prior,
#          BF10 = 1/BF01) %>% 
#   round(3)

#do the bayesian statistics here:
#objectiveFrequency = 0.05 # H0[0,objectiveFrequency] vs H1[objectiveFrequency,1]
#bayesPlot =  paste(plotPrefix,"bayes","pdf",sep=".")
#adapted from https://mvuorre.github.io/post/2017/bayes-factors-with-brms/
# if(length(suffixes)==2) { #second row is control, i.e. prior
#   alpha0=1; beta0=1; # pseudo count, #of success, #of failures
#   alpha1=alpha0+infoTabNew[2,]$featureMol #pseudo count adjusted prior # of success
#   beta1=beta0+infoTabNew[2,]$totalMol-2*infoTabNew[2,]$featureMol #psuedo count adjusted prior # of failure; two times bc molecule counted twice
#   d=data.frame(s = infoTabNew[1,]$featureMol, k = infoTabNew[1,]$totalMol - infoTabNew[1,]$featureMol) # s success of k tries
#   Prior=brms::set_prior(paste("beta(",alpha1,",", beta1,")",sep=""), class = "b", lb = 0, ub = 1)
#   m=brms::brm(s | trials(k) ~ 0 + intercept, family = binomial(link="identity"), prior = Prior, data = d, sample_prior = TRUE, iter = 1e4, cores = 1)
#   samples=brms::posterior_samples(m, "b")
#   h=brms::hypothesis(m, paste("intercept >", minAlleleFraction)) #specifies H0
#   BF01=h$hypothesis$Evid.Ratio #BF01, additional evidence of H0 vs H1
#   Star=h$hypothesis$Star
#   #Sample is 18
#   infoTabNew$minAlleleFraction=minAlleleFraction #19
#   infoTabNew$controlTotalMol=infoTabNew[2,]$totalMol #20
#   infoTabNew$controlFeatureMol=infoTabNew[2,]$featureMol #21
#   infoTabNew$BF01=BF01 #22
#   infoTabNew$Star=Star #23
#   pdf(file=bayesPlot)
#   plot(h)
#   dev.off()
# } else {
#   infoTabNew$minAlleleFraction=minAlleleFraction
#   infoTabNew$BF01=NA
#   infoTabNew$Star=NA
#   infoTabNew$controlTotalMol=NA
#   infoTabNew$controlFeatureMol=NA
# }

#infoTabNew=infoTabNew[c(names(infoTabNew)[1:17],names(infoTabNew)[22:23],names(infoTabNew)[18:21])]
#infoTabOutput=infoTabNew[,c(1:17,20:21,19,22:23,18)]
