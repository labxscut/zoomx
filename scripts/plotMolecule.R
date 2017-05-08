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

#for(p in c("optparse", "rjson", "ggplot2", "RColorBrewer")) myrequire(p)
#for(p in c("GenomicRanges", "ggbio")) myrequire(p,repo="Bioc")
for(p in c("optparse", "rjson", "ggplot2", "RColorBrewer","GenomicRanges", "ggbio")) require(p,character.only=T)

#genome plot packages: ggbio, Sushi, genoplotR, circlize 

####Options####
option_list2 <- list(
  make_option(c("-f", "--flushSD"), type="integer", default=5000,
              help="maximum of flushing sd [default %default]"),
  make_option(c("-p", "--gapSize"), type="integer", default=100000,
              help="2 times the boundary padding bp [default %default]"),
  make_option(c("-s", "--scaleSize"), type="integer", default=10000,
              help="boundary rounded to bp [default %default]"),
  make_option(c("-b", "--breakCI"), type="numeric", default=0.95,
              help="limit of CI to estimate break point [default %default]"),
  make_option(c("-u", "--isUCSC"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)

parser <- OptionParser(usage = "%prog [options] plotDir statsJSON", option_list=option_list2)
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)
cat("-Info: invoking command:",commandArgs(),"\n")
breakCI=cmd$options$breakCI
verbose=cmd$options$noQuiet
debug=cmd$options$debug
isUCSC=cmd$options$isUCSC
gapSize=cmd$options$gapSize
scaleSize=cmd$options$scaleSize
flushSD=cmd$options$flushSD
if(length(cmd$args)!=2 & !debug){ print_help(parser); quit() }
plotDir=cmd$args[1]
statsFile=cmd$args[2]

####Debug####
if(debug){ 
  #options(warn=2)
  isUCSC=FALSE     #input has no chr in chromosome name
  #isUCSC=TRUE
  breakCI=0.95
  gapSize = 100000     #boundary padding is 10K/2 
  scaleSize = 10000    #boundary rounded to 10K
  flushSD = 5000
  plotDir=paste("~/ws/crcbrain/02_ZoomX","P05378_7175_mt.wg12.bam.dist.E25_P05378_7176_nt.wg12.bam.dist.E10.somatic.curate.Plot_1",sep="/")
  statsFile=paste("~/ws/crcbrain/02_ZoomX","P05378_7175_mt.wg12.bam.gridScan.cfg.json",sep="/") 
} else { 
  options(warn=-1)
}

####Main####

###LoadInput###
options(stringsAsFactors = FALSE)
statsTab=fromJSON(paste(readLines(statsFile), collapse=""))
setwd(plotDir)
regCols=c("chr1","start1","end1","chr2","start2","end2","molcnt","IM","FF","RR") # in refined
bedpeCols=c("chr1","st1","ed1","chr2","st2","ed2","hp")
bedpeClass=c("character","numeric","numeric","character","numeric","numeric","numeric")
molCols=c("chr1","st1","ed1","barID","molSize","molLen","barMolNum","molID","molHP","molPS","molBX")
molColsOld=c("chr1","st1","ed1","barID","molSize","molLen","barMolNum","molID","molHP","molPS")
refined=tryCatch(read.table("refined.bedpe",sep="\t",col.names=regCols),error=function(e) {NULL})
refinedIM=tryCatch(read.table("refined.IM.bedpe",col.names=bedpeCols, colClasses=bedpeClass),error=function(e) {NULL})
refinedFF=tryCatch(read.table("refined.FF.bedpe",col.names=bedpeCols, colClasses=bedpeClass),error=function(e) {NULL})
refinedRR=tryCatch(read.table("refined.RR.bedpe",col.names=bedpeCols, colClasses=bedpeClass),error=function(e) {NULL})
refinedReg1=tryCatch(read.table("refined.reg1.mol",col.names=molCols),error=function(e) {t=read.table("refined.reg1.mol",col.names=molColsOld);t[["molBx"]]=NA;t}, finally={NULL})
refinedReg2=tryCatch(read.table("refined.reg2.mol",col.names=molCols),error=function(e) {t=read.table("refined.reg2.mol",col.names=molColsOld);t[["molBx"]]=NA;t}, finally={NULL})
infoTab=refined

sharedBarID = unique(refinedReg1$barID[refinedReg1$barID %in% refinedReg2$barID])
mergedMolList = lapply(sharedBarID, function(b) { sel1=which(refinedReg1$barID==b)
                                  sel2=which(refinedReg2$barID==b)
                                  mdat=cbind(refinedReg1[sel1[1],c(1,2,3)], refinedReg2[sel2[1],c(1,2,3)])
                                  colnames(mdat)=bedpeCols[1:6]
                                  mdat$id=paste("molID",b,sep=":")
                                  mdat$type="molecule"
                                  mdat$barcode=refinedReg1[sel1[1],dim(refinedReg1)[2]]
                                  mdat
                                })
mergedMol = do.call("rbind", mergedMolList)
  
chr1 = unique(mergedMol$chr1)[1]
chr2 = unique(mergedMol$chr2)[1]
refinedIMid=fromJSON(refined$IM)
refinedFFid=fromJSON(refined$FF)
refinedRRid=fromJSON(refined$RR)

#following assume for mergedMol if (chr1==chr2) => st1 <= st2 
mergedIM = do.call("rbind", 
                   apply(refinedIM[,1:6], 1, FUN = function(r, chr1) { 
                     if( r[1]!=chr1 | (r[4]==chr1 & as.numeric(r[2])>as.numeric(r[5])) ) 
                                                                            { t=r[1:3]; r[1:3]=r[4:6]; r[4:6]=t; }
                                                                        as.data.frame(t(r)) 
                                                                     }, chr1) )
if(!is.null(mergedIM)) { #set_colClass(mergedIM, c("character","integer","integer","character","integer","integer"));
                         mergedIM$id = paste("rdpIM",refinedIMid,sep=":"); mergedIM$type="rdpIM"; mergedIM$barcode=refinedIM$hp}

mergedFF = do.call("rbind", 
                   apply(refinedFF[,1:6], 1, FUN = function(r, chr1) { 
                     if( r[1]!=chr1 | (r[4]==chr1 & as.numeric(r[2])>as.numeric(r[5])) ) { t=r[1:3]; r[1:3]=r[4:6]; r[4:6]=t; }
                                                                                  as.data.frame(t(r)) }, chr1))
if(!is.null(mergedFF)) { #set_colClass(mergedFF, c("character","integer","integer","character","integer","integer")); 
                         mergedFF$id = paste("rdpFF",refinedFFid,sep=":"); mergedFF$type="rdpFF"; mergedFF$barcode=refinedFF$hp} 

mergedRR = do.call("rbind", 
                   apply(refinedRR[,1:6], 1, FUN = function(r, chr1) { 
                     if( r[1]!=chr1 | (r[4]==chr1 & as.numeric(r[2])>as.numeric(r[5])) ) { t=r[1:3]; r[1:3]=r[4:6]; r[4:6]=t; } 
                   as.data.frame(t(r)) }, chr1))
if(!is.null(mergedRR)) { #set_colClass(mergedRR, c("character","integer","integer","character","integer","integer")); 
                         mergedRR$id = paste("rdpRR",refinedRRid,sep=":"); mergedRR$type="rdpRR"; mergedRR$barcode=refinedRR$hp } 

if(!isUCSC & nrow(refined)>0) {refined$chr1=paste("chr",refined$chr1,sep="")}
if(!isUCSC & nrow(refined)>0) {refined$chr2=paste("chr",refined$chr2,sep="")}
if(!isUCSC & nrow(refinedReg1)>0) {refinedReg1$chr1=paste("chr",refinedReg1$chr1,sep="")}
if(!isUCSC & nrow(refinedReg2)>0) {refinedReg2$chr1=paste("chr",refinedReg2$chr1,sep="")}
if(!isUCSC & length(mergedMol)>0) {mergedMol$chr1=paste("chr",mergedMol$chr1,sep="")
                               mergedMol$chr2=paste("chr",mergedMol$chr2,sep="")}
if(!isUCSC & length(mergedIM)>0) {mergedIM$chr1=paste("chr",mergedIM$chr1,sep="")
                                 mergedIM$chr2=paste("chr",mergedIM$chr2,sep="")}
if(!isUCSC & length(mergedFF)>0) {mergedFF$chr1=paste("chr",mergedFF$chr1,sep="")
                                 mergedFF$chr2=paste("chr",mergedFF$chr2,sep="")}
if(!isUCSC & length(mergedRR)>0) {mergedRR$chr1=paste("chr",mergedRR$chr1,sep="")
                                 mergedRR$chr2=paste("chr",mergedRR$chr2,sep="")}

Type=c("molecule","rdpFF","rdpRR","rdpIM")
strandInfo=c(".","1","-1",".")
strandDf=c("1","-1",".")
strandGR=c("+","-","*")
supportType=c("molecule","readPair","readPair","readPair")
supportGR=c("molecule","readPair")
supportColors=brewer.pal(3,"Set1")[1:length(supportGR)]
names(supportColors)=supportGR
strandColors=brewer.pal(3,"Set2")
names(strandColors)=strandGR

mergedAll=rbind(mergedMol,mergedIM,mergedFF,mergedRR)
mergedAll$strand1=strandInfo[match(mergedAll$type, Type)]
mergedAll$strand2=strandInfo[match(mergedAll$type, Type)]
mergedAll$supportType=mergedAll$type
colnames(mergedAll)=c("chrom1","start1","end1","chrom2","start2","end2","name","score","barcode","strand1","strand2","supportType")
mergedAll$score=0

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

convertedAll = mergedAll
convertedAll = set_colClass(convertedAll, c("character","numeric","numeric","character","numeric","numeric","character","numeric","character","character","character","character"))
convertedAll$chrom1 = "chrA"
convertedAll$start1 = convertedAll$start1 - scaledReg1St + 1
convertedAll$end1 = convertedAll$end1 - scaledReg1St + 150
convertedAll$chrom2 = "chrA"
convertedAll$start2 = convertedAll$start2 - scaledReg2St + shiftReg2St + 1
convertedAll$end2 = convertedAll$end2 - scaledReg2St + shiftReg2St +150
convertedAll$supportType = supportType[match(convertedAll$supportType,Type)]
#convertedAll$supportType = factor(convertedAll$supportType, levels=supportGR)
convertedAll$strand1 = strandGR[match(convertedAll$strand1,strandDf)]
convertedAll$strand2 = strandGR[match(convertedAll$strand2,strandDf)]
#convertedAll$strand1 = factor(strandGR[match(convertedAll$strand1,strandDf)], levels=strandGR)
#convertedAll$strand2 = factor(strandGR[match(convertedAll$strand2,strandDf)], levels=strandGR)
flushSDReg1 = c(sd(convertedAll$start1),sd(convertedAll$end1))
flushSDReg2 = c(sd(convertedAll$start2),sd(convertedAll$end2))
leftBreakReg1 = if(min(flushSDReg1)<flushSD) TRUE else NA
leftBreakReg2 = if(min(flushSDReg2)<flushSD) TRUE else NA
if(!is.na(leftBreakReg1)) leftBreakReg1 = flushSDReg1[1] < flushSDReg1[2]
if(!is.na(leftBreakReg2)) leftBreakReg2 = flushSDReg2[1] < flushSDReg2[2]
# leftBreak: TRUE iff left edge of Reg1 is break, FALSE if right edge, NA if break is middle
breakReg1 = round(if(is.na(leftBreakReg1)) { median(c(convertedAll$start1,convertedAll$end1)) } else if (leftBreakReg1) { quantile(convertedAll$start1,1-breakCI) } else { quantile(convertedAll$end1,breakCI) })
breakReg1CILow = min(convertedAll$start1)
breakReg1CIHigh = max(convertedAll$end1)
breakReg2 = round(if(is.na(leftBreakReg2)) { median(c(convertedAll$start2,convertedAll$end2)) } else if (leftBreakReg2) { quantile(convertedAll$start2,1-breakCI) } else { quantile(convertedAll$end2,breakCI) })
breakReg2CILow = min(convertedAll$start2)
breakReg2CIHigh = max(convertedAll$end2)

cat("gr...\n")
#print(convertedAll[convertedAll$end1>convertedAll$start2,c(1,2,3,4,5,6)])
#print(mergedAll[convertedAll$end1>convertedAll$start2,c(1,2,3,4,5,6)])
gr <- GRanges(seqnames = c(convertedAll$chrom1, convertedAll$chrom2), 
              ranges = IRanges(start = c(convertedAll$start1, convertedAll$start2), 
                               end = c(convertedAll$end1, convertedAll$end2)),
              #strand = factor(c(as.character(convertedAll$strand1),as.character(convertedAll$strand2)), levels=strandGR),
              strand = c(convertedAll$strand1,convertedAll$strand2),
              value = rep(0,2*nrow(convertedAll)), 
              score = c(convertedAll$score, convertedAll$score), 
              supportType = c(convertedAll$supportType, convertedAll$supportType), 
              pair = c(convertedAll$name, convertedAll$name))
cat("grComp...\n")
print(convertedAll)
grComp <- GRanges(seqnames = convertedAll$chrom1,
                  ranges = IRanges(start = convertedAll$end1, end=convertedAll$start2),
                  strand = rep(NA,nrow(convertedAll)),
                  value = rep(0,nrow(convertedAll)),
                  score = rep(0,nrow(convertedAll)),
                  supportType = convertedAll$supportType,
                  pair = convertedAll$name)
tickReg1 = c(shiftReg12[1]+gapSize/2, shiftReg12[2]-gapSize/2, breakReg1) 
tickReg2 = c(shiftReg12[3]+gapSize/2, shiftReg12[4]-gapSize/2, breakReg2)
tickLab1 = c(boundaryReg12[1]+gapSize/2, boundaryReg12[2]-gapSize/2, boundaryReg12[1]+breakReg1) 
tickLab2 = c(boundaryReg12[3]+gapSize/2, boundaryReg12[4]-gapSize/2, boundaryReg12[3]+breakReg2-shiftReg12[2])
tickLab1 = paste(c(chr1,"",""),tickLab1,sep=":")
tickLab2 = paste(c(chr2,"",""),tickLab2,sep=":")

####WriteOutput####
plotPrefix=paste(refined[1:6],collapse="_")
mergedAllFile = paste(plotPrefix,"mol",sep=".")
write.table(mergedAll,file=mergedAllFile,col.names=TRUE,row.names=FALSE,sep="\t",quote=FALSE)

infoFile = paste(plotPrefix,"info","txt",sep=".")
infoTab = cbind(infoTab[,seq(7)], as.data.frame(list(breakpoint1=boundaryReg12[1]+unname(breakReg1[1]),
                                                      breakpoint2=boundaryReg12[3]+unname(breakReg2[1])-shiftReg12[2],
                                                      SVType="UNK",
                                                      AlleleFraction=infoTab$molcnt/statsTab$Cem,
                                                      ReadPair=length(refinedIMid)+length(refinedFFid)+length(refinedRRid))))
write.table(infoTab,file=infoFile,col.names=FALSE,row.names=FALSE,sep="\t",quote=FALSE)
#chr1  start1	end1	chr2	start2	end2	molcnt	breakpoint1	breakpoint2	SVType	AlleleFraction	ReadPair

jsonFile=file(paste(plotPrefix,"json",sep="."))
writeLines(toJSON(list(reg1=chr1,
                       break1=infoTab$breakpoint1,
                       break1CILow=boundaryReg12[1]+unname(breakReg1CILow[1]),
                       break1CIHigh=boundaryReg12[1]+unname(breakReg1CIHigh[1]),
                       reg2=chr2,
                       break2=infoTab$breakpoint2, 
                       break2CILow=boundaryReg12[3]+unname(breakReg2CILow[1])-shiftReg12[2],
                       break2CIHigh=boundaryReg12[3]+unname(breakReg2CIHigh[1])-shiftReg12[2],
                       sharedBarcode=mergedMol$barcode)),
                       con=jsonFile)
close(jsonFile)

cat("plotFile=",paste(plotDir, paste(plotPrefix,"pdf",sep="."),sep="/"),"\n")
pdf(file=paste(plotPrefix,"pdf",sep="."), width=16, height=12)
p = ggplot()+geom_arrowrect(gr, aes(color=supportType, fill=strand, group=pair))
p = p + geom_segment(grComp, aes(color=supportType, group=pair))
p = p + scale_color_manual(values=supportColors, na.value="black")
p = p + scale_fill_manual(values=strandColors, na.value="black")
p = p + scale_x_continuous("Genomic Region Pair", breaks=tickReg1, labels=tickLab1,
                           limits=c(shiftReg12[1],shiftReg12[4]),
                           sec.axis=sec_axis(~., breaks=tickReg2, labels=tickLab2))
p = p + geom_vline(xintercept=breakReg1) + geom_vline(xintercept=breakReg2)
print(p)
dev.off()
cat("plotMoleculeDone in",proc.time()[3]-ptm[3],"seconds! \n")

