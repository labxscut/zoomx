#!/usr/bin/env Rscript
version="REPLACE_WITH_COMMIT_OR_VERSION"
#library("BSgenome.Hsapiens.UCSC.hg19")
#genome = BSgenome.Hsapiens.UCSC.hg19
#seqlengths(genome)
#library("Sushi")
library("optparse")
library("rjson")
library("ggbio")
library("ggplot2")
library("GenomicRanges")
library("RColorBrewer")

####Externalities####
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

#genome plot packages: ggbio, Sushi, genoplotR, circlize 

option_list2 <- list(
  make_option(c("-f", "--flushSD"), type="integer", default=5000,
              help="maximum of flushing sd [default %default]"),
  make_option(c("-p", "--gapSize"), type="integer", default=100000,
              help="2 times the boundary padding bp [default %default]"),
  make_option(c("-s", "--scaleSize"), type="integer", default=10000,
              help="boundary rounded to bp [default %default]"),
  make_option(c("-u", "--isUCSC"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)

parser <- OptionParser(usage = "%prog [options] plotDir", option_list=option_list2)
cat("-Info: zoomx vesion:",version,"\n")
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)
verbose=cmd$options$noQuiet
debug=cmd$options$debug
isUCSC=cmd$options$isUCSC
gapSize=cmd$options$gapSize
scaleSize=cmd$options$scaleSize
flushSD=cmd$options$flushSD
ptm=proc.time() #time keeper
if(length(cmd$args)!=1 & !debug){ print_help(parser); quit() }
plotDir=cmd$args[1]

if(debug){ 
  #options(warn=2)
  #isUCSC=FALSE     #input has no chr in chromosome name
  isUCSC=TRUE
  gapSize = 100000     #boundary padding is 10K/2 
  scaleSize = 10000    #boundary rounded to 10K
  flushSD = 5000
  #plotDir=paste("~/ws/10x6","testPlot-1",sep="/")
  #plotDir=paste("~/ws/10x6","plots/NA12878_WGS_phased_possorted_bam.bam.dis.E25.refined.Plot-2",sep="/")
  #plotDir=paste("~/ws/10x6","P00526_2721_mt.st.10X_v2.bam.dis.E25.refined.somatic.Plot-4",sep="/")
  plotDir=paste("~/ws/10x6","P00526_2721_mt.st.10X_v2.bam.dis.E25.refined.somatic.Plot-2",sep="/")
} else { 
  options(warn=-1)
}
options(stringsAsFactors = FALSE)
setwd(plotDir)
regCols=c("chr1","start1","end1","chr2","start2","end2","molcnt","IM","FF","RR") # in refined
bedpeCols=c("chr1","st1","ed1","chr2","st2","ed2","hp")
bedpeClass=c("character","numeric","numeric","character","numeric","numeric","numeric")
molCols=c("chr1","st1","ed1","barID","molSize","molLen","barMolNum","molID","molHP","molPS")
refined=tryCatch(read.table("refined.bedpe",sep="\t",col.names=regCols),error=function(e) {NULL})
refinedIM=tryCatch(read.table("refined.IM.bedpe",col.names=bedpeCols, colClasses=bedpeClass),error=function(e) {NULL})
refinedFF=tryCatch(read.table("refined.FF.bedpe",col.names=bedpeCols, colClasses=bedpeClass),error=function(e) {NULL})
refinedRR=tryCatch(read.table("refined.RR.bedpe",col.names=bedpeCols, colClasses=bedpeClass),error=function(e) {NULL})
refinedReg1=tryCatch(read.table("refined.reg1.mol",col.names=molCols),error=function(e) {NULL})
refinedReg2=tryCatch(read.table("refined.reg2.mol",col.names=molCols),error=function(e) {NULL})
if(!isUCSC & nrow(refined)>0) {refined$chr1=paste("chr",refined$chr1,sep="")
             refined$chr2=paste("chr",refined$chr2,sep="")}
if(!isUCSC & nrow(refinedIM)>0) {refinedIM$chr1=paste("chr",refinedIM$chr1,sep="")
             refinedIM$chr2=paste("chr",refinedIM$chr2,sep="")}
if(!isUCSC & nrow(refinedFF)>0) {refinedFF$chr1=paste("chr",refinedFF$chr1,sep="")
             refinedFF$chr2=paste("chr",refinedFF$chr2,sep="")}
if(!isUCSC & nrow(refinedRR)>0) {refinedRR$chr1=paste("chr",refinedRR$chr1,sep="")
             refinedRR$chr2=paste("chr",refinedRR$chr2,sep="")}
if(!isUCSC & nrow(refinedReg1)>0) {refinedReg1$chr1=paste("chr",refinedReg1$chr1,sep="")}
if(!isUCSC & nrow(refinedReg2)>0) {refinedReg2$chr1=paste("chr",refinedReg2$chr1,sep="")}

sharedBarID = unique(refinedReg1$barID[refinedReg1$barID %in% refinedReg2$barID])
mergedMolList = lapply(sharedBarID, function(b) { sel1=which(refinedReg1$barID==b)
                                  sel2=which(refinedReg2$barID==b)
                                  mdat=cbind(refinedReg1[sel1[1],c(1,2,3)], refinedReg2[sel2[1],c(1,2,3)])
                                  colnames(mdat)=bedpeCols[1:6]
                                  mdat$id=paste("molID",b,sep=":")
                                  mdat$type="molecule"
                                  mdat
                                })
mergedMol = do.call("rbind", mergedMolList)
chr1 = unique(mergedMol$chr1)[1]
chr2 = unique(mergedMol$chr2)[1]
refinedIMid=fromJSON(refined$IM)
refinedFFid=fromJSON(refined$FF)
refinedRRid=fromJSON(refined$RR)

mergedIM = do.call("rbind", 
                   apply(refinedIM[,1:6], 1, FUN = function(r, chr1) { 
                                                                        if( r[1]!=chr1 | (r[1]==chr1 & as.numeric(r[2])>as.numeric(r[5])) ) 
                                                                            { t=r[1:3]; r[1:3]=r[4:6]; r[4:6]=t; }
                                                                        as.data.frame(t(r)) 
                                                                     }, chr1) )
if(!is.null(mergedIM)) { #set_colClass(mergedIM, c("character","integer","integer","character","integer","integer"));
                         mergedIM$id = paste("rdpIM",refinedIMid,sep=":"); mergedIM$type="rdpIM"}

mergedFF = do.call("rbind", 
                   apply(refinedFF[,1:6], 1, FUN = function(r, chr1) { 
                     if(r[1]!=chr1 | (r[1]==chr1 & as.numeric(r[2])>as.numeric(r[5])) ) { t=r[1:3]; r[1:3]=r[4:6]; r[4:6]=t; }
                                                                                  as.data.frame(t(r)) }, chr1))
if(!is.null(mergedFF)) { #set_colClass(mergedFF, c("character","integer","integer","character","integer","integer")); 
                         mergedFF$id = paste("rdpFF",refinedFFid,sep=":"); mergedFF$type="rdpFF"} 

mergedRR = do.call("rbind", 
                   apply(refinedRR[,1:6], 1, FUN = function(r, chr1) { 
                     if(r[1]!=chr1 | (r[1]==chr1 & as.numeric(r[2])>as.numeric(r[5])) ) { t=r[1:3]; r[1:3]=r[4:6]; r[4:6]=t; } 
                   as.data.frame(t(r)) }, chr1))
if(!is.null(mergedRR)) { #set_colClass(mergedRR, c("character","integer","integer","character","integer","integer")); 
                         mergedRR$id = paste("rdpRR",refinedRRid,sep=":"); mergedRR$type="rdpRR" } 

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
colnames(mergedAll)=c("chrom1","start1","end1","chrom2","start2","end2","name","score","strand1","strand2","supportType")
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
convertedAll = set_colClass(convertedAll, c("character","numeric","numeric","character","numeric","numeric","character","numeric","character","character","character"))
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
breakCI=0.95
breakReg1 = round(if(is.na(leftBreakReg1)) { median(c(convertedAll$start1,convertedAll$end1)) } else if (leftBreakReg1) { quantile(convertedAll$start1,1-breakCI) } else { quantile(convertedAll$end1,breakCI) })
breakReg2 = round(if(is.na(leftBreakReg2)) { median(c(convertedAll$start2,convertedAll$end2)) } else if (leftBreakReg2) { quantile(convertedAll$start2,1-breakCI) } else { quantile(convertedAll$end2,breakCI) })

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
tickLab1 = paste(c(chr1,chr1,chr1),tickLab1,sep=":")
tickLab2 = paste(c(chr2,chr2,chr2),tickLab2,sep=":")
                    
plotFile=paste(refined[1:6],collapse="-")
cat("plotFile=",paste(plotDir, paste(plotFile,"pdf",sep="."),sep="/"),"\n")
pdf(file=paste(plotFile,"pdf",sep="."), width=16, height=12)
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

#p = p + geom_rect(grComp, aes(color=supportType,group=pair))
#p = p + geom_chevron(grComp, aes(color=supportType,group=pair))

#print(p)
#
#print(p)
#grA=c(gr,grComp)
#p = ggplot()+geom_arrowrect(grA, aes(color=supportType,fill=supportType,group=pair))
#print(p)
#
#tmp=convertedAll[convertedAll$supportType=="molecule",]
#plotBedpe(tmp,chrom="chrA",chromstart=shiftReg12[1],chromend=shiftReg12[4],plottype="lines",
#          color=c("black"), lwd=1)
#tmp=convertedAll[convertedAll$supportType!="molecule",]
#plotBedpe(tmp,chrom="chrA",chromstart=shiftReg12[1],chromend=shiftReg12[4],plottype="loops", flip=TRUE, heights=rnorm(nrow(tmp)),
#          color=c("red","green","blue"), lwd=1)
#dev.off()
#
#tmp1=convertedAll[,c(1,2,3,7,8,9,11)]; names(tmp1)=c("chrom","start","stop","gene","score","strand","supportType")
#tmp2=convertedAll[,c(4,5,6,7,8,10,11)]; names(tmp2)=c("chrom","start","stop","gene","score","strand","supportType")
#convertedAll2=rbind(tmp1,tmp2)
#tmp=convertedAll2; tmp$type="exon"
#pdf(file="refined.pdf")
#plotGenes(tmp,chrom="chrA",chromstart=shiftReg12[1],chromend=shiftReg12[4],type=tmp$type,plotgenetype="box",
#          col=as.factor(tmp$supportType),arrowlength=0)
#dev.off()

# #swith if region 1 and 2 not in order
# #switchReg(mergedMol, chr1)
# #refinedMol1 = refinedReg1[refinedReg1$barID %in% sharedBarID,]  #59 shared barcode in reg1
# #refinedMol2 = refinedReg2[refinedReg2$barID %in% sharedBarID,]  #59 shared barcode in reg2
# convertedMol1 = refinedMol1
# 
# convertedMol2 = refinedMol2
#  
# convertedMol = rbind(convertedMol1[,c(1,2,3,4)],convertedMol2[,c(1,2,3,4)])
# convertedMol[,5] = 1               #score for types: 1 for mol
# convertedMol[,6] = 1                #
# convertedMol[,7] = "exon"
# colnames(convertedMol) = c("chrom","start","stop","gene","score","strand","type")
# convertedMol$gene = paste("molID",convertedMol$gene,sep=":")
# sharedBarID = paste("molID",sharedBarID,sep="-")
# 
# pdf(file="refined.pdf")
# tmp=convertedMol
# plotGenes(tmp,chrom="chrA",chromstart=1,chromend=shiftReg2Ed, 
#           colorby=tmp$score, colorbycol= SushiColors(3))
# #plotGenes(tmp,chrom="chr1",chromstart=1,chromend=shiftReg2Ed,
# #          types=tmp$type, colorby=tmp$score,
# #          colorbycol= SushiColors(5),colorbyrange=c(0,1.0),
# #          labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box")
# dev.off()
# pdf(file="refined.pdf")
# tmp=Sushi_transcripts.bed
# plotGenes(tmp,chrom="chr15",chromstart=72798936,chromend=73100889)
# #plotGenes(tmp,chrom="chr15",chromstart=72798936,chromend=73100889,
#           #types=tmp$type, colorby=tmp$score,
#           #colorbycol= SushiColors(5),colorbyrange=c(0,1.0),
#           #labeltext=TRUE,maxrows=50,height=0.4,plotgenetype="box")
# dev.off()
# #labelgenome(chrom=unique(refined$chr1), chromstart=1, chromend=scaledReg1Ed-shiftReg1,n=3,scale="Kb")
# #labelgenome(chrom=unique(refined$chr2), chromstart=scaledReg2St-shiftReg2, 
# #            chromend=scaledReg2Ed-shiftReg2, n=3, scale="Kb")
# 
# title(main="Junction", #col.main="red", 
#       #sub="Molecule", col.sub="blue", 
#       xlab="Genomic Region Pair", ylab="Read/Molecule Support",
#       #col.lab="green", 
#       cex.lab=0.75)
# axis(1, c(1, scaledReg1Ed-shiftReg1), 
#      paste(unique(refined$chr1), c(scaledReg1St, scaledReg1Ed), sep=":") )
# axis(1, c(scaledReg2St-shiftReg2, scaledReg2Ed-shiftReg2), 
#      paste(unique(refined$chr2), c(scaledReg2St, scaledReg2Ed), sep=":") )
# 
# 
# 
# ###
# 
# # add in IM, FF, RR
# 
# 
# refinedRegLen1 = refined$end1 - refined$start1
# refinedRegLen2 = refined$end2 - refined$start2
# refinedMolLen1 = max(refinedMol1Bound)-min(refinedMol1Bound)
# refinedMolLen2 = max(refinedMol2Bound)-min(refinedMol2Bound)
# 
# 
# refinedMol1[with(refinedMol1, order(st1)),]
# refinedMol2[with(refinedMol2, order(st1)),]
# 
# merged = 
# # make an artificial scaffold
# # gapSize = 100K
# # Len1 = End1 - Start1
# # Len2 = End2 - Start2
# # ChrA = [ 0 ~ Len1 ] + gapSize + [ Len1+gapSize+1 ~ Len1+gapSize+Len2 ] 
# # 
# 
# plotMol=function(chrom, chromstart, chromend, beddata, sharedBarID){
#   beddata$row=as.integer(as.factor(sharedBarID))
#   beddata$color=maptocolors(beddata$row,col=SushiColors(6))
#   plotBed(beddata = beddata,
#           chrom = chrom, 
#           chromstart = chromstart, 
#           chromend = chromend, 
#           type = "circles",
#           color = beddata$color,
#           rownumber = beddata$row,
#           row = "given",
#           plotbg = "grey95",
#           rowlabels = unique(beddata$barID),
#           rowlabelcol = unique(beddata$color) )
#   labelgenome(chrom,chromstart,chromend,scale="Mb")
# }
# 
# plotMol(chrom=refined$chr1,chromstart=refined$start1,chromend=refined$end1,beddata=refinedMol1)
# plotMol(chrom=refined$chr2,chromstart=refined$start2,chromend=refined$end2,beddata=refinedMol2)
# 
# 
# 
# 
# 
# 
# par(mfrow=c(3,1))
# chrom="chr14"
# chromstart=96000000
# chromend=97000000
# plotMol(chrom,chromstart,chromend,raw)
# chromstart=106000000
# chromend=107000000
# plotMol(chrom,chromstart,chromend,raw)
# chromstart=116000000
# chromend=117000000
# plotMol(chrom,chromstart,chromend,raw)
# 
# chromstart=106000000
# chromend=106200000
# plotMol(chrom,chromstart,chromend,raw)
# 
# chrom="chr11"
# chromstart=4295000
# chromend=4360000
# plotMol(chrom,chromstart,chromend,raw)
# 
# labelgenome(chrom,chromstart,chromend,n=3,scale="Mb")
# 
# 
# 
# data(package="Sushi")
