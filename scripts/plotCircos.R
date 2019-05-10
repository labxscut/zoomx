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

for(p in c("optparse", "rjson", "ggplot2", "RColorBrewer","cowplot","GenomicRanges", "ggbio", "GenomeInfoDb")) require(p,character.only=T)

option_list2 <- list(
  make_option(c("-o", "--outPrefix"), default="none",
              help="none=$moleSet1File1;$moleSet2File1;... [default %default]"),
  make_option(c("-l", "--eventLabels"), default="none",
              help="none=$moleSet1File1;$moleSet2File1;... [default %default]"),
  make_option(c("-t", "--titleText"), default="none",
              help="none=$outPrefix  [default %default]"),
  make_option(c("-r", "--reference"), default="hg19",
              help="use hg38 for GRCH 38 [default %default]"),
  make_option(c("-u", "--isUCSC"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-c", "--maxCnt"), default=75,
              help="maximum cnt in legend [default %default]"),
  make_option(c("-m", "--isMale"), action="store_true", default=FALSE,
              help="show Y or chrY, [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)

parser <- OptionParser(usage = "%prog [options] moleSet1File1,moleSet1File2,...;moleSet2File1,moleSet2File2,...;...",
                       option_list=option_list2)
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)
if(length(cmd$args)!=1 & !cmd$options$debug){ print_help(parser); quit() }
moleFiles=strsplit(unlist(strsplit(cmd$args[1],";")),",")

outPrefix=cmd$options$outPrefix
eventLabels=strsplit(cmd$options$eventLabels,";")[[1]]
maxCnt=cmd$options$maxCnt
isMale=cmd$options$isMale
titleText=cmd$options$titleText
reference=cmd$options$reference
isUCSC=cmd$options$isUCSC
debug=cmd$options$debug
noQuiet=cmd$options$noQuiet

if(outPrefix=="none") outPrefix="zoomx.circus.plot"
if(titleText=="none") { titleText=lapply(moleFiles,"[[",1) } else { titleText=list(strsplit(cmd$options$titleText,";")[[1]]) }
if(eventLabels=="none") eventLabels = c("IntraChr","InterChr")
  
if(debug){ 
  maxCnt = 75
  isMale = FALSE
  setwd("~/ws/crcbrain/03_Batch")
  outPrefix = list("tmp.somatic")
  moleFiles = list(c("P05378_mt.wg12.hg38x.bam.dist.E10_P05378_nt.wg12.hg38x.bam.dist.E10.somatic",
                     "P05683_mt.wg12.hg38x.bam.tra.E10_P05683_nt.wg12.hg38x.bam.tra.E10.somatic"))
  titleText = list("tmp.somatic")
  eventLabels = c("IntraChr","InterChr")
  reference = "hg38"
  isUCSC = TRUE
  #setwd("~/ws/crcbrain/02_ZoomX")
  #outPrefix = list("zoomx-Figure4")
  #moleFiles = list(c("P05378_7175_mt.wg12.bam.dist.E25_P05378_7176_nt.wg12.bam.dist.E10.somatic.curate","P05378_7175_mt.wg12.bam.tra.E25_P05378_7176_nt.wg12.bam.tra.E10.somatic.curate"))
  #titleText = list("MetB7175")
  #eventLabels = c("IntraChr","InterChr")
  #setwd("~/ws/10x6")
  #moleFiles = list(c("P05378_7175_mt.wg12.bam.dist.E25_P05378_7176_nt.wg12.bam.dist.E10.somatic","P05378_7175_mt.wg12.bam.dist.E25_P05378_7176_nt.wg12.bam.dist.E10.germline","P05378_7175_mt.wg12.bam.tra.E25_P05378_7176_nt.wg12.bam.tra.E10.somatic","P05378_7175_mt.wg12.bam.tra.E25_P05378_7176_nt.wg12.bam.tra.E10.germline"),
  #  c("P00526_2721_mt.st.10X_v2.bam.dist.E25_P00526_2386_nt.10X_v2.bam.dist.E10.somatic","P00526_2721_mt.st.10X_v2.bam.tra.E25_P00526_2386_nt.10X_v2.bam.tra.E10.somatic","P00526_2721_mt.st.10X_v2.bam.dist.E25_P00526_2386_nt.10X_v2.bam.dist.E10.germline","P00526_2721_mt.st.10X_v2.bam.tra.E25_P00526_2386_nt.10X_v2.bam.tra.E10.germline"),
  #  c("P00526_2725_mt.10X_v2.bam.dist.E25_P00526_2386_nt.10X_v2.bam.dist.E10.somatic","P00526_2725_mt.10X_v2.bam.tra.E25_P00526_2386_nt.10X_v2.bam.tra.E10.somatic","P00526_2725_mt.10X_v2.bam.dist.E25_P00526_2386_nt.10X_v2.bam.dist.E10.germline","P00526_2725_mt.10X_v2.bam.tra.E25_P00526_2386_nt.10X_v2.bam.tra.E10.germline"))
  #titleText = list("MetB7175","MetR2721","MetL2725")
  #eventLabels = c("IntraChr(Tumour)","InterChr(Tumour)","IntraChr(NormTiss)","InterChr(NormTiss)")
}

#options control subplots in ggbio
colScheme="Dark2"
circleTextSize=3.5
circleLineSize=1
mainTextSize=10
circleScaleSize=2
processedArgs=list("outPrefix"=outPrefix,"moleFiles"=moleFiles,"eventLabels"=eventLabels,"maxCnt"=maxCnt,"isMale"=isMale,"colScheme"=colScheme)
#handle both 10 and 11 columns input and convert as 9 column input
#regCols10=c("chrom1","start1","end1","chrom2","start2","end2","molcnt","IM","FF","RR") # in refined
#regCols11=c("chrom1","start1","end1","chrom2","start2","end2","molcnt","IM","FF","RR","BX") # in refined
regCols13=c("chrom1","start1","end1","chrom2","start2","end2","name","molcnt","FR","RF","FF","RR","MI") # in refined
combinedCols=c("chrom1","start1","end1","chrom2","start2","end2","molcnt","rdpcnt","type") 

print(processedArgs)
options(stringsAsFactors = FALSE)
ref=GenomeInfoDb::fetchExtendedChromInfoFromUCSC(reference)
nref=ifelse(isMale,24,23)
refinfo=ref$UCSC_seqlength[1:nref]
refnames=ifelse(rep(isUCSC,nref), ref$UCSC_seqlevel[1:nref], ref$NCBI_seqlevel[1:nref])
names(refinfo)=refnames
refsub=GRanges(seqnames=names(refinfo),
                ranges=IRanges(start=1, end=refinfo),
                seqlengths=refinfo)
colors=RColorBrewer::brewer.pal(length(eventLabels), colScheme)
names(colors) = eventLabels

plotCircos = function(idx, moleFiles, eventLabels, refsub) {
  message("processing fileSet",idx,":",moleFiles)
  raw = lapply(seq_along(moleFiles[[idx]]), FUN=function(x) { 
    tmp=read.table(moleFiles[[idx]][x], header=T, sep="\t")
    if(nrow(tmp)>0){
      #if(length(names(tmp))==10) tmp$BX=NA
      #names(tmp)=regCols11
      names(tmp) = regCols13
      frcnt = sapply(seq(nrow(tmp)), function(y) { length(fromJSON(tmp$FR[y])) })
      rfcnt = sapply(seq(nrow(tmp)), function(y) { length(fromJSON(tmp$RF[y])) })
      ffcnt = sapply(seq(nrow(tmp)), function(y) { length(fromJSON(tmp$FF[y])) })
      rrcnt = sapply(seq(nrow(tmp)), function(y) { length(fromJSON(tmp$RR[y])) })
      micnt = sapply(seq(nrow(tmp)), function(y) { length(fromJSON(tmp$MI[y])) })
      tmp$rdpcnt = frcnt+ffcnt+rrcnt+ffcnt
      tmp = tmp[,c("chrom1","start1","end1","chrom2","start2","end2","molcnt","rdpcnt")]
      tmp$type = eventLabels[x]
      names(tmp) = combinedCols
    } else {
      tmp=NULL
    }                                                    
    tmp
  })
  rawCombined = do.call("rbind",raw)
  
  if(is.null(rawCombined)) return(ggbio())

  convertedAll = rawCombined
  convertedAll$start1 = pmax(1,convertedAll$start1)
  convertedAll$start2 = pmax(1,convertedAll$start2)
  convertedAll$molcnt = pmin(convertedAll$molcnt,maxCnt)
  convertedAll = convertedAll[ ( convertedAll$chrom1 %in% refnames & convertedAll$chrom2 %in% refnames ), ]
  if(!isMale) convertedAll = convertedAll[ !( convertedAll$chrom1 %in% c("Y","chrY") | convertedAll$chrom2 %in% c("Y","chrY") ), ]
  
  if(nrow(convertedAll) == 0) return(ggbio())
  print(convertedAll)

  gr <- GRanges(seqnames = convertedAll$chrom1, 
                ranges = IRanges(start = convertedAll$start1, 
                                 end = convertedAll$end1),
                score = convertedAll$molcnt,
                linkTo = GRanges(seqnames = convertedAll$chrom2,
                                 ranges = IRanges(start = convertedAll$start2, 
                                                  end = convertedAll$end2),
                                 seqlengths = seqlengths(refsub)),
                Type = convertedAll$type,
                seqlengths = seqlengths(refsub))
  genome(gr) = reference

  grmp <- GRanges( seqnames = c(convertedAll$chrom1, convertedAll$chrom2),
                   ranges = IRanges(start = c(convertedAll$start1, convertedAll$start2), 
                                    end = c(convertedAll$end1, convertedAll$end2)),
                   Molecule =  round(rep(convertedAll$molcnt, 2) / 5) * 5, #ceiling to nearest 5
                   Type = rep(convertedAll$type, 2),
                   seqlengths = seqlengths(refsub))
  genome(grmp) = reference

  p <- ggbio() +
    circle(gr, geom = "link", linked.to = "linkTo", aes(color = Type), size = circleLineSize) +
    scale_color_manual(values=colors[match(eventLabels,names(colors))], guide=guide_legend(nrow=1)) +
    circle(grmp, geom = "point", aes(size = Molecule)) +
    scale_radius(breaks=seq(min(grmp$Molecule),max(grmp$Molecule),by=5),range=c(3,10),
                 guide=guide_legend(title="Number of Molecules", nrow=1)) +
    circle(refsub, geom = "ideo", fill = "gray70") +
    #circle(refsub, geom = "scale", size = circleScaleSize) +
    #theme(axis.text=element_blank(), axis.ticks=element_blank()) + #Scale
    circle(refsub, geom = "text", aes(label = seqnames), vjust = 0, size = circleTextSize) +
    ggtitle(titleText[[idx]]) +
    theme(plot.title = element_text(hjust=0.5, size=mainTextSize)) +
    theme(legend.position="bottom", legend.box="vertical", legend.direction = "horizontal") + #box-v, item-h
    theme(legend.margin = margin(c(0,0,0,0), unit="in")) + #reduces legend margin
    theme(legend.box.just="top") #this moves legend to the top of its bounding box
  p
}
#p=plotCircos(1,moleFiles,eventLabels,refsub)

#options control cowplot
legendHeight = 0.1 #control legend area height
cowLabelVjust = 8
cowLabelTextSize = mainTextSize
plots = lapply(seq(length(moleFiles)), plotCircos, moleFiles, eventLabels, refsub)
plotGrobs1=lapply(plots,ggbio:::Grob)
legend <- get_legend(plotGrobs1[[1]])
plotGrobs2=lapply(plots,function(p) { ggbio:::Grob(p+theme(legend.position="none")) })
if(length(plotGrobs2)>1) { printLabels="AUTO" } else { printLabels=NULL }
prow = plot_grid( plotlist = plotGrobs2, labels = printLabels, nrow = 1 )
lrow = ggdraw() + draw_grob(legend) #this only draws the legend, all aes has to be predefined.
pall = plot_grid( prow, lrow, ncol=1, rel_heights=c(1,legendHeight))
pdfFile=paste(outPrefix,"pdf",sep=".")
pdf(pdfFile,width=5.6*length(plots),height=5.6*(1+legendHeight)) #assuming length 

print(pall)

dev.off()


