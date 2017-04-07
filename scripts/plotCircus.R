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

for(p in c("optparse", "rjson", "ggplot2", "RColorBrewer","cowplot")) myrequire(p)
for(p in c("GenomicRanges", "ggbio", "GenomeInfoDb")) myrequire(p, repo="bioC")

option_list2 <- list(
  make_option(c("-o", "--outPrefix"), default="none",
              help="none=$moleSet1File1;$moleSet2File1;... [default %default]"),
  make_option(c("-l", "--labels"), default="none",
              help="none=$moleSet1File1;$moleSet2File1;... [default %default]"),
  make_option(c("-t", "--titleText"), default="none",
              help="none=$outPrefix  [default %default]"),
  make_option(c("-c", "--maxCnt"), default=75,
              help="maximum cnt in legend [default %default]"),
  #make_option(c("-w", "--withLegend"), action="store_true", default=FALSE,
  #            help="show legends, [default %default]"),
  #make_option(c("-m", "--isMale"), action="store_true", default=FALSE,
  #            help="show Y or chrY, [default %default]"),
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
#moleFiles=strsplit(cmd$args[1],";")[[1]]
moleFiles=strsplit(unlist(strsplit(cmd$args[1],";")),",")

outPrefix=cmd$options$outPrefix
labels=strsplit(cmd$options$labels,",")[[1]]
maxCnt=cmd$options$maxCnt
#isMale=cmd$options$isMale
#withLegend=cmd$options$withLegend
titleText=cmd$options$titleText
debug=cmd$options$debug
noQuiet=cmd$options$noQuiet

if(outPrefix=="none") outPrefix="zoomx.circus.plot"
if(titleText=="none") titleText=lapply(moleFiles,"[[",1)
if(debug){ 
  setwd("~/ws/10x6")
  outPrefix = list("zoomx-Figure4")
  moleFiles = list(c("P05378_7175_mt.wg12.bam.dist.E25_P05378_7176_nt.wg12.bam.dist.E10.somatic","P05378_7175_mt.wg12.bam.dist.E25_P05378_7176_nt.wg12.bam.dist.E10.germline","P05378_7175_mt.wg12.bam.tra.E25_P05378_7176_nt.wg12.bam.tra.E10.somatic","P05378_7175_mt.wg12.bam.tra.E25_P05378_7176_nt.wg12.bam.tra.E10.germline"),
    c("P00526_2721_mt.st.10X_v2.bam.dist.E25_P00526_2386_nt.10X_v2.bam.dist.E10.somatic","P00526_2721_mt.st.10X_v2.bam.tra.E25_P00526_2386_nt.10X_v2.bam.tra.E10.somatic","P00526_2721_mt.st.10X_v2.bam.dist.E25_P00526_2386_nt.10X_v2.bam.dist.E10.germline","P00526_2721_mt.st.10X_v2.bam.tra.E25_P00526_2386_nt.10X_v2.bam.tra.E10.germline"),
    c("P00526_2725_mt.10X_v2.bam.dist.E25_P00526_2386_nt.10X_v2.bam.dist.E10.somatic","P00526_2725_mt.10X_v2.bam.tra.E25_P00526_2386_nt.10X_v2.bam.tra.E10.somatic","P00526_2725_mt.10X_v2.bam.dist.E25_P00526_2386_nt.10X_v2.bam.dist.E10.germline","P00526_2725_mt.10X_v2.bam.tra.E25_P00526_2386_nt.10X_v2.bam.tra.E10.germline"))
  titleText = list("MetB7175","MetR2721","MetL2725")
  labels = c("IntraChr(MetTiss)","InterChr(MetTiss)","IntraChr(NormTiss)","InterChr(NormTiss)")
  maxCnt = 75
}

colScheme="Dark2"
pdfScale=1.4
circleTextSize=3.5
circleLineSize=1*pdfScale
mainTextSize=10*pdfScale
scaleSize=2
isMale=FALSE
processedArgs=list("outPrefix"=outPrefix,"moleFiles"=moleFiles,"labels"=labels,"maxCnt"=maxCnt,"isMale"=isMale,"colScheme"=colScheme)
#handle both 10 and 11 columns input and convert as 9 column input
regCols10=c("chrom1","start1","end1","chrom2","start2","end2","molcnt","IM","FF","RR") # in refined
regCols11=c("chrom1","start1","end1","chrom2","start2","end2","molcnt","IM","FF","RR","BX") # in refined
combinedCols=c("chrom1","start1","end1","chrom2","start2","end2","molcnt","rdpcnt","type") 

print(processedArgs)
options(stringsAsFactors = FALSE)
hg19=fetchExtendedChromInfoFromUCSC("hg19")
hg19info=hg19$UCSC_seqlength[1:24]
names(hg19info)=hg19$NCBI_seqlevel[1:24]
hg19sub=GRanges(seqnames=names(hg19info),
                ranges=IRanges(start=1, end=hg19info),
                seqlengths=hg19info)
colors=brewer.pal(length(labels), colScheme)
names(colors)=labels

#plotCircus(1,moleFiles,labels)
plotCircus = function(idx, moleFiles, labels) {
  raw=lapply(seq_along(moleFiles[[idx]]), FUN=function(x) { 
    tmp=read.table(moleFiles[[idx]][x], header=T, sep="\t")
    if(nrow(tmp)>0){
      if(length(names(tmp))==10) tmp$BX=NA
      names(tmp)=regCols11
      imcnt = sapply(seq(nrow(tmp)), function(y) { length(fromJSON(tmp[y,8])) })
      ffcnt = sapply(seq(nrow(tmp)), function(y) { length(fromJSON(tmp[y,9])) })
      rrcnt = sapply(seq(nrow(tmp)), function(y) { length(fromJSON(tmp[y,10])) })    
      tmp$rdpcnt=imcnt+ffcnt+rrcnt
      tmp=tmp[,c(1:7,12)]
      tmp$type=labels[x]
      names(tmp)=combinedCols
    } else {
      tmp=NULL
    }                                                    
    tmp
  })
  rawCombined=do.call("rbind",raw)

  convertedAll=rawCombined
  convertedAll$start1=pmax(1,convertedAll$start1)
  convertedAll$start2=pmax(1,convertedAll$start2)
  convertedAll$molcnt=pmin(convertedAll$molcnt,maxCnt)
  #convertedAll$rdpcnt=pmin(convertedAll$rdpcnt,maxCnt)
  if(!isMale) convertedAll=convertedAll[ !(convertedAll$chrom1 %in% c("Y","chrY") | convertedAll$chrom2 %in% c("Y","chrY")), ]
  
  print(convertedAll)

  gr <- GRanges(seqnames = convertedAll$chrom1, 
                ranges = IRanges(start = convertedAll$start1, 
                                 end = convertedAll$end1),
                score = convertedAll$molcnt,
                linkTo = GRanges(seqnames = convertedAll$chrom2,
                                 ranges = IRanges(start = convertedAll$start2, 
                                                  end = convertedAll$end2),
                                 seqlengths=seqlengths(hg19sub)),
                Type = convertedAll$type,
                seqlengths=seqlengths(hg19sub))
  genome(gr) = "hg19"

  grmp <- GRanges( seqnames = c(convertedAll$chrom1, convertedAll$chrom2),
                   ranges = IRanges(start = c(convertedAll$start1, convertedAll$start2), 
                                    end = c(convertedAll$end1, convertedAll$end2)),
                   Molecule =  rep(convertedAll$molcnt, 2),
                   Type = rep(convertedAll$type, 2),
                   seqlengths=seqlengths(hg19sub))
  genome(grmp) = "hg19"

  p <- ggbio() +
    circle(gr, geom = "link", linked.to = "linkTo", aes(color = Type), size=circleLineSize) +
    scale_color_manual(values=colors[match(labels,names(colors))]) +
    circle(grmp, geom = "point", aes(size=Molecule)) +
    circle(hg19sub, geom = "ideo", fill = "gray70") +
    #circle(hg19sub, geom = "scale", size = scaleSize) +
    #theme(axis.text=element_blank(), axis.ticks=element_blank()) + #Scale
    circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = circleTextSize) +
    ggtitle(titleText[[idx]]) +
    theme(plot.title = element_text(hjust=0.5, size=mainTextSize)) +
    #theme(legend.direction = "horizontal", legend.box="vertical", legend.position="bottom") + 
    theme(legend.box.just="top", legend.position="bottom", legend.box="vertical") + 
    theme(legend.text = element_text(size = mainTextSize-2)) + 
    theme(legend.title = element_text(size = mainTextSize-2)) +
    theme(legend.box.spacing = unit(c(0,0,0,0), "pt")) +
    theme(plot.margin = unit(c(0,0,0,0), "pt")) +
    theme(panel.spacing = unit(c(0,0,0,0), "pt"))
  
  p
}

cowLabelVjust=8*pdfScale
cowLabelTextSize=mainTextSize/pdfScale
plots=lapply(seq(length(moleFiles)),plotCircus,moleFiles,labels)
plotGrobs1=lapply(plots,ggbio:::Grob)
legend <- get_legend(plotGrobs1[[1]])
plotGrobs2=lapply(plots,function(p) { ggbio:::Grob(p+theme(legend.position="none")) })
#par(mar=c(0,0,0,0))
prow <- plot_grid( plotlist = plotGrobs2,
                   align = 'v',
                   labels = LETTERS[seq(length(plotGrobs2))],
                   label_size = mainTextSize,
                   hjust = -0.5,
                   vjust = cowLabelVjust,
                   nrow = 1
)
prow = prow + draw_grob(legend,x=0,y=-0.37,width=1,height=1) +
       theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, "pt")) +
       theme(panel.spacing = margin(t = 0, r = 0, b = 0, l = 0, "pt")) +
       theme(legend.box.spacing = unit(c(0,0,0,0), "pt"))

pdfFile=paste(outPrefix,"pdf",sep=".")
pdf(pdfFile,width=8.5*pdfScale)
print(prow)
dev.off()


