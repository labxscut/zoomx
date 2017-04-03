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

for(p in c("optparse", "rjson", "ggplot2")) myrequire(p)
for(p in c("GenomicRanges", "ggbio", "GenomeInfoDb")) myrequire(p, repo="bioC")

option_list2 <- list(
  make_option(c("-o", "--outPrefix"), default="none",
              help="none=$moleFile1 [default %default]"),
  make_option(c("-l", "--labels"), default="none",
              help="none=$moleFile1,$moleFile2,... [default %default]"),
  make_option(c("-m", "--maxCnt"), default=75,
              help="maximum cnt in legend [default %default]"),
  make_option(c("-f", "--isFemale"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)

parser <- OptionParser(usage = "%prog [options] moleFiles", option_list=option_list2)
cat("-Info: invoking command:",commandArgs(),"\n")
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)

outPrefix=cmd$options$outPrefix
labels=strsplit(cmd$options$outPrefix,",")[[1]]
maxCnt=cmd$options$maxCnt
isFemale=cmd$options$isFemale
debug=cmd$options$debug
noQuiet=cmd$options$noQuiet

if(length(cmd$args)!=1 & !debug){ print_help(parser); quit() }
moleFiles=strsplit(cmd$args[1],",")[[1]]

if(debug){ 
  setwd("~/ws/10x6")
  outPrefix = "P00526_2725_mt.10X_v2.bam.circus"
  moleFiles = c("P00526_2725_mt.10X_v2.bam.dist.E25_P00526_2386_nt.10X_v2.bam.dist.E10.somatic","P00526_2725_mt.10X_v2.bam.tra.E25_P00526_2386_nt.10X_v2.bam.tra.E10.somatic","P00526_2725_mt.10X_v2.bam.dist.E25_P00526_2386_nt.10X_v2.bam.dist.E10.germline","P00526_2725_mt.10X_v2.bam.tra.E25_P00526_2386_nt.10X_v2.bam.tra.E10.germline")
  labels = c("Somatic-IntraChr","Somatic-InterChr","Germline-IntraChr","Germline-InterChr")
  maxCnt = 100
  isFemale = TRUE
} else { 
  #options(warn=-1)
}
textSize=4
subwidth=10
options(stringsAsFactors = FALSE)
hg19=fetchExtendedChromInfoFromUCSC("hg19")
hg19info=hg19$UCSC_seqlength[1:24]
names(hg19info)=hg19$NCBI_seqlevel[1:24]
hg19sub=GRanges(seqnames=names(hg19info),
                ranges=IRanges(start=1, end=hg19info),
                seqlengths=hg19info)
pdfFile=paste(outPrefix,"pdf",sep=".")

#handle both 10 and 11 columns input and convert as 9 column input
regCols10=c("chrom1","start1","end1","chrom2","start2","end2","molcnt","IM","FF","RR") # in refined
regCols11=c("chrom1","start1","end1","chrom2","start2","end2","molcnt","IM","FF","RR","BX") # in refined
combinedCols=c("chrom1","start1","end1","chrom2","start2","end2","molcnt","rdpcnt","type")

raw=lapply(seq_along(moleFiles), FUN=function(x) { tmp=read.table(moleFiles[x], header=T, sep="\t")
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

pdf(pdfFile,width=subwidth,height=subwidth/2)

convertedAll=rawCombined
convertedAll$molcnt=pmin(convertedAll$molcnt,maxCnt)
convertedAll$rdpcnt=pmin(convertedAll$rdpcnt,maxCnt)
if(isFemale) convertedAll=convertedAll[ !(convertedAll$chrom1 %in% c("Y","chrY") | convertedAll$chrom2 %in% c("Y","chrY")),  ]
gr <- GRanges(seqnames = as.character(convertedAll$chrom1), 
              ranges = IRanges(start = convertedAll$start1, 
                               end = convertedAll$end1),
              score = convertedAll$molcnt,
              linkTo = GRanges(seqnames = as.character(convertedAll$chrom2),
                               ranges = IRanges(start = convertedAll$start2, 
                                                  end = convertedAll$end2),
                               seqlengths=seqlength(hg19sub)),
              type = convertedAll$type,
              seqlengths=seqlength(hg19sub))

grmp <- GRanges( seqnames = as.character(c(convertedAll$chrom1, convertedAll$chrom2)),
                 ranges = IRanges(start = c(convertedAll$start1, convertedAll$start2), 
                                  end = c(convertedAll$end1, convertedAll$end2)),
                 molecule =  rep(convertedAll$molcnt, 2),
                 type = rep(convertedAll$type, 2),
                 seqlengths=seqlength(hg19sub)) 

#grrp <- GRanges( seqnames = c(convertedAll$chrom1, convertedAll$chrom2),
#                 ranges = IRanges(start = c(convertedAll$start1, convertedAll$start2), 
#                                  end = c(convertedAll$end1, convertedAll$end2)),
#                 readpair =  rep(convertedAll$rdpcnt, 2),
#                 type = rep(convertedAll$type, 2),
#                 seqlengths=seqlengths(hg19sub)) 

p <- ggbio() +
  circle(gr, geom = "link", linked.to = "linkTo", aes(color = type)) +
  circle(grmp, geom = "point", aes(size=molecule)) +
  #circle(grrp, geom = "point", aes(size=readpair)) +
  circle(hg19sub, geom = "ideo", fill = "gray70") +
  circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = textSize/2) +
  theme(legend.text = element_text(size = textSize/2))

print(p)
dev.off()
