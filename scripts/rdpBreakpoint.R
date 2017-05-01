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

for(p in c("devtools", "optparse", "rjson", "ggplot2", "RColorBrewer", "cowplot", "igraph", "raster","httr","plyr","stringr")) myrequire(p)
#for(p in c("XML")) myrequire(p,version="3.96-0.2",repo="http://cran.us.r-project.org")
for(p in c("GenomicRanges", "ggbio", "Rsamtools")) myrequire(p,repo="Bioc")
for(p in c("primerTree")) myrequire(p,repo="https://github.com/chaelir/primerTree.git")
#use following to test if primerTree works
#mammals_16S = search_primer_pair(name='A','CGGTTGGGGTGACCTCGGA','GCTGTTATCCCTAGGGTAACT',num_aligns=1000, total_primer_specificity_mismatch=3)


####Options####
option_list2 <- list(
  make_option(c("-d", "--workDir"), default=".",
              help="working directory [default %default]"),
  make_option(c("-s", "--slopSize"), type="integer", default=1000,
              help="2 times the boundary padding bp [default %default]"),
  make_option(c("-p", "--padding"), default="XXXXXXXXXX",
              help="signature pattern for breakpoint [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)

parser <- OptionParser(usage = "%prog [options] makeContig infoFile refFile", option_list=option_list2)
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)
if(length(cmd$args)!=2 & !debug){ print_help(parser); quit() }
message("-Info: invoking command:",commandArgs(),"\n")
infoFile = cmd$args[1]
refFile = cmd$args[2]
workDir = cmd$options$workDir
slopSize = cmd$options$slopSize
padding = cmd$options$padding

####Debug####
if(debug){
  workDir="/Users/charlie/Box/01_InProgress/0110_manuscripts/011000_zoomx_analysis_lcx/SET-DPP10-Fusion/"
  infoFile="chr2_116365000_116390000_chr9_131450000_131465000.mol"
  refFile="~/ws/hg/hg19/human_g1k_v37.fasta"
  padding = "XXXXXXXXXX"
  slopSize = 1000
} else { 
  options(warn=-1)
}

setwd(workDir)
colScheme="Dark2"
supportTypes=c("molecule","rdpFF","rdpRR","rdpIM")
colors=brewer.pal(length(supportTypes), colScheme)
names(colors)=supportTypes
info=read.table(infoFile, header=T, stringsAsFactors = F)
strandDf=c("1","-1",".")
strandGR=c("+","-","*")
rdpInfo=info[grepl("rdp",info[["name"]]),]
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

#only trust the start location of mapping
hasSupportTypes=sapply(supportTypes, function(x) { any(grepl(x,info$supportType)) } )

#estimate breakpoint based on each supportType group (same supportType group)
#break1IsLeft=rep(NA,length(supportTypes))
#break2IsLeft=rep(NA,length(supportTypes))
if(tail(strsplit(refFile,'\\.')[[1]],1)=="gz") {
  refSeq=getSeq(FaFile(gzfile(refFile)))
} else {
  refSeq=getSeq(FaFile(refFile))
}

# contruct a graph of |rdps|*|rdps| vertices
# c1 & c2  => c
# T  &  F  => F
# T  &  T  => T

#we first group by type
#  in each group by type, 
for(i in seq_along(supportTypes)) {
  #
  s = supportTypes[i]
  rdps = info[info$supportType==s,]
  if(nrow(rdps)==0) next
  rdpRanges1 = IRanges(start=rdps$start1, end=rdps$start1+1)
  rdpRangePairs1 = expand.grid(rdpRanges1,rdpRanges1) #
  rdpRangePairOverlaps1 = matrix(poverlaps(rdpRangePairs1[,1], rdpRangePairs1[,2], maxgap=slopSize/2), nrow=length(rdpRanges1) )
  rdpRanges2 = IRanges(start=rdps$start2, end=rdps$start2+1)
  rdpRangePairs2 = expand.grid(rdpRanges2,rdpRanges2) #
  rdpRangePairOverlaps2 = matrix(poverlaps(rdpRangePairs2[,1], rdpRangePairs2[,2], maxgap=slopSize/2), nrow=length(rdpRanges2) )
  rdpRangePairOverlaps=rdpRangePairOverlaps1 & rdpRangePairOverlaps2
  rdpRangePairOverlapsGraph=graph_from_incidence_matrix(rdpRangePairOverlaps) #|rdps| by |rdps| graph
  rdpRangePairOverlapsGraphComponents=components(rdpRangePairOverlapsGraph)   #will assume 2*|rdps| nodes
  for(j in seq_len(rdpRangePairOverlapsGraphComponents$no)) { 
    sel=which(rdpRangePairOverlapsGraphComponents$membership[1:length(rdpRanges1)]==j) #only need first half
    selrdps=rdps[sel,]
    #if less than three rdps, we cannot determine direction
    #if more, we determine direction by find small difference size
    break1IsLeft=NA; break2IsLeft=NA;
    start1Diff=diff(sort(selrdps$start1))
    if(nrow(rdps)>2 & start1Diff[1]<tail(start1Diff,1)) break1IsLeft = TRUE
    if(nrow(rdps)>2 & start1Diff[1]>tail(start1Diff,1)) break1IsLeft = FALSE
    start2Diff=diff(sort(rdps$start2))
    if(nrow(rdps)>2 & start2Diff[1]<tail(start2Diff,1)) break2IsLeft = TRUE
    if(nrow(rdps)>2 & start2Diff[1]>tail(start2Diff,1)) break2IsLeft = FALSE
  }

  
}

#determine breakpoint automatically, we will just use the start of mapping


IMbr1=min(info[info$supportType=="rdpIM",]$start1)
IMbr2=min(info[info$supportType=="rdpIM",]$start2)
RRbr1=min(info[info$supportType=="rdpRR",]$start1)
RRbr2=min(info[info$supportType=="rdpRR",]$start2)
FFbr1=max(info[info$supportType=="rdpFF",]$start1)
FFbr2=max(info[info$supportType=="rdpFF",]$start2)

g1=ggplot(rd1GR) + geom_arrowrect(aes(y=name,group=name,fill=type,group.selfish=TRUE)) + theme(legend.position="bottom", legend.box="vertical", legend.direction = "horizontal") #box-v, item-h
g2=ggplot(rd2GR) + geom_arrowrect(aes(y=name,group=name,fill=type,group.selfish=TRUE))
legendHeight=0.3 #control legend area height
plots=list(g1,g2)
plotGrobs1=lapply(plots,ggbio:::Grob)
legend <- cowplot::get_legend(plotGrobs1[[1]])
plotGrobs2=lapply(plots,function(p) { ggbio:::Grob(p+theme(legend.position="none")) })
if(length(plotGrobs2)>1) { labels="AUTO" } else {labels=NULL}
prow <- cowplot::plot_grid( plotlist = plotGrobs2,
                   labels = labels,
                   nrow = 1
)
lrow <- cowplot::ggdraw() + cowplot::draw_grob(legend) #this only draws the legend, all aes has to be predefined.
pall = cowplot::plot_grid( prow, lrow, ncol=1, rel_heights=c(1,legendHeight))
pall



# chr2: 116376677 | chr2: 116368067 x chr 9: 131457051 | chr 9: 131457138 - chr 9: 131458625 x chr2: 116376628
# chr2: 116376677 x chr 9: 131457051 - chr 9: 131458625 x chr2: 116376628
# chr2: 116368067 x chr 9: 131457138 - chr 9: 131458625 x chr2: 116376628



min(end(rd2GR))
max(end(rd2GR))


