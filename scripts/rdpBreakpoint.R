#!/usr/bin/env Rscript
setwd("/Users/charlie/Box/01_InProgress/0110_manuscripts/011000_zoomx_analysis_lcx/SET-DPP10-Fusion/")
info=read.table("chr2_116365000_116390000_chr9_131450000_131465000.mol", header=T, stringsAsFactors = F)
library(ggbio)
library(GenomicRanges)
library(cowplot)
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
                             strand=rdpInfo$strand2)

g1=ggplot(rd1GR) + geom_arrowrect(aes(y=name,group=name,fill=name,group.selfish=TRUE)) + theme(legend.position="bottom", legend.box="vertical", legend.direction = "horizontal") #box-v, item-h
g2=ggplot(rd2GR) + geom_arrowrect(aes(y=name,group=name,fill=name,group.selfish=TRUE))
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

#this maybe automated if we take the end with minimal variaiton
IMbr1=min(info[info$supportType=="rdpIM",]$start1)
IMbr2=min(info[info$supportType=="rdpIM",]$start2)
RRbr1=min(info[info$supportType=="rdpRR",]$start1)
RRbr2=min(info[info$supportType=="rdpRR",]$start2)
FFbr1=max(info[info$supportType=="rdpFF",]$start1)
FFbr2=max(info[info$supportType=="rdpFF",]$start2)

# chr2: 116376677 | chr2: 116368067 x chr 9: 131457051 | chr 9: 131457138 - chr 9: 131458625 x chr2: 116376628

# chr2: 116376677 x chr 9: 131457051 - chr 9: 131458625 x chr2: 116376628
# chr2: 116368067 x chr 9: 131457138 - chr 9: 131458625 x chr2: 116376628



min(end(rd2GR))
max(end(rd2GR))