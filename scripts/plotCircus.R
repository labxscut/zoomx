library("ggbio")
library("rjson")
library("ggplot2")
library("GenomicRanges")
data("CRC", package = "biovizBase")
setwd("~/ws/10x6")
prefix=c("P00526_2721_mt.st.10X_v2.bam","P00526_2725_mt.10X_v2.bam")
tag=c("dis","tra")
sample=c("MetR2721","MetL2725")
type=c("intra","inter")
textSize=12
options(stringsAsFactors = FALSE)
regCols10=c("chrom1","start1","end1","chrom2","start2","end2","molcnt","IM","FF","RR") # in refined
regCols11=c("chrom1","start1","end1","chrom2","start2","end2","molcnt","IM","FF","RR","BX") # in refined

pairs=expand.grid(prefix=prefix,tag=tag)
pairs$prefix=as.character(pairs$prefix)
pairs$tag=as.character(pairs$tag)
raw=apply(pairs, 1, FUN=function(r) { tmp=read.table(paste(r[1],r[2],"E25.refined.somatic",sep="."), header=T, sep="\t")
                                      if(nrow(tmp)>0){
                                        tmp$sample=sample[match(r[1],prefix)]
                                        tmp$tag=type[match(r[2],tag)]
                                        tmp$id=paste(tmp$tag, seq_len(nrow(tmp)), sep="")
                                      } else {
                                        tmp=NULL
                                      } 
                                    tmp
                                    })
rawCombined=do.call("rbind",raw)
combined=rawCombined[rawCombined$sample=="MetL2725",]
if(dim(combined)[2]==13) {
  colnames(combined)=c(regCols10, "sample", "type", "id") 
} else {
  colnames(combined)=c(regCols11, "sample", "type", "id") }

convertedAll=combined
#convertedAll=combined[combined$type=="intra",][1,]  #seqlengths must be specified for both link and linked.to GRanges 
subwidth=6
pdf("zoomx-Figure6D.pdf",width=subwidth,height=subwidth/2)
gr <- GRanges(seqnames = convertedAll$chrom1, 
              ranges = IRanges(start = convertedAll$start1, 
                               end = convertedAll$end1),
              score = convertedAll$molcnt,
              linkTo = GRanges(seqnames = convertedAll$chrom2,
                               ranges = IRanges(start = convertedAll$start2, 
                                                  end = convertedAll$end2),
                               seqlengths=seqlengths(hg19sub)),
              type = convertedAll$type,
              seqlengths=seqlengths(hg19sub))

grp <- GRanges( seqnames = c(convertedAll$chrom1, convertedAll$chrom2),
                ranges = IRanges(start = c(convertedAll$start1, convertedAll$start2), 
                                 end = c(convertedAll$end1, convertedAll$end2)),
                molecule =  rep(convertedAll$molcnt, 2),
                type = rep(convertedAll$type, 2),
                seqlengths=seqlengths(hg19sub)) 

p <- ggbio() +
  circle(gr, geom = "link", linked.to = "linkTo", aes(color = type)) +
  circle(grp, geom = "point", aes(size=molecule)) +
  circle(hg19sub, geom = "ideo", fill = "gray70") +
  circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = textSize) +
  theme(legend.text = element_text(size = textSize))

print(p)
dev.off()
