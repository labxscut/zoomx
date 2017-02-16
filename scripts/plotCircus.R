library("ggbio")
setwd("~/ws/10x6")
prefix=c("P00526_2721_mt.st.10X_v2.bam","P00526_2725_mt.10X_v2.bam")
tag=c("dis","tra")
sample=c("MetR2721","MetL2725")
type=c("intra","inter")
options(stringsAsFactors = FALSE)
regCols=c("chr1","start1","end1","chr2","start2","end2","molcnt","IM","FF","RR") # in refined
library("rjson")

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
combined=do.call("rbind",raw)
colnames(combined)=c(regCols, "sample", "tag", "id")
combined$chr1 = paste("chr",combined$chr1,sep="")
combined$chr2 = paste("chr",combined$chr2,sep="")

convertedAll=combined
colnames(convertedAll)[1:6]=c("chrom1","start1","end1","chrom2","start2","end2")

gr <- GRanges(seqnames = c(convertedAll$chrom1, convertedAll$chrom2), 
              ranges = IRanges(start = c(convertedAll$start1, convertedAll$start2), 
                               end = c(convertedAll$end1, convertedAll$end2)),
              #strand = factor(c(as.character(convertedAll$strand1),as.character(convertedAll$strand2)), levels=strandGR),
              #strand = c(convertedAll$strand1,convertedAll$strand2),
              #value = rep(0,2*nrow(convertedAll)), 
              score = c(convertedAll$score, convertedAll$score), 
              supportType = c(convertedAll$supportType, convertedAll$supportType), 
              pair = c(convertedAll$name, convertedAll$name))