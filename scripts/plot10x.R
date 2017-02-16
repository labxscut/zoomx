#!/usr/bin/env Rscript
library(optparse)
library(rhdf5)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

#library(RJSONIO)
#library(IRanges)
#library(GenomicRanges)
#library(ff)
#library(ffbase)
#library(vcd)
#library(MASS)
#options(bitmapType='cairo')

option_list <- list(
  make_option(c("-o", "--output"), default="parse10xh5.pdf", help="output file name, [default: %default] \n"),
  make_option(c("-t", "--title"), default="10X GemCode Statistics", help="output title, [default: %default] \n")
  )

parser <- OptionParser(
  usage="%prog [options] parse10x.a1.h5,parse10x.a2.h5,...\n Plotting multiple 10x background statistics", 
	option_list=option_list)

cat("-Info: invoking command:",commandArgs(),"\n")

args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, 
                 positional_arguments = TRUE)
parse10xh5_text = cmd$args[1]
master_pdf_text=cmd$options$output; 
master_title_text=cmd$options$title;
if(length(cmd$args)!=1){ 
	print_help(parser); 
	#parse10xh5_text="P04979_CR1.wg10.bam.multiFragBar.a30000.h5,P04979_CR1.wg10.bam.multiFragBar.a50000.h5,P04979_CR1.wg10.bam.multiFragBar.a70000.h5"
  parse10xh5_text="LS411N_CR1.wg12.bam.multiFragBar.h5.a0"
	master_pdf_text="P04979_CR1.wg10.bam.multiFragBar.h5.a0.pdf"
	master_title_text="P04979_CR1 10X Sequencing Statistics"
}

parse10xh5=strsplit(parse10xh5_text,',')[[1]] 
master_pdf=master_pdf_text
master_title=master_title_text
file.remove(master_pdf) #restart master_pdf

#initiate three plot handles for T(b), H(l) and G(d)
### PRETTYME: set plotting size and output names here ###
textsize=20
halfsize=10
quadsize=5
sidesize=1
fonttype=3
fontsize=10
colors=rainbow(length(parse10xh5))
names(colors)=parse10xh5
pdf(file=master_pdf, title=master_title, height=quadsize, width=halfsize+sidesize)

#DEBUG:
#  h="P04979_CR1.wg10.bam.a30000.h5"
fragLenDistr=list()
rdpSpaceDistr=list()
barCount=list()
for(hi in seq_along(parse10xh5)) {
	h = parse10xh5[hi]
  #htag = substr(h,1,15)
  htag = h
  fragLenDistr[[htag]]=as.data.frame(h5read(h,"S",bit64conversion='double'))
  fragLenDistr[[htag]][["sample"]]=htag
	rdpSpaceDistr[[htag]]=as.data.frame(h5read(h,"D",bit64conversion='double'))
	rdpSpaceDistr[[htag]][["sample"]]=htag
	barCount[[htag]]=as.data.frame(h5read(h,"T",bit64conversion='double'))
	barCount[[htag]][["sample"]]=htag
  cat("adding plot data",htag,"\n")
}
fragLenDistr=do.call("rbind",fragLenDistr)
rdpSpaceDistr=do.call("rbind",rdpSpaceDistr)
barCount=do.call("rbind",barCount)

fragLenDistrPlot = ggplot()
fragLenDistrPlot = fragLenDistrPlot + geom_point(data=fragLenDistr, aes(x=index,y=values,group=as.factor(sample), color=as.factor(sample), shape=as.factor(sample)))
fragLenDistrPlot = fragLenDistrPlot + geom_line(data=fragLenDistr, aes(x=index,y=values,group=as.factor(sample), color=as.factor(sample), linetype=as.factor(sample)))
fragLenDistrPlot = fragLenDistrPlot + ggtitle("Fragment Length Distribution (Density): S(l)")
fragLenDistrPlot = fragLenDistrPlot + labs(x="Fragmeng Length: l (Number of co-fragment reads)",y="Density: S") 
fragLenDistrPlot = fragLenDistrPlot + theme(legend.position="bottom")
fragLenDistrPlot = fragLenDistrPlot + xlim(c(1,NA))
grid.arrange(fragLenDistrPlot, ncol=1)
fragLenDistrPlot = fragLenDistrPlot + xlim(c(1,100))
grid.arrange(fragLenDistrPlot, ncol=1)

rdpSpaceDistrPlot = ggplot()
rdpSpaceDistrPlot = rdpSpaceDistrPlot + geom_point(data=rdpSpaceDistr, aes(x=index,y=values,group=sample, color=as.factor(sample), shape=as.factor(sample)))
rdpSpaceDistrPlot = rdpSpaceDistrPlot + geom_line(data=rdpSpaceDistr, aes(x=index,y=values,group=sample, color=as.factor(sample), linetype=as.factor(sample)))
rdpSpaceDistrPlot = rdpSpaceDistrPlot + ggtitle("Fragment Read Pair Spacing Distribution (Density): D(d)")
rdpSpaceDistrPlot = rdpSpaceDistrPlot + labs(x="Fragmeng Read Pair Spacing: d (Distance between adjacent co-fragment read pairs)",y="Density: D")
rdpSpaceDistrPlot = rdpSpaceDistrPlot + theme(legend.position="bottom")
grid.arrange(rdpSpaceDistrPlot, ncol=1)
rdpSpaceDistrPlot = rdpSpaceDistrPlot + xlim(c(0,1000))
grid.arrange(rdpSpaceDistrPlot, ncol=1)

barCountPlot = ggplot(data=barCount, aes(as.factor(sample), block0_values))
barCountPlot = barCountPlot + geom_violin()
barCountPlot = barCountPlot + labs(y="Barcode Cound Density: T", x="Samples")
barCountPlot = barCountPlot + theme(axis.text.x = element_text(angle=45,family="Helvetica", face="plain", colour="black", size=10))
grid.arrange(barCountPlot, ncol=1)
barCountPlot = barCountPlot + ylim(c(0,0.00001))
grid.arrange(barCountPlot, ncol=1)

dev.off()
quit()

#chrs=cmd$options$chrs
#DEBUG: parse10xh5="NA12878_WGS.chr22"; setwd('~/ws/10x/'); chrs="chr22"
#barrp=fromJSON(paste(parse10xh5,".barrp.json",sep="")) # {barcode:{r1:[(c1,p1),...], r2:[(c2,p2),...]}}
#brcount=read.table.ffdf(file=paste(parse10xh5,".brcount.dat",sep=""))
#cat("read hist_brcount\n")
#brspace=read.table.ffdf(file=paste(parse10xh5,".brspace.dat",sep=""))
#cat("read hist_brspace\n")
#brinsert=read.table.ffdf(file=paste(parse10xh5,".brinsert.dat",sep=""))
#cat("read hist_brinsert\n")
#mbrcount=read.table.ffdf(file=paste(parse10xh5,".mbrcount.dat",sep=""))
#cat("read hist_mbrcount\n")
#mbrspace=read.table.ffdf(file=paste(parse10xh5,".mbrspace.dat",sep=""))
#cat("read hist_mbrspace\n")
#mbrinsert=read.table.ffdf(file=paste(parse10xh5,".mbrinsert.dat",sep=""))
#cat("read hist_mbrinsert\n")
#
#ref=fromJSON(paste(parse10xh5,".ref.json",sep=""))
#barstat_pdf=paste(parse10xh5,".barstat.pdf",sep="")
#fragstat_pdf=paste(parse10xh5,".fragstat.pdf",sep="")
##master_png=paste(parse10xh5,".10x2.png",sep="")
#master_RData=paste(parse10xh5,".10x2.RData",sep="") 

#aethetic settings
textsize=20
halfsize=10
quadsize=5
sidesize=1
fonttype=3
fontsize=10

#utility functions
cp2g=function(chr,pos){
  ref$A[chr]+pos
}

g2cp=function(gin){
  chr=ref$chrname[tail(which(ref$acclen<cp2g('chr22',100)),1)] #if gin at acclen, give right
  list(chr=chr,pos=gin-ref$A[chr])
}

if( chrs %in% c("A","a","All","ALL","all") ){
  chrs = ref$chrname
} else {
  chrs = strsplit(chrs,",")[[1]]
}

#pdf(file=master_pdf, title="10X Statistics",height=4*quadsize,width=quadsize+sidesize)
pdf(file=barstat_pdf, title="10X Read-Barcode Analysis Statistics",width=2*textsize)
#png(file=master_png, title="10X Statistics")
par(mfrow=c(3,4))
par(mar=c(1.5, 4.5, 1.5, 0.5)) #bottom, left, top, right

#4. data structure for some statistics
# barrp: barcode read {b1: {c1:[c11,c12,...],p1:[p11,p12,....],c2:[c21,c22,...],p2:[p21,p22,...]}, b2: ...}

#TODO these functions are all whole genome based
#4.1 brcount
cat("Total Barcodes:",length(brcount$V1),"\n")
cat("Total Read Pairs (exclude hanging, low qual, multiple aligned):",length(brinsert$V1),"\n")
#brcount=lapply(barrp,function(x){length(x$c1)})
#summary(unlist(brcount))
#hist_brcount=hist(unlist(brcount),breaks=seq(0,max(unlist(brcount))),xlim=c(0,50),xlab="read pair count on barcode (pooled)",main="barcode-read pair count histo: hist_brcount")  
hist_brcount=hist.ff(brcount$V1,breaks=seq(0,max.ff(brcount$V1),1),freq=F,xlim=c(0,500),xlab="count",main=paste("hist_brcount: histogram of read pairs counts given the same barcode, total",length(brcount$V1)))  #1-1
hist_brcount_100=hist.ff(brcount$V1,breaks=seq(0,max.ff(brcount$V1),100),freq=F,xlim=c(0,5000),xlab="count",main=paste("hist_brcount_100: histogram of read pairs counts given the same barcode, total",length(brcount$V1)))  #1-2
#print(hist_brcount)
cat("done hist_brcount\n")

#4.2 brspace
#brspace=lapply(barrp,function(x){  diff(x$p1)[which(diff(x$c1)==0)]  })
#summary(unlist(brspace))
#hist_brspace=hist(unlist(brspace),breaks=seq(0,max(unlist(brspace))),xlim=c(0,500),xlab="read pair space on barcode (pooled)",main="barcode-read pair space histo: hist_brspace")
hist_brspace=hist.ff(brspace$V1,breaks=seq(0,max.ff(brspace$V1)),freq=F,xlim=c(0,500),xlab="space",main=paste("hist_brspace: histogram of read pairs spacing given the same barcode, total",length(brspace$V1))) #1-3
hist_brspace_10=hist.ff(brspace$V1,breaks=seq(0,max.ff(brspace$V1),10),freq=F,xlim=c(0,500),xlab="space",main=paste("hist_brspace_10: histogram of read pairs spacing given the same barcode, total",length(brspace$V1))) #1-4
#print(hist_brspace)
cat("done hist_brspace\n")

#4.3 brinsert
#brinsert=lapply(barrp,function(x){  abs((x$p1-x$p2)[which(x$c1==x$c2)])  })
#summary(unlist(brinsert))
#hist_brinsert=hist(unlist(brinsert),breaks=seq(0,max(unlist(brinsert))),xlim=c(0,500),xlab="read pair insert on barcode (pooled)",main="barcode-read pair insert histo: hist_brinsert") 
hist_brinsert=hist.ff(brinsert$V1,breaks=seq(0,max.ff(brinsert$V1)),freq=F,xlim=c(0,500),xlab="insert",main=paste("hist_brinsert: histogram of read pairs insert size given the same barcode, total",length(brinsert$V1))) #2-1
hist_brinsert_10=hist.ff(brinsert$V1,breaks=seq(0,max.ff(brinsert$V1),10),freq=F,xlim=c(0,500),xlab="insert",main=paste("hist_brinsert_10: histogram of read pairs insert size given the same barcode, total",length(brinsert$V1))) #2-2
#print(hist_brinsert)
cat("done hist_brinsert\n")

# bcount_single_percent=sum(unlist(bcount==1))/length(bcount) #singleton=66.5%
# #mostly one barcode two fragments; 370432 barcode -> 775300 fragment
# #?how to plot barcode-wise fragment size? --> F
# bsize=lapply(bar,function(x){ unlist(fsize[x+1]) }) #some have substantial more reads
# #?how to plot barcode-wise fragment spacing --> H?
# bpart=lapply(bar,function(x){ gin=sapply(frag[x+1],"[[",1); chrpos=lapply(gin,g2cp);
#                               part=sapply(tail(seq_along(chrpos),-1),
#                               function(x){ if(chrpos[[x]]$chr==chrpos[[x-1]]$chr) 
#                                       chrpos[[x]]$pos-chrpos[[x-1]]$pos else NA }) })
# brcount=lapply(bar,function(x){ gin=sapply(frag[x+1],"[[",1); length(unlist(gin)) })
# plot_brcount=hist(unlist(brcount),breaks=c(seq(0,100),max(unlist(brcount))),xlim=c(0,101),xlab="read on barcode",main="barcode-wise read count: plot_bcount")
# #bpart_tab=table(bpart)
# #bpart contains consecutive barcode spacing if on the same chromosome
# #plot_bpart=plot(as.numeric(bpart_tab[,1]),as.numeric(bpart_tab[,2]),xlab="fragment spacing on barcode",main="barcode-wise fragment spacing: plot_bpart")
# plot_bpart=hist(unlist(bpart),breaks=c(seq(0,100000,100),max(unlist(bpart))),xlim=c(0,100001),xlab="fragment spacing on barcode",main="barcode-wise fragment spacing: plot_bpart")
# #?how to plot barcode-wise read spacing distribution  --> G
# bread=lapply(bar,function(x){ unlist(fread[x+1]) })
# plot_bread=hist(unlist(bread),breaks=c(seq(-1,1000,10),max(unlist(bread),na.rm=T)),xlim=c(-1,1001),xlab="fragment spacing on barcode",main="barcode-wise read spacing: plot_bread")

#4.3 fragdata: genome region - coverage -> N
# frange=GRanges()
# for(chr in ref$chrname){ #lower bound in A and upper bound in B
#   fstart=lapply(frag,function(x){ if(x[[3]][1]>=ref$A[chr] & x[[3]][1]<ref$B[chr]) x[[3]][1]-ref$A[chr] else NA  })
#   fend=lapply(frag,function(x){ if(tail(x[[3]],1)>=ref$A[chr] & tail(x[[3]],1)<ref$B[chr]) tail(x[[3]],1)-ref$A[chr]+88 else NA })
#   fvalid=!is.na(fstart) & !is.na(fend)
#   if(length(fvalid)>0)
#     frange=c(frange,GRanges(seqname=chr,ranges=IRanges(start=fstart[fvalid],end=fend[fvalid])))
#   #Sadly, IRanges not handle 64-bit integer, need to convert to Genomic Range
# } #DONE, above is a slow algorithm; for all chromosome conversion use g2cp and construct GRanges from g2cp
# ftriple=lapply(frag,function(x){ schrpos=g2cp(x[[3]][1])
#                                  echrpos=g2cp(tail(x[[3]],1))
#                                  list(seqname=schrpos$chr,start=schrpos$pos,end=echrpos$pos) })
# frange=GRanges(seqname=sapply(ftriple,"[[",1),
#                 ranges=IRanges(start=sapply(ftriple,"[[",2),end=sapply(ftriple,"[[",3)))
# fcvg=unlist(lapply(coverage(frange),as.vector))

#print(fcvg[which(fcvg!=0)])
#summary(fcvg[fcvg!=0])
# fcvg_tab=as.data.frame(table(fcvg[fcvg!=0]))
# plot_fcvg=plot(as.numeric(fcvg_tab[,1]),as.numeric(fcvg_tab[,2]),xlim=c(1,500),xlab="fragment coverage",ylab="Frequency",main="genome-wide fragment coverage: plot_fcvg")

hist_mbrcount=hist.ff(mbrcount$V1,breaks=seq(0,max.ff(mbrcount$V1)),freq=F,xlim=c(0,500),xlab="count",main=paste("hist_mbrcount: histogram of read pairs counts given the same barcode, total",length(mbrcount$V1)))  #2-3
hist_mbrcount_100=hist.ff(mbrcount$V1,breaks=seq(0,max.ff(mbrcount$V1),100),freq=F,xlim=c(0,5000),xlab="count",main=paste("hist_mbrcount_100: histogram of read pairs counts given the same barcode, total",length(mbrcount$V1))) #2-4
#print(hist_mbrcount)
cat("done hist_mbrcount\n")

hist_mbrspace=hist.ff(mbrspace$V1,breaks=seq(0,max.ff(mbrspace$V1)),freq=F,xlim=c(0,500),xlab="space",main=paste("hist_mbrspace: histogram of read pairs spacing given the same barcode, total",length(mbrspace$V1))) #3-1
#hist.ff(as.ff(mbrspace$V1[1:1000]),breaks=seq(0,max.ff(as.ff(mbrspace$V1[1:1000]))),freq=F,xlim=c(0,500),xlab="space",main=paste("hist_mbrspace: histogram of read pairs spacing given the same barcode, total",length(mbrspace$V1[1:1000]))) #3-1
#lines(density(mbrspace$V1[1:1000],adjust=2))
#lines(density(mbrspace$V1,adjust=2))
#This is the G function we will use
#fit a density to it

hist_mbrspace_10=hist.ff(mbrspace$V1,breaks=seq(0,max.ff(mbrspace$V1),10),freq=F,xlim=c(0,500),xlab="space",main=paste("hist_mbrspace_10: histogram of read pairs spacing given the same barcode, total",length(mbrspace$V1))) #3-2

#print(hist_mbrspace)
cat("done hist_mbrspace\n") 
hist_mbrinsert=hist.ff(mbrinsert$V1,breaks=seq(0,max.ff(mbrinsert$V1)),freq=F,xlim=c(0,500),xlab="insert",main=paste("hist_mbrinsert: histogram of read pairs insert size given the same barcode, total",length(mbrinsert$V1))) #3-3
hist_mbrinsert_10=hist.ff(mbrinsert$V1,breaks=seq(0,max.ff(mbrinsert$V1),10),freq=F,xlim=c(0,500),xlab="insert",main=paste("hist_mbrinsert_10: histogram of read pairs insert size given the same barcode, total",length(mbrinsert$V1))) #3-4
#print(hist_mbrinsert)
cat("done hist_mbrinsert\n")
dev.off()

#pdf(file=master_pdf, title="10X Statistics",height=4*quadsize,width=quadsize+sidesize)
pdf(file=fragstat_pdf, title="10X Frag-Barcode Analysis Statistics",width=halfsize)
#png(file=master_png, title="10X Statistics")
par(mfrow=c(3,1))
par(mar=c(1.5, 4.5, 1.5, 0.5)) #bottom, left, top, right
alpha=c(30000,50000,70000)
#list of files to load:
#parse10xh5.mfrag.dat.alpha                  #mfrag size
#parse10xh5.mbarcode.dat.alpha               #mbarcode for mfrag
#parse10xh5.frag.dat.alpha                   #frag size
#parse10xh5.barcode.dat.alpha                #barcode for frag
mfragsize=list()
for(a in paste(alpha)){
  mfragsize[[a]]=read.table.ffdf(file=paste(parse10xh5,".mfrag.dat.",a,sep=""))
  hist_mfragsize=hist.ff(mfragsize[[a]]$V1,breaks=seq(0,max.ff(mfragsize[[a]]$V1)),freq=F,xlim=c(0,500),xlab="fragment size",main=paste("hist_mfragsize: histogram of fragment size given the same multi-fragment barcode, total",length(mfragsize[[a]]$V1))) #3-3
}
cat("done hist_mfragsize\n")
dev.off()

#pdf(file=master_pdf, title="10X Statistics",height=4*quadsize,width=quadsize+sidesize)
pdf(file=fragstat_pdf, title="10X Count-Barcode Analysis Statistics",width=halfsize)
#png(file=master_png, title="10X Statistics")
par(mfrow=c(3,2))
par(mar=c(1.5, 4.5, 1.5, 0.5)) #bottom, left, top, right
alpha=c(30000,50000,70000)
#list of files to load:
#parse10xh5.mfrag.dat.alpha                  #mfrag size
#parse10xh5.mbarcode.dat.alpha               #mbarcode for mfrag
#parse10xh5.frag.dat.alpha                   #frag size
#parse10xh5.barcode.dat.alpha                #barcode for frag
mfragsize=list()
for(a in paste(alpha)){
  mfragsize[[a]]=read.table.ffdf(file=paste(parse10xh5,".mfrag.dat.",a,sep=""))
  hist_mfragsize=hist.ff(mfragsize[[a]]$V1,breaks=seq(0,max.ff(mfragsize[[a]]$V1)),freq=F,xlim=c(0,500),xlab="fragment size",main=paste("hist_mfragsize: histogram of fragment size given the same multi-fragment barcode, total",length(mfragsize[[a]]$V1))) #3-3
}
cat("done hist_mfragsize\n")
mfragbar=list()
for(a in paste(alpha)){
  mfragbar[[a]]=read.table.ffdf(file=paste(parse10xh5,".mbarcode.dat.",a,sep=""))
  mfragbar_density=sort(as.vector(table.ff(mfragbar[[a]]$V1)),decreasing=T)
  mfragbar_density=mfragbar_density/sum(mfragbar_density)
  box_mfragbar=plot(mfragbar_density,xlab="barcode",main=paste("box_mfragbar: boxplot of multi-frag barcode density, total",length(mfragbar_density)))
}
cat("done box_mfragbar\n")
dev.off()

cat("save to", master_RData, "\n")
save.image(master_RData)

# fragment size? first peak? second peak? everything is cutted 12 or 14bp barcode 
# barcode: ==; genomic sequence: --; cut: |
#   theory for the first peak: cut for pair-end sequencing the almost same fragment
#                R1 ==----->
#                     |     |
#                     <-----== R2
#   collary: create first peak at 0bp
#   theory for the second peak: cut for sequencing exactly the same fragment
#                R1 ==----->
#                     |   |
#                   <-----== R2                 
#   collary: on read 2/1 create second peak at 12bp
#
# hanging read percentage?
#   unmapped read were filtered out (lowering the percentage based on flags)
#   many of only one end of paired read exist (25%)


#5 for all chromosome scan, iterate above for all chromosomes

#doing following in parallelized scheme:
# pick a grid (u,v): quad sect into four possibilities (u-,v-), (u-,v+), (u+,v-), (u+,v+)
# for each possibility: find frags in R(u) and R(v) and their bar B(u) and B(v)
# find barcode in intersect of B(u) and B(v)
# calculated the score S(u,v)

#eventually combine all scores S(u,v), pick the significant ones

#4.4 pick size of delta by compare Gb and NTb/T

#6 what are our noises; what are our signals
#  6.1 give signal if most barcode 

