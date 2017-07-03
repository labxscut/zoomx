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

#for(p in c("devtools", "optparse", "rjson", "ggplot2", "RColorBrewer", "cowplot", "igraph", "raster","httr","plyr","stringr")) myrequire(p)
for(p in c("devtools", "optparse", "rjson", "ggplot2", "RColorBrewer", "cowplot", "igraph", "raster","httr","plyr","stringr")) require(p,character.only=T)
#for(p in c("XML")) myrequire(p,version="3.96-0.2",repo="http://cran.us.r-project.org")
for(p in c("GenomicRanges", "ggbio", "Rsamtools", "Biostrings","BSgenome.Hsapiens.UCSC.hg19","BSgenome.Hsapiens.UCSC.hg38")) require(p,character.only=T)
#for(p in c("GenomicRanges", "ggbio", "Rsamtools", "Biostrings","BSgenome.Hsapiens.UCSC.hg19","BSgenome.Hsapiens.UCSC.hg38")) myrequire(p,repo="Bioc")
#for(p in c("primerTree")) myrequire(p,repo="https://github.com/chaelir/primerTree.git")

#use following to test if primerTree works
#mammals_16S = search_primer_pair(name='A','CGGTTGGGGTGACCTCGGA','GCTGTTATCCCTAGGGTAACT',num_aligns=1000, total_primer_specificity_mismatch=3)
#devtools::source_url("https://gist.github.com/chaelir/raw/callPrimer3.R")

####Options####
option_list2 <- list(
  make_option(c("-o", "--outPrefix"), default="none",
              help="output prefix [default %default]"),
  make_option(c("-r", "--reference"), default="hg37",
              help="use hg19 for UCSC build [default %default]"),
  make_option(c("-s", "--breakSlopSize"), type="integer", default=1000,
              help="2 times the boundary padding bp [default %default]"),
  make_option(c("-p", "--padding"), default="-",
              help="signature pattern for breakpoint [default %default]"),
  make_option(c("-m", "--includeMolecule"), action="store_true", default=FALSE,
              help="including molecules in the plot [default %default]"),
  make_option(c("-q", "--noQuiet"), action="store_true", default=FALSE,
              help="show verbose, [default %default]"),
  make_option(c("-a", "--debug"), action="store_true", default=FALSE,
              help="save debug, [default %default]")
)

parser <- OptionParser(usage = "%prog [options] workDir", option_list=option_list2)
args <- commandArgs(trailingOnly = TRUE)
cmd = parse_args(parser, args, print_help_and_exit = TRUE, positional_arguments = TRUE)
if(length(cmd$args)!=1 & !cmd$options$debug){ print_help(parser); quit() }
message("-Info: invoking command:",commandArgs(),"\n")
debug = cmd$options$debug
workDir = cmd$args[1]
reference = cmd$options$reference
breakSlopSize = cmd$options$breakSlopSize
padding = cmd$options$padding
outPrefix = cmd$options$outPrefix
includeMolecule = cmd$options$includeMolecule

####Debug####
if(debug){
  #workDir="~/ws/zoomx/02_ZoomX/NA12878_from10X_wg12.bam.cis.E40.refined.Plot_1"
  workDir="~/ws/zoomx/02_ZoomX/NA12878_from10X_wg12.bam.cis.E40.refined.Plot_2"
  #workDir="~/ws/zoomx/02_ZoomX/NA12878_from10X_wg12.bam.cis.E40.refined.Plot_112"
  includeMolecule=TRUE
  outPrefix="none"
  reference="hg37"
  padding = "-"
  breakSlopSize = 1000
  moleEndFlushSize = 500
} else { 
  options(warn=-1)
}

setwd(workDir)
moleFiles=list.files(".",".mol$")
moleFile=moleFiles[!grepl("refined",moleFiles)][1]
callFile="refined.bedpe"
callInfo=read.table(callFile,sep="\t") 

if(outPrefix=="none") outPrefix=moleFile
colScheme="Dark2"
supportTypes=c("rdpFF","rdpRR","rdpIM","molecule")
colors=RColorBrewer::brewer.pal(length(supportTypes), colScheme)
names(colors)=supportTypes
info=read.table(moleFile, header=T, stringsAsFactors = F)
#retrieve sequences from reference
if(reference=="hg38") {
  refSeq = BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
} else { #fall back to hg19/hg37
  refSeq = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
}
strandDf=c("1","-1",".")
strandGR=c("+","-","*")
#rdpInfo=info[grepl("rdp",info[["name"]]),]
rdpInfo=info
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

if(!includeMolecule) { #analyzing only
  rd1GR = rd1GR[rd1GR$type!="molecule"]
  rd2GR = rd2GR[rd1GR$type!="molecule"]
} else {
  sel = (end(ranges(rd2GR))!=end(ranges(rd1GR)) & start(ranges(rd2GR))!=start(ranges(rd1GR))) | rd1GR$type!="molecule"
}

#analyzing only non-duplicates
rd1GR = rd1GR[sel]
rd2GR = rd2GR[sel]

if(length(rd1GR)==0 | length(rd2GR)==0) {
  message("no supporting read pairs and/or molecules, quitting")
  quit()
}

g1 = ggplot2::ggplot() + 
  ggbio::geom_arrowrect(data=rd1GR,ggplot2::aes(y=name,group=name,fill=type,group.selfish=TRUE),xlab=as.character(GenomicRanges::seqnames(rd1GR)[1])) + 
  ggplot2::theme(legend.position="bottom", legend.box="vertical", legend.direction = "horizontal") #box-v, item-h
g2 = ggplot2::ggplot() + 
  ggbio::geom_arrowrect(data=rd2GR,ggplot2::aes(y=name,group=name,fill=type,group.selfish=TRUE),xlab=as.character(GenomicRanges::seqnames(rd2GR)[1]))
legendHeight=0.3 #control legend area height
plots=list(g1,g2)
plotGrobs1=lapply(plots,ggbio:::Grob)
legend <- cowplot::get_legend(plotGrobs1[[1]])
plotGrobs2=lapply(plots,function(p) { ggbio:::Grob(p+ggplot2::theme(legend.position="none")) })
if(length(plotGrobs2)>1) { labels="AUTO" } else {labels=NULL}
prow <- cowplot::plot_grid( plotlist = plotGrobs2,
                            labels = labels,
                            nrow = 1
)
lrow <- cowplot::ggdraw() + cowplot::draw_grob(legend) #this only draws the legend, all aes has to be predefined.
pall = cowplot::plot_grid( prow, lrow, ncol=1, rel_heights=c(1,legendHeight))
message("saving support plots to:", paste(outPrefix,".pdf",sep="") )
cowplot::save_plot(paste(outPrefix,".pdf",sep=""), pall, base_height = 10, base_aspect_ratio = 2.1)

#convert support groups to contigs

#estimate breakpoint based on each supportType group (same supportType group)
#break1IsLeft=rep(NA,length(supportTypes))
#break2IsLeft=rep(NA,length(supportTypes))

# contruct a graph of |rdps|*|rdps| vertices
# c1 & c2  => c
# T  &  F  => F
# T  &  T  => T

#we first group by type
#  in each group by type, 
#  debug not finished
contigs=list()
for(i in seq_along(supportTypes)) {
  # i=1; i=2; i=3; i=4
  s = supportTypes[i]
  if( s != "molecule") { #rdp support
    if(s=="rdpFF") {
      strand1 = "+"; strand2 = "+"
    } else if(s=="rdpRR") {
      strand1 = "-"; strand2 = "-"
    } else {
      strand1 = c("-","+"); strand2 = c("+","-")
    }
    rdps = info[info$supportType==s,]
    if(nrow(rdps)==0) next #no such rdp support
    rdpRanges1 = IRanges(start=rdps$start1, end=rdps$start1+1)
    rdpRangePairs1 = expand.grid(rdpRanges1,rdpRanges1) #
    rdpRangePairOverlaps1 = matrix(poverlaps(rdpRangePairs1[,1], rdpRangePairs1[,2], maxgap=breakSlopSize/2), nrow=length(rdpRanges1) )
    rdpRanges2 = IRanges(start=rdps$start2, end=rdps$start2+1)
    rdpRangePairs2 = expand.grid(rdpRanges2,rdpRanges2) #
    rdpRangePairOverlaps2 = matrix(poverlaps(rdpRangePairs2[,1], rdpRangePairs2[,2], maxgap=breakSlopSize/2), nrow=length(rdpRanges2) )
    rdpRangePairOverlaps=rdpRangePairOverlaps1 & rdpRangePairOverlaps2
    rdpRangePairOverlapsGraph=graph_from_incidence_matrix(rdpRangePairOverlaps) #|rdps| by |rdps| graph
    rdpRangePairOverlapsGraphComponents=components(rdpRangePairOverlapsGraph)   #will assume 2*|rdps| nodes
    for(j in seq_len(rdpRangePairOverlapsGraphComponents$no)) { 
      sel=which(rdpRangePairOverlapsGraphComponents$membership[1:length(rdpRanges1)]==j) #only need first half
      selrdps=rdps[sel,]
      #if(nrow(selrdps)<2) next #require at least two rdps 
      #if less than three rdps, we cannot determine direction
      #if more, we determine direction by find small difference size
      #break1IsLeft is NA if neither Left or Right, is False if Right
      direction1="U"; direction2="U";
      if(nrow(selrdps)>2) { #rdp >=3, determine direction by read mapping density
        start1Diff=diff(sort(selrdps$start1)) #sorted
        if(start1Diff[1]<tail(start1Diff,1)) direction1 = "L"
        if(start1Diff[1]>tail(start1Diff,1)) direction1 = "R"
        start2Diff=diff(sort(selrdps$start2)) #sorted
        if(start2Diff[1]<tail(start2Diff,1)) direction2 = "L"
        if(start2Diff[1]>tail(start2Diff,1)) direction2 = "R"
      }
      sel1 = paste(strand1,direction1) %in% c("+ R", "- L", "+ U", "- U") #possible read1 configuration
      sel2 = paste(strand2,direction2) %in% c("+ R", "- L", "+ U", "- U") #possible read2 configuration
      sel = sel1 & sel2                                                   #only select possible ones
      strand1 = strand1[sel]; strand2 = strand2[sel];                     
      direction1 = ifelse(strand1=="+","R","L"); direction2 = ifelse(strand2=="+","R","L"); 
      #once selected convert U to L/R based on +/-
      #direction1 = direction1[sel]; direction2 = direction2[sel]; 
      start1 = (strand1=="+")*max(selrdps$start1) + (strand1=="-")*min(selrdps$start1); end1=start1
      start2 = (strand2=="+")*max(selrdps$start2) + (strand2=="-")*min(selrdps$start2); end2=start2
      for(k in seq_along(start1)) { # handles the two possibility of "IM" type read pairs
        c=list(chrom1=selrdps$chrom1[1],start1=start1[k],end1=end1[k],
               chrom2=selrdps$chrom2[1],start2=start2[k],end2=end2[k],
               name=paste("contig","type",i,"group",j,"possibility",k,sep="-"), score=nrow(selrdps),
               strand1=strand1[k],strand2=strand2[k],direction1=direction1[k],direction2=direction2[k],
               type=s, support=nrow(selrdps))
        contigs[[length(contigs)+1]]=c
      }
    }
  } else { #molecule support, no clustering needed
    next #too complex for molecule, unsupported for now
    moles = info[info$supportType=="molecule",]
    moleLeftFlush1 = sd(moles$start1) < moleEndFlushSize #first seg left flush <500
    moleRightFlush1 = sd(moles$end1) < moleEndFlushSize #first seg right flush <500
    moleLeftFlush2 = sd(moles$start2) < moleEndFlushSize #second seg left flush <500
    moleRightFlush2 = sd(moles$end2) < moleEndFlushSize #second seg right flush <500
    #deal 16 combinations of flushing
    if(moleLeftFlush1 & moleRightFlush1 & !moleLeftFlush2 & !moleRightFlush2 ) { 
      # a short translocated fragment, derive both start and end contigs
      # mol2L - mol1L - mol1R - mol2L
      # r(mol2R) - mol1L - mol1R - r(mol2L)
      c1=list(chrom1=)
    } else if(!moleLeftFlush1 & !moleRightFlush1 & moleLeftFlush2 & moleRightFlush2 ) {
      # a short translocated fragment, derive both start and end contigs
      
    } else if(moleLeftFlush1 & !moleRightFlush1 & moleLeftFlush2 & !moleRightFlush2) { 
      # derive one contig
      
    } else if(moleLeftFlush1 & !moleRightFlush1 & !moleLeftFlush2 & moleRightFlush2) { 
      # derive one contig
      
    } else if(!moleLeftFlush1 & moleRightFlush1 & moleLeftFlush2 & !moleRightFlush2) {
      # derive one contig
    } else if(!moleLeftFlush1 & moleRightFlush1 & !moleLeftFlush2 & moleRightFlush2) {
      # derive one contig
    } else {
      message("Warning: unsupported molecule configuration")
      next
    }
  }
}

#determine breakpoint automatically, we will just use the start of mapping

# a manual procedure example of rdpRR, as in contig2
# diff(info[info$supportType=="rdpRR",]$start1) #find outlier in differ, 1 in this case, exlude it
# diff(info[info$supportType=="rdpRR",]$start1[2:7]) #find the flushing ends, left in this case
# min(info[info$supportType=="rdpRR",]$start1[2:7]) #find the leftmost of flushing ends, 116376679
# posSort=sort(info[info$supportType=="rdpRR",]$start2, index.return=T) #if ends not sorted, sort first 
# excludeIndex=posSort$ix[which(diff(posSort$x)>500)+1] #find outlier (>500) by differ, 1 in this case
# Here we can use the results from start1 ends
# diff(info[info$supportType=="rdpRR",]$start2[-excludeIndex]) #find the flushing ends, left in this case
# min(info[info$supportType=="rdpRR",]$start2[2:7]) #find the leftmost of flushing ends, 131457138
# posSort=sort(info[info$supportType=="rdpIM",]$start1, index.return=T) #if ends not sorted, sort first 
# excludeIndex=posSort$ix[which(diff(posSort$x)>500)+1] #nothing to be excluded
# diff(sort(info[info$supportType=="rdpIM",]$start1))  #find the flushing ends, left in this case
# min(sort(info[info$supportType=="rdpIM",]$start1))

# determine from molecule
# find flushing ends
#   sd(info[info$supportType=="molecule",]$start1) #first seg left flush <500  
#   sd(info[info$supportType=="molecule",]$end1) #first seg right flush <500
#   sd(info[info$supportType=="molecule",]$start2) #second seg left flush <500  
#   sd(info[info$supportType=="molecule",]$end2) #second seg right flush <500
# if one end bloth flush, determine short segment by taking [min(start),max(end)]
#   moleStart2 = min(info[info$supportType=="molecule",]$start2)
#   moleEnd2 = max(info[info$supportType=="molecule",]$end2)
#   #deduce the two breakpoints in long segment by taking the median (median is not a bad guess)
#   moleStart1 = round(median(info[info$supportType=="molecule",]$start1+info[info$supportType=="molecule",]$end1)/2)-1500 #116376483 #KOD-hifi 1-5kb; Taq 500bp; Phusion 20kb
#   moleEnd1 = round(median(info[info$supportType=="molecule",]$start1+info[info$supportType=="molecule",]$end1)/2)+1500
#   in this case, a maxSlopSize is bounded by the size of short segment 
#   maxSlopSize = max(info[info$supportType=="molecule",]$end2-info[info$supportType=="molecule",]$start2)
    #all 4 contig combinations: 
    #              moleStart1 (R) x moleStart2 (L)
    #              moleStart1 (R) x moleEnd2 (R)
    #              moleStart2 (L) x moleEnd1 (L)
    #              moleEnd2 (R)   x moleEnd1 (L)

# chr2: 116376677 | chr2: 116368067 x chr9: 131457051 | chr 9: 131457138 - chr 9: 131458625 x chr2: 116376628

# chr2: 116376628 x chr9: 131457051 # [contig1: FFbr1 + right x FFbr2 + right] 
# chr2: 116368067 x chr9: 131457138 # [contig2: RRbr1 - left x RRbr2 - right] 
# chr2: 116376677 x chr2: 131457051 # [contig3: IMbr1 * left x IMbr2 * unknown]
#
# chr2: 116368067 x chr9: 131457138 - chr9: 131458625 x chr2: 116376628 [RRbr1 - RRbr2]

#NOTE breakSlopSize needs to be capped;
#direction is where the breakpoint is located relative to the given the coordinate
#L is left, R is right, U is unknown
# contig1 = list(chrom1 = "chr2", 
#                start1 = 116376628, 
#                end1 = 116376628,
#                strand1 = "+",
#                direction1 = "R",
#                chrom2 = "chr9",
#                start2 = 131458625,
#                end2 = 131458625, 
#                strand2 = "+",
#                direction2 = "R")
# contig2 = list(chrom1 = "chr2", 
#                start1 = 116368067,
#                end1 = 116368067,
#                strand1 = "-",
#                direction1 = "L",
#                chrom2 = "chr9",
#                start2 = 131458503,
#                end2 = 131458503,
#                strand2 = "-",
#                direction2 = "L"
#                )
# contig3 = list(chrom1 = "chr2", 
#                start1 = 116376677,
#                end1 = 116376677,
#                strand1 = "*",
#                direction1 = "L",
#                chrom2 = "chr9",
#                start2 = 131457051,
#                end2 = 131457051,
#                strand2 = "*",
#                direction2 = "U"
# )
#contigs = list(contig1,contig2,contig3)

#print(contigs)

contigsDF = Reduce(rbind, lapply(contigs,as.data.frame)) #save contigs
contigsDF = cbind(contigsDF, callInfo)
write.table(contigsDF, file=paste(outPrefix,"ctg",sep="."), sep="\t", quote=F, row.names=F, col.names=F)

contigSeqs = x <- vector(mode = "list", length = length(contigs))
#tmpSlopSize = min(breakSlopSize,maxSlopSize)
tmpSlopSize = breakSlopSize

formatCoordinate=function(chrom1,start1,end1,chrom2,start2,end2,type,support){
  paste(chrom1,start1,end1,chrom2,start2,end2,type,support,sep="|")
}

formContig=function(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,tag,type,support){
  contigSeqInfo = formatCoordinate(seg1chr, seg1start, seg1end, seg2chr, seg2start, seg2end,type,support)
  contigName = paste(tag,contigSeqInfo,sep="")
  contigSeg1 = refSeq[[seg1chr]][seg1start:seg1end]
  if(seg1rc) contigSeg1 = reverseComplement(contigSeg1)
  contigSeg2 = refSeq[[seg2chr]][seg2start:seg2end]
  if(seg2rc) contigSeg2 = reverseComplement(contigSeg2)
  setNames(list(xscat(contigSeg1, DNAString(x=padding), contigSeg2)), contigName)
}

retrieveContig = function(i, contigs, breakSlopSize, padding, refSeq) {
  contig = contigs[[i]]
  contigSeqs = NULL
  #contigInfo = NULL
  
  #since the breakpoint is in the middle
  #+ must be associated with read with break point to the Right
  #- must be associated with read with break point to the Left 
  #by paired-end read sequencing only (except palindromes)
  #only following patterns are expected
  #only need to consider
  #(-,L) x (+,R) 
  #(+,R) x (-,L)
  #(+,R) x (+,R)
  #(-,L) x (-,L)
  #(*,L) x (*,L) expands to (-,L) x (-,L)
  #(*,R) x (*,R) expands to (+,R) x (+,R)
  #(*,L) x (*,U) expands to (-,L) x (-,L) or (-,L) x (+,R)
  #(*,U) x (*,L) expands to (-,L) x (-,L) or (+,R) x (-,L)
  #(*,R) x (*,U) expands to (+,R) x (+,R) or (+,R) x (-,L)
  #(*,U) x (*,R) expands to (-,L) x (+,R) or (+,R) x (+,R)
  #(*,U) x (*,U) expands to (+,R) x (+,R), (+,R) x (-,L), (-,L) x (+,R) and (-,L) x (-,L)
  #let's say the two segment reference sequences are
  #     seg1        seg2
  # 5-ACACGTGTGA-ATATCACAGC-3
  # 3-TGTGCACACT-TATAGTGTCG-5
  #(+,R1) x (-,L2) implies contig1
  # 5-ACACGTGTGA-ATATCACAGC-3
  #(-,L1) x (+,R2) implies contig2
  # 5-ATATCACAGC-ACACGTGTGA-3
  #(+,R1) x (+,R2) implies
  # 5-ACACGTGTGA-GCTGTGATAT-3
  #(+,R2) x (+,R1) implies
  # 5-ATATCACAGC-AGTGTGCACA-3
  #(-,L1) x (-,L2) implies
  # 5-TCACACGTGT-ATATCACAGC-3
  #(-,L2) x (-,L1) implies
  # 5-GCTGTGATAT-ACACGTGTGA-3
  #contig:
  #  [seg1] - read => x <= read - [seg2]; 
  #reads are part of segments
  #read order is important, let's say the order is (r1,r2) and (r2,r1).
  #we assign the order to (-)/(-), (-)/(+) and (+)/(+), 
  #then we have (1-,2-),(2-,1-),(1-,2+),(1+,2-),(1+,2+),(2+,1+) six cases
  #NOTE: for L-R pairs, mapping orders are obvious; for L-L, R-R pairs, 1 and 2 are added to dnote their mapping order
  
  #generate (+,R) x (-,L) contig for compatible modes 
  # [seg1] | R+ => x <= L- | [seg2]
  if(paste(contig$strand1,contig$direction1) %in% c("+ R","+ U","* R","* U")  
     && paste(contig$strand2,contig$direction2) %in% c("- L","- U","* L","* U") ) {
    seg1chr=contig$chrom1; seg1start=contig$start1-breakSlopSize+1; seg1end=contig$end1; seg1rc=FALSE
    seg2chr=contig$chrom2; seg2start=contig$start2; seg2end=contig$end2+breakSlopSize-1; seg2rc=FALSE
    contigSeqs = append(contigSeqs, formContig(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,paste("contig-",i,"(R+,L-)|",sep=""),contig$type,contig$support))
    #contigInfo = append(contigInfo, list(chrom1=contig$chrom1, start1=contig$start1, end1=contig$end1, 
  }
  #generate (-,L) x (+,R) contig for compatible modes
  # [seg1] | R+ => x <= L- | [seg2]
  if(paste(contig$strand1,contig$direction1) %in% c("- L","- U","* L","* U") 
     && paste(contig$strand2,contig$direction2) %in% c("+ R","+ U","* R","* U") ) {
    seg1chr=contig$chrom2; seg1start=contig$start2-breakSlopSize+1; seg1end=contig$end2; seg1rc=FALSE
    seg2chr=contig$chrom1; seg2start=contig$start1; seg2end=contig$end1+breakSlopSize-1; seg2rc=FALSE
    contigSeqs = append(contigSeqs, formContig(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,paste("contig-",i,"(L-,R+)|",sep=""),contig$type,contig$support))
  }
  #generate (+,R1) x (+,R2) contig for compatible modes
  # [seg1] | R1+ => x <= rc(R2+) | [rc(seg2)]
  if(paste(contig$strand1,contig$direction1) %in% c("+ R","+ U","* R","* U")  
     && paste(contig$strand2,contig$direction2) %in% c("+ R","+ R","* R","* U") ) {
    seg1chr=contig$chrom1; seg1start=contig$start1-breakSlopSize+1; seg1end=contig$end1; seg1rc=FALSE
    seg2chr=contig$chrom2; seg2start=contig$start2; seg2end=contig$end2+breakSlopSize-1; seg2rc=TRUE
    contigSeqs = append(contigSeqs, formContig(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,paste("contig-",i,"(R1+,R2+)|",sep=""),contig$type,contig$support))
  }
  #generate (+,R2) x (+,R1) contig for compatible modes
  # [seg1] | R2+ => x <= rc(R1+) | [rc(seg2)]
  if(paste(contig$strand1,contig$direction1) %in% c("+ R","+ U","* R","* U")  
     && paste(contig$strand2,contig$direction2) %in% c("+ R","+ R","* R","* U") ) {
    seg1chr=contig$chrom2; seg1start=contig$start2-breakSlopSize+1; seg1end=contig$end2; seg1rc=FALSE
    seg2chr=contig$chrom1; seg2start=contig$start1; seg2end=contig$end1+breakSlopSize-1; seg2rc=TRUE
    contigSeqs = append(contigSeqs, formContig(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,paste("contig-",i,"(R2+,R1+)|",sep=""),contig$type,contig$support))
  }
  #generate (-,L1) x (-,L2) contig for compatible modes
  # [rc(seg1)] | rc(L1-) => x <= L2- | [seg2]
  if(paste(contig$strand1,contig$direction1) %in% c("- L","- U","* L","* U")  
     && paste(contig$strand2,contig$direction2) %in% c("- L","- U","* L","* U") ) {
    seg1chr=contig$chrom1; seg1start=contig$start1-breakSlopSize+1; seg1end=contig$end1; seg1rc=TRUE
    seg2chr=contig$chrom2; seg2start=contig$start2; seg2end=contig$end2+breakSlopSize-1; seg2rc=FALSE
    contigSeqs = append(contigSeqs, formContig(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,paste("contig-",i,"(L1-,L2-)|",sep=""),contig$type,contig$support))
  }
  #generate (-,L2) x (-,L1) contig for compatible modes
  # [rc(seg1)] | rc(L2-) => x <= L1- | [seg2]
  if(paste(contig$strand1,contig$direction1) %in% c("- L","- U","* L","* U")  
     && paste(contig$strand2,contig$direction2) %in% c("- L","- U","* L","* U") ) {
    seg1chr=contig$chrom2; seg1start=contig$start2-breakSlopSize+1; seg1end=contig$end2; seg1rc=TRUE
    seg2chr=contig$chrom1; seg2start=contig$start1; seg2end=contig$end1+breakSlopSize-1; seg2rc=FALSE
    contigSeqs = append(contigSeqs, formContig(seg1chr,seg1start,seg1end,seg1rc,seg2chr,seg2start,seg2end,seg2rc,paste("contig-",i,"(L2-,L1-)|",sep=""),contig$type,contig$support))
  }
  
  contigSeqs
}

#contigSeqs[[1]] = retrieveContig(1, contigs, breakSlopSize, padding, refSeq)

#need to work on this
contigSeqs = DNAStringSet(unlist(lapply(seq_along(contigs),"retrieveContig", contigs, breakSlopSize, padding, refSeq)))

#output to fasta
message("write contigs to:", paste(outPrefix,"fas",sep="."))
writeXStringSet(contigSeqs, file=paste(outPrefix,"fas",sep="."), format="fasta")

#quit()
#old proof-of-idea scripts
#IMbr1=min(info[info$supportType=="rdpIM",]$start1)
#IMbr2=min(info[info$supportType=="rdpIM",]$start2)
#RRbr1=min(info[info$supportType=="rdpRR",]$start1)
#RRbr2=min(info[info$supportType=="rdpRR",]$start2)
#FFbr1=max(info[info$supportType=="rdpFF",]$start1)
#FFbr2=max(info[info$supportType=="rdpFF",]$start2)

#only trust the start location of mapping
#hasSupportTypes=sapply(supportTypes, function(x) { any(grepl(x,info$supportType)) } )

