#!/bin/bash

if [ -z $5 ]; then
	echo "usage: getMolecule.sh inMoleFile inPrefix lineNum outDir"
  echo "example:"
  echo "  getMolecule.sh NA12878.wg12.chr22.bam.dis.E40.refined NA12878.wg12.chr22.bam 1 PlotDir cisSize"
  echo "  #this gets singlue molecules supports 1st junction listed in NA12878.wg12.chr22.bam.dis.E40.refined"
  echo "  #from NA12878.wg12.chr22.bam and store them in PlotDir"
  exit
fi 

echo "extracting molecules ..."
echo "inMolFile=$1"
inMolFile=$1
echo "inPrefix=$2"
inPrefix=$2
source $inPrefix
inPrefix=$datDir/
echo "lineNum=$3"
lineNum=$3
headLineNum=$(expr $lineNum + 1)
tailLineNum=1
outDir=$4
echo "outDir=$outDir"_"$lineNum"
outDir="$outDir"_"$lineNum"
cisSize=$5
echo "cisSize=$5"
#inMolFile=results/NA12878.wg12.chr22-22M-24M.bam.dist.E45.refined
#inPrefix=data/NA12878.wg12.chr22-22M-24M.bam
#lineNum=1
#outDir=plots/NA12878.wg12.chr22-22M-24M.bam.dist.E45.refined.Plot_1
#cisSize=50000
mkdir $outDir
echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum >$outDir/refined.bedpe"
head -n $headLineNum $inMolFile | tail -n $tailLineNum >$outDir/refined.bedpe
echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f1,2,3 >$outDir/refined.reg1.bed"
head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f1,2,3 >$outDir/refined.reg1.bed
echo -n "found reg1="
cat $outDir/refined.reg1.bed
echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f4,5,6 >$outDir/refined.reg2.bed"
head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f4,5,6 >$outDir/refined.reg2.bed
echo -n "found reg2="
cat $outDir/refined.reg2.bed
echo "bedtools intersect -u -a $inPrefix.multiFragBar.bed.a0 -b $outDir/refined.reg1.bed >$outDir/refined.reg1.mol"
bedtools intersect -u -a $inPrefix.multiFragBar.bed.a0 -b $outDir/refined.reg1.bed >$outDir/refined.reg1.mol
echo "bedtools intersect -u -a $inPrefix.multiFragBar.bed.a0 -b $outDir/refined.reg2.bed >$outDir/refined.reg2.mol"
bedtools intersect -u -a $inPrefix.multiFragBar.bed.a0 -b $outDir/refined.reg2.bed >$outDir/refined.reg2.mol
###
###
FRcol=8
RFcol=9
FFcol=10
RRcol=11
echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$FRcol | grep -o '[0-9]\+' | sed '/^$/d' | sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.FR"
head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$FRcol | grep -o '[0-9]\+' | sed '/^$/d' |  sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.FR
echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$RFcol | grep -o '[0-9]\+' | sed '/^$/d' | sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.RF"
head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$RFcol | grep -o '[0-9]\+' | sed '/^$/d' |  sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.RF
echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$FFcol | grep -o '[0-9]\+' | sed '/^$/d' |  sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.FF"
head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$FFcol | grep -o '[0-9]\+' | sed '/^$/d' |  sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.FF
echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$RRcol | grep -o '[0-9]\+' | sed '/^$/d' |  sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.RR"
head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$RRcol | grep -o '[0-9]\+' |  sed '/^$/d' | sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.RR
###
echo "sed -n -f $outDir/refined.rindex.FR $inPrefix.geq$cisSize.FR.bedpe >$outDir/refined.FR.bedpe"
sed -n -f $outDir/refined.rindex.FR $inPrefix.geq$cisSize.FR.bedpe >$outDir/refined.FR.bedpe
echo "sed -n -f $outDir/refined.rindex.RF $inPrefix.geq$cisSize.RF.bedpe >$outDir/refined.RF.bedpe"
sed -n -f $outDir/refined.rindex.RF $inPrefix.geq$cisSize.RF.bedpe >$outDir/refined.RF.bedpe
echo "sed -n -f $outDir/refined.rindex.FF $inPrefix.geq$cisSize.FF.bedpe >$outDir/refined.FF.bedpe"
sed -n -f $outDir/refined.rindex.FF $inPrefix.geq$cisSize.FF.bedpe >$outDir/refined.FF.bedpe
echo "sed -n -f $outDir/refined.rindex.RR $inPrefix.geq$cisSize.RR.bedpe >$outDir/refined.RR.bedpe"
sed -n -f $outDir/refined.rindex.RR $inPrefix.geq$cisSize.RR.bedpe >$outDir/refined.RR.bedpe
echo "getMolecule.sh Done"
