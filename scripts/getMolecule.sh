#!/bin/bash

if [ -z $6 ]; then
	echo "usage: getMolecule.sh inMoleFile inPrefix lineNum outDir suffix"
  echo "example:"
  echo "  getMolecule.sh NA12878.wg12.chr22.bam.dis.E40.refined NA12878.wg12.chr22.bam 1 PlotDir cisSize case/control"
  echo "  #this gets singlue molecules supports 1st junction listed in NA12878.wg12.chr22.bam.dis.E40.refined"
  echo "  #from NA12878.wg12.chr22.bam and store them in PlotDir"
  exit
fi 

echo "extracting molecules ..."
echo "inMolFile=$1"
inMolFile=$1
echo "inPrefix=$2"
inPrefix=$2
echo "lineNum=$3"
lineNum=$3
headLineNum=$(expr $lineNum + 1)
tailLineNum=1
outDir=$4
echo "outDir=$outDir"_"$lineNum"
outDir="$outDir"_"$lineNum"
cisSize=$5
echo "cisSize=$5"
suffix=$6
echo "suffix=$6"
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
echo "bedtools intersect -u -a $inPrefix.multiFragBar.bed.a0 -b $outDir/refined.reg1.bed >$outDir/refined.reg1.mol.$suffix"
bedtools intersect -u -a $inPrefix.multiFragBar.bed.a0 -b $outDir/refined.reg1.bed >$outDir/refined.reg1.mol.$suffix
echo "bedtools intersect -u -a $inPrefix.multiFragBar.bed.a0 -b $outDir/refined.reg2.bed >$outDir/refined.reg2.mol.$suffix"
bedtools intersect -u -a $inPrefix.multiFragBar.bed.a0 -b $outDir/refined.reg2.bed >$outDir/refined.reg2.mol.$suffix
###
### only produce these files for case suffix
### first eight columns are fixed now as specified in bedpe format
### 
if [[ $suffix = "case" ]]; then
  echo "#generating rdp for" $suffix
  FRcol=9
  RFcol=10
  FFcol=11
  RRcol=12
  echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$FRcol | grep -o '[0-9]\+' | sed '/^$/d' | sed -e 's/$/p/g' | tr '\n' ';' >$outDir/refined.rindex.FR.$suffix"
  head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$FRcol | grep -o '[0-9]\+' | sed '/^$/d' |  sed -e 's/$/p/g' | tr '\n' ';' >$outDir/refined.rindex.FR.$suffix
  echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$RFcol | grep -o '[0-9]\+' | sed '/^$/d' | sed -e 's/$/p/g' | tr '\n' ';' >$outDir/refined.rindex.RF.$suffix"
  head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$RFcol | grep -o '[0-9]\+' | sed '/^$/d' |  sed -e 's/$/p/g' | tr '\n' ';' >$outDir/refined.rindex.RF.$suffix
  echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$FFcol | grep -o '[0-9]\+' | sed '/^$/d' |  sed -e 's/$/p/g' | tr '\n' ';' >$outDir/refined.rindex.FF.$suffix"
  head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$FFcol | grep -o '[0-9]\+' | sed '/^$/d' |  sed -e 's/$/p/g' | tr '\n' ';' >$outDir/refined.rindex.FF.$suffix
  echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$RRcol | grep -o '[0-9]\+' | sed '/^$/d' |  sed -e 's/$/p/g' | tr '\n' >$outDir/refined.rindex.RR.$suffix"
  head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f$RRcol | grep -o '[0-9]\+' |  sed '/^$/d' | sed -e 's/$/p/g' | tr '\n' ';' >$outDir/refined.rindex.RR.$suffix
  ###
  echo "sed -n -f $outDir/refined.rindex.FR.$suffix $inPrefix.geq$cisSize.FR.bedpe >$outDir/refined.FR.bedpe.$suffix"
  sed -n -f $outDir/refined.rindex.FR.$suffix $inPrefix.geq$cisSize.FR.bedpe >$outDir/refined.FR.bedpe.$suffix
  echo "sed -n -f $outDir/refined.rindex.RF.$suffix $inPrefix.geq$cisSize.RF.bedpe >$outDir/refined.RF.bedpe.$suffix"
  sed -n -f $outDir/refined.rindex.RF.$suffix $inPrefix.geq$cisSize.RF.bedpe >$outDir/refined.RF.bedpe.$suffix
  echo "sed -n -f $outDir/refined.rindex.FF.$suffix $inPrefix.geq$cisSize.FF.bedpe >$outDir/refined.FF.bedpe.$suffix"
  sed -n -f $outDir/refined.rindex.FF.$suffix $inPrefix.geq$cisSize.FF.bedpe >$outDir/refined.FF.bedpe.$suffix
  echo "sed -n -f $outDir/refined.rindex.RR.$suffix $inPrefix.geq$cisSize.RR.bedpe >$outDir/refined.RR.bedpe.$suffix"
  sed -n -f $outDir/refined.rindex.RR.$suffix $inPrefix.geq$cisSize.RR.bedpe >$outDir/refined.RR.bedpe.$suffix
  ### find common barcodes
  (cut -f11 $outDir/refined.reg1.mol.$suffix | while read pat; do grep "$pat" $outDir/refined.reg2.mol.$suffix; done;) > $outDir/refined.jct.mol.$suffix
fi
echo "getMolecule.sh Done"
