#!/bin/bash

if [ -z $4 ]; then
	echo "usage: getMolecule.sh inMoleFile inPrefix lineNum outDir"
  echo "example:"
  echo "  getMolecule.sh P00526_2721_mt.st.10X_v2.bam.tra.E25.refined P00526_2721_mt.st.10X_v2.bam 1 testPlot"
  echo "  #this gets singlue molecules supports 1st junction listed in P00526_2721_mt.st.10X_v2.bam.tra.E25.refined"
  echo "  #from P00526_2721_mt.st.10X_v2.bam and store them in testPlot"
  exit
fi 

echo "extracting molecules ..."
echo "inMolFile=$1"
inMolFile=$1
echo "inPrefix=$2"
inPrefix=$2
echo "lineNum=$3"
let headLineNum=$3+1
let tailLineNum=1
echo "outDir=$4-$3"
outDir=$4-$3
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
echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f8 | grep -o '[0-9]*' | sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.IM"
head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f8 | grep -o '[0-9]*' | sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.IM
echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f9 | grep -o '[0-9]*' | sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.FF"
head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f9 | grep -o '[0-9]*' | sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.FF
echo "head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f10 | grep -o '[0-9]*' | sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.RR"
head -n $headLineNum $inMolFile | tail -n $tailLineNum | cut -f10 | grep -o '[0-9]*' | sed -e 's/$/p/g' | tr '\n' ';' | sed -e 's/p;$/p;\n/g' >$outDir/refined.rindex.RR
###
echo "sed -n -f $outDir/refined.rindex.IM $inPrefix.geq50K.IM.bedpe >$outDir/refined.IM.bedpe"
sed -n -f $outDir/refined.rindex.IM $inPrefix.geq50K.IM.bedpe >$outDir/refined.IM.bedpe
echo "sed -n -f $outDir/refined.rindex.FF $inPrefix.geq50K.FF.bedpe >$outDir/refined.FF.bedpe"
sed -n -f $outDir/refined.rindex.FF $inPrefix.geq50K.FF.bedpe >$outDir/refined.FF.bedpe
echo "sed -n -f $outDir/refined.rindex.RR $inPrefix.geq50K.RR.bedpe >$outDir/refined.RR.bedpe"
sed -n -f $outDir/refined.rindex.RR $inPrefix.geq50K.RR.bedpe >$outDir/refined.RR.bedpe
echo "getMolecule.sh Done"
