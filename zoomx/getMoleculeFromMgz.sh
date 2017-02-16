#!/bin/bash

if [ -z $3 ]; then
	echo "usage: getMolecule.sh Mgz line bam ref"
  echo "example:"
  echo "  getMolecule.sh NA12878_WGS_phased_possorted_bam.bam.a0.grid13K.M_data.jct.gz 2 NA12878_WGS_phased_possorted_bam.bam /home/lixia/ws/hg/hg19/human_g1k_v37.ucsc.fasta"
  echo "  Mgz=NA12878_WGS_phased_possorted_bam.bam.a0.grid13K.M_data.jct.gz" 
  echo "  line=2"
  echo "  bam=NA12878_WGS_phased_possorted_bam.bam"
  echo "  ref=/home/lixia/ws/hg/hg19/human_g1k_v37.ucsc.fasta"
  exit
fi 

Mgz=$1
line=$2
bam=$3
ref=$4
bed=$bam.multiFragBar.bed.a0
h5=$bam.multiFragBar.h5.a0
gap=$ref.acen+gap+hc.bed
grid=$ref.bed.w10Ks5K
out=$Mgz.$line
Xgz=$bam.a0.grid13K.X_data.gz
Tgz=$bam.a0.grid13K.tagcol2tagidx.gz

GridFile=$grid
BamFile=$bam
echo GridFile=$GridFile
echo BamFile=$BamFile
MDataFile=$Mgz
echo MDataFile=$MDataFile
XDataFile=$Xgz
echo XDataFile=$XDataFile
MoleculeFile=$bed
echo MoleculeFile=$MoleculeFile
Tagcol2TagidxFile=$Tgz
echo Tagcol2TagidxFile=$Tagcol2TagidxFile

#find all tags shared by grid1 and grid2
base0grid1=`zcat $MDataFile | head -n $line | tail -n 1 | cut -f1`
base0grid2=`zcat $MDataFile | head -n $line | tail -n 1 | cut -f2`
let base1grid1=$base0grid1+1
let base1grid2=$base0grid2+1
cat <(zcat $XDataFile | grep "^$base1grid1[.0]*\b") <(zcat $XDataFile | grep "^$base1grid2[.0]*\b") | cut -f2 | sort -n | uniq -d | sed 's/\.0//g' >$out.tagcol

#output 
sedcmd=$(cat $out.tagcol | awk '{ print $1+1 }' | tr '\n' 'p' | sed 's/p/p;/g') 
sed -n $sedcmd <(zcat $Tagcol2TagidxFile) | awk '{ print $1+1 }' >$out.tagidx
awk 'NR==FNR{a[$1]=++i;next} { if ($4 in a) {print $0}}' $out.tagidx $MoleculeFile >$out.bed
echo "find result in $out.bed"
echo "done getting the molecules having overlapping tag/barcodes" 

#echo GridIndex0File=$1
#GridIndex0Tag=${GridIndex0File%.gi}
#echo GridIndex0Tag=$GridIndex0Tag
#sedcmd=$(cat $out.tagcol | awk '{ print $1+1 }' | tr '\n' 'p' | sed 's/p/p;/g') 
#sed -n $sedcmd <(zcat $Tagcol2TagidxFile) | awk '{ print $1+1 }' >$BamFile.a$atag.$GridIndex0Tag.overlap.tagidx_1b
#sed -n $sedcmd $GridFile >$BamFile.multiFragBar.a$atag.$GridIndex0Tag.overlap.bed
#echo "done getting involved grid regions"
#GridIndex01=(`cat $GridIndex0File`)
#GirdIndex0=(`cat $GridIndex0File`)
#for ind in ${GridIndex01[@]}; do
#  (( ind++ ))
#awk 'NR==FNR{a[$1]=++i;next} { if ($1 in a) {print $2}}' $GridIndex0File <(zcat $XDataFile) | sort -n | uniq -c | sort -n -k 1 -r | awk '{ print $2+1 }' | sort -n >$BamFile.a$atag.$GridIndex0Tag.overlap.tagcol_1b
#echo "done getting tagcols for molecules overlapping grids" 
#sedcmd=$(cat $BamFile.a$atag.$GridIndex0Tag.overlap.tagcol_1b | tr '\n' 'p' | sed 's/p/p;/g') 
#sed -n $sedcmd <(zcat $Tagcol2TagidxFile) | awk '{ print $1+1 }' >$BamFile.a$atag.$GridIndex0Tag.overlap.tagidx_1b
#echo "done convertting tagcols to tagidx for molecules overlapping grids" 
#awk 'NR==FNR{a[$1]=++i;next} { if ($4 in a) {print $0}}' $BamFile.a$atag.$GridIndex0Tag.overlap.tagidx_1b $MoleculeFile | bedtools intersect -u -a - -b $BamFile.multiFragBar.a$atag.$GridIndex0Tag.overlap.bed >$BamFile.multiFragBar.a$atag.$GridIndex0Tag.overlap.molecules
#echo "done getting the molecules overlapping grids" 
