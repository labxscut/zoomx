#!/bin/bash

if [ -z $2 ]; then
  echo "usage:"
  echo "  validateCall.sh chr:start-end x.refined/somatic"
  echo "  e.g. validateCall.sh 22:22505000-23275000 NA12878_WGS_phased_possorted_bam.bam.dis.E25.refined"
  exit
fi

reg=$1
zoomx=$2
dir="ref"

pacbio="$dir/NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam"
moleculo="$dir/NA12878.moleculo.bwa-mem.20140110.bam"
regArray=($(echo $reg | tr ":-" "\t"))
rstring=$(date +%s) #random string from date
echo "===validateRegion==="
mkdir ~/tmp
tmpfile="$HOME/tmp/zoomx_$rstring"
tmpfile2="$HOME/tmp/zoomx_2$rstring"
touch $tmpfile
touch $tmpfile2
echo "tmpfile=$tmpfile"
echo "tmpfile2=$tmpfile2"
echo $reg | tr ":-" "\t" >$tmpfile
echo chr$reg | tr ":-" "\t" >$tmpfile2
cat $tmpfile
cat $tmpfile2 #ucsc format

echo "===pacbioSplitRead==="
samtools view -o $pacbio.reg.bam -bh $pacbio $reg
bedtools bamtobed -bedpe -i $pacbio.reg.bam >$pacbio.reg.bam.bedpe
bedtools intersect -u -a $pacbio.reg.bam.bedpe -b $tmpfile

echo "===moleculoSplitRead==="
samtools view -o $moleculo.reg.bam -bh $moleculo $reg
bedtools bamtobed -bedpe -i $moleculo.reg.bam >$moleculo.reg.bam.bedpe
bedtools intersect -u -a $moleculo.reg.bam.bedpe -b $tmpfile

echo "===cloneByClone==="
file="$ref/NA12878_Kidd2008_hg19.bed"
bedtools intersect -u -a $file -b $tmpfile2

echo "===arrayCGH==="
file="$ref/NA12878_Conrad2010_hg19.bed"
bedtools intersect -u -a $file -b $tmpfile2

echo "===conventionalWGS==="
file="$ref/NA12878_Pacbio_Moleculo_Lumpy.bedpe"
bedtools intersect -u -a $file -b $tmpfile2

echo "===Hybrid==="
file="$ref/NA12878_1KG_SV.vcf.hg19.bed"
bedtools intersect -u -a $file -b $tmpfile2

#echo "===10xCalls==="
#file="NA12878_WGS_210_large_sv_calls.bedpe"
#bedtools intersect -u -a $file -b $tmpfile2
#
#echo "===10xCandidates==="
#file="NA12878_WGS_210_large_sv_candidates.bedpe"
#bedtools intersect -u -a $file -b $tmpfile2

