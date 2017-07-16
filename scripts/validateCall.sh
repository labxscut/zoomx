#!/bin/bash

if [ -z $2 ]; then
  (>&2 echo "usage:")
  (>&2 echo "  validateCall.sh chr:start-end")
  (>&2 echo "  e.g. validateCall.sh chr22:22525000-23260000 call.bedpe")
  exit
fi
#input will be chr

bedpe=$2
reg=`echo $1 | sed s/chr//g`
refdir="ref"
restxt=""

pacbio="NA12878.pacbio_fr_MountSinai.bwa-sw.20140211.bam"
moleculo="NA12878.moleculo.bwa-mem.20140110.bam"
schadt="schadt_NA12878_pacbio.bam"
regArray=($(echo $reg | tr ":-" "\t"))
rstring=$(date +%s) #random string from date

(>&2 echo "===validateRegion===")
(>&2 echo $reg)
mkdir ~/tmp
tmpfile1=$HOME/tmp/$bedpe.tmp1
tmpfile2=$HOME/tmp/$bedpe.tmp2
tmpfile3=$HOME/tmp/$bedpe.tmp3
tmpfile4=$HOME/tmp/$bedpe.tmp4
(>&2 echo "tmpfile1=$tmpfile1")
(>&2 echo "tmpfile2=$tmpfile2")
(>&2 echo "tmpfile3=$tmpfile3")
(>&2 echo "tmpfile4=$tmpfile4")
cat $bedpe | sed 's/chr//g' >$tmpfile1
cat $bedpe >$tmpfile2                                 #ucsc
echo $reg | tr ":-" "\t" >$tmpfile3
echo chr$reg | tr ":-" "\t" >$tmpfile4                #ucsc
restxt=$reg
(>&2 cat $tmpfile1)
(>&2 cat $tmpfile2)
(>&2 cat $tmpfile3)
(>&2 cat $tmpfile4)

#svReg pacbioSplit moleculoSplit Kidd2008 Conrad2010 Layer2014 Sudmant2016

(>&2 echo "===pacbioSplitRead===")
samtools view -o $HOME/tmp/$pacbio.reg.bam -bh $refdir/$pacbio $reg
bedtools bamtobed -bedpe -i $HOME/tmp/$pacbio.reg.bam >$HOME/tmp/$pacbio.reg.bam.bedpe
full=`bedtools pairtopair -type both -is -a $HOME/tmp/$pacbio.reg.bam.bedpe -b $tmpfile1`
cnt=`bedtools pairtopair -type both -is -a $HOME/tmp/$pacbio.reg.bam.bedpe -b $tmpfile1 | wc -l`
restxt="$restxt\t$cnt"
(>&2 echo $full)

(>&2 echo "===schadtSplitRead===")
samtools view -o $HOME/tmp/$schadt.reg.bam -bh $refdir/$schadt chr$reg
bedtools bamtobed -bedpe -i $HOME/tmp/$schadt.reg.bam >$HOME/tmp/$schadt.reg.bam.bedpe
full=`bedtools pairtopair -type both -is -a $HOME/tmp/$schadt.reg.bam.bedpe -b $tmpfile2`
cnt=`bedtools pairtopair -type both -is -a $HOME/tmp/$schadt.reg.bam.bedpe -b $tmpfile2 | wc -l`
restxt="$restxt\t$cnt"
(>&2 echo $full)

(>&2 echo "===moleculoSplitRead===")
samtools view -o $HOME/tmp/$moleculo.reg.bam -bh $refdir/$moleculo $reg
bedtools bamtobed -bedpe -i $HOME/tmp/$moleculo.reg.bam >$HOME/tmp/$moleculo.reg.bam.bedpe
full=`bedtools pairtopair -type both -is -a $HOME/tmp/$moleculo.reg.bam.bedpe -b $tmpfile1`
cnt=`bedtools pairtopair -type both -is -a $HOME/tmp/$moleculo.reg.bam.bedpe -b $tmpfile1 | wc -l`
restxt="$restxt\t$cnt"
(>&2 echo $full)

(>&2 echo "===cloneByClone===")
file="$refdir/NA12878_Kidd2008_hg19.bed"
full=`bedtools intersect -u -a $file -b $tmpfile4`
cnt=`bedtools intersect -u -a $file -b $tmpfile4 | wc -l`
restxt="$restxt\t$cnt"
(>&2 echo $full)

(>&2 echo "===arrayCGH===")
file="$refdir/NA12878_Conrad2010_hg19.bed"
full=`bedtools intersect -u -a $file -b $tmpfile4`
cnt=`bedtools intersect -u -a $file -b $tmpfile4 | wc -l`
restxt="$restxt\t$cnt"
(>&2 echo $full)

(>&2 echo "===conventionalWGS===")
file="$refdir/NA12878_Pacbio_Moleculo_Lumpy.bedpe"
full=`bedtools intersect -u -a $file -b $tmpfile4`
cnt=`bedtools intersect -u -a $file -b $tmpfile4 | wc -l`
restxt="$restxt\t$cnt"
(>&2 echo $full)

(>&2 echo "===layerPB===")
file="$refdir/NA12878_Pacbio_Moleculo_Lumpy.bedpe"
cnt=`bedtools intersect -u -a $file -b $tmpfile4 | grep -o "PB,[0-9]\+"`
restxt="$restxt\t$cnt"
(>&2 echo $cnt)

(>&2 echo "===layerMC===")
file="$refdir/NA12878_Pacbio_Moleculo_Lumpy.bedpe"
cnt=`bedtools intersect -u -a $file -b $tmpfile4 | grep -o "MC,[0-9]\+"`
restxt="$restxt\t$cnt"
(>&2 echo $cnt)

(>&2 echo "===Hybrid===")
file="$refdir/NA12878_1KG_SV.vcf.hg19.bed"
full=`bedtools intersect -u -a $file -b $tmpfile4`
cnt=`bedtools intersect -u -a $file -b $tmpfile4 | wc -l`
restxt="$restxt\t$cnt"
(>&2 echo $full)

echo -e $restxt

#echo "===10xCalls==="
#file="NA12878_WGS_210_large_sv_calls.bedpe"
#bedtools intersect -u -a $file -b $tmpfile2
#
#echo "===10xCandidates==="
#file="NA12878_WGS_210_large_sv_candidates.bedpe"
#bedtools intersect -u -a $file -b $tmpfile2

