# copy grch38 fasta file into folder

#starting with hg38 files

bedtools slop -g human_g1k_v37.ucsc.fasta.size -i human_g1k_v37.ucsc.fasta.gap.bed -b 500000 | bedtools merge -i - >human_g1k_v37.ucsc.fasta.gap.bed.slop500K
bedtools slop -g human_g1k_v37.fasta.size -i human_g1k_v37.fasta.gap.bed -b 500000 | bedtools merge -i - >human_g1k_v37.fasta.gap.bed.slop500K

samtools faidx hg19/human_g1k_v37.ucsc.fasta
samtools faidx hg19/human_g1k_v37.fasta

dir=~/ws/hg/hg19
#sentieon grids
bedtools makewindows -w 99 -s 100 -b $dir/human_g1k_v37.ucsc.fasta.bed >$dir/human_g1k_v37.ucsc.fasta.bed.w99s100.sentieon

dir=~/ws/hg/hg19x

#non-overlapping windows
bedtools makewindows -w 100 -s 100 -b $dir/human_g1k_v37.fasta.bed >$dir/human_g1k_v37.fasta.bed.w100s100
bedtools makewindows -w 100 -s 100 -b $dir/human_g1k_v37.ucsc.fasta.bed >$dir/human_g1k_v37.ucsc.fasta.bed.w100s100
bedtools makewindows -w 10000 -s 10000 -b $dir/human_g1k_v37.fasta.bed >$dir/human_g1k_v37.fasta.bed.w10Ks10K
bedtools makewindows -w 10000 -s 10000 -b $dir/human_g1k_v37.ucsc.fasta.bed >$dir/human_g1k_v37.ucsc.fasta.bed.w10Ks10K
#overlapping windows
bedtools makewindows -w 10000 -s 5000 -b $dir/human_g1k_v37.fasta.bed >$dir/human_g1k_v37.fasta.bed.w10Ks5K
bedtools makewindows -w 10000 -s 5000 -b $dir/human_g1k_v37.ucsc.fasta.bed >$dir/human_g1k_v37.ucsc.fasta.bed.w10Ks5K
bedtools makewindows -w 100000 -s 50000 -b $dir/human_g1k_v37.fasta.bed >$dir/human_g1k_v37.fasta.bed.w100Ks50K
bedtools makewindows -w 100000 -s 50000 -b $dir/human_g1k_v37.ucsc.fasta.bed >$dir/human_g1k_v37.ucsc.fasta.bed.w100Ks50K
bedtools makewindows -w 1000000 -s 500000 -b $dir/human_g1k_v37.fasta.bed >$dir/human_g1k_v37.fasta.bed.w1Ms500K
bedtools makewindows -w 1000000 -s 500000 -b $dir/human_g1k_v37.ucsc.fasta.bed >$dir/human_g1k_v37.ucsc.fasta.bed.w1Ms500K
bedtools makewindows -w 10000000 -s 5000000 -b $dir/human_g1k_v37.fasta.bed >$dir/human_g1k_v37.fasta.bed.w10Ms5M
bedtools makewindows -w 10000000 -s 5000000 -b $dir/human_g1k_v37.ucsc.fasta.bed >$dir/human_g1k_v37.ucsc.fasta.bed.w10Ms5M
bedtools makewindows -w 100000000 -s 50000000 -b $dir/human_g1k_v37.fasta.bed >$dir/human_g1k_v37.fasta.bed.w100Ms50M
bedtools makewindows -w 100000000 -s 50000000 -b $dir/human_g1k_v37.ucsc.fasta.bed >$dir/human_g1k_v37.ucsc.fasta.bed.w100Ms50M
