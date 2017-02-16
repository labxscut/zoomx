file="P04979_CR1.wg10.bam"
frag="$file.multiFragBar.bed.a50000"
grid="$HOME/ws/hg/hg19/human_g1k_v37.fasta.bed.w100"
bedtools intersect -a $grid -b $frag -loj | grep -v [-][1] | cut -f 1,2,3,7 | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$1}'
