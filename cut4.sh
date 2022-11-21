#!/bin/bash
## Assess replicate reproducibility ##

binLen=500

cat filenames | while read i; 
do
## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ./bam/${i}_mm10_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n > ./bam/${i}_mm10_bowtie2.fragmentsCount.bin$binLen.bed &
done
