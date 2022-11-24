#!/bin/bash
## Assess replicate reproducibility ##

cat filenames | while read i; 
do
## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.
awk -v w=500 '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ./bed/${i}_mm10_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n > ./bed/${i}_mm10_bowtie2.fragmentsCount.bin500.bed &
done
