#!/bin/bash

cat filenames | while read i; 
do
## Extract the 9th column from the alignment sam file which is the fragment length
nohup samtools view \
-F 0x04 ./bam/${i}_mm10_bowtie2.sam | awk \
-F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk \
-v OFS="\t" '{print $2, $1/2}' > ./fragmentLen/${i}_fragmentLen.txt &

done
