#!/bin/bash
## Spike-in calibration ##

cat filenames | while read i; 
do
bedtools genomecov -bg -i ./bed/${i}_mm10_bowtie2.fragments.bed \
-g /home/yangjiajun/downloads/genome/mm10_GRCm38/ucsc_fa/mm10.chrom.sizes > ./bedgraph/${i}_mm10_bowtie2.fragments.normalized.bedgraph &
done
