#!/bin/bash
## File format conversion ##

cat filenames | while read i; 
do
## Filter and keep the mapped read pairs
nohup samtools view -bS -F 0x04 ./bam/${i}_mm10_bowtie2.sam > ./bam/${i}_mm10_bowtie2.mapped.bam &&

## Convert into bed file format
bedtools bamtobed -i ./bam/${i}_mm10_bowtie2.mapped.bam -bedpe > ./bed/${i}_mm10_bowtie2.bed &&

## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
awk '$1==$4 && $6-$2 < 1000 {print $0}' ./bed/${i}_mm10_bowtie2.bed > ./bed/${i}_mm10_bowtie2.clean.bed &&

## Only extract the fragment related columns
cut -f 1,2,6 ./bed/${i}_mm10_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n > ./bed/${i}_mm10_bowtie2.fragments.bed &
done
