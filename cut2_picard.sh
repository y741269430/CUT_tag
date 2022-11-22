#!/bin/bash
## Remove duplicates ##

cat filenames | while read i; 
do
## Sort by coordinate
nohup picard SortSam -I ./bam/${i}_mm10_bowtie2.sam \
-O ./bam/${i}_bowtie2.sorted.sam \
-SO coordinate &&

## mark duplicates
picard MarkDuplicates -I ./bam/${i}_bowtie2.sorted.sam \
-O ./bam/${i}_bowtie2.sorted.dupMarked.sam \
-M ./picard_summary/${i}_picard.dupMark.txt &&

## remove duplicates
picard MarkDuplicates -I ./bam/${i}_bowtie2.sorted.sam \
-O ./bam/${i}_bowtie2.sorted.rmDup.sam \
--REMOVE_DUPLICATES true \
-M ./picard_summary/${i}_picard.rmDup.txt &
done
