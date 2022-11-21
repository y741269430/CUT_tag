#!/bin/bash
## mapping mm10 (bowtie2) ##

mm10="/home/yangjiajun/downloads/genome/mm10_GRCm38/bowtie2_idx/mm10"

cat filenames | while read i; 
do
nohup bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 \
-I 10 -X 700 -p 8 -x ${mm10} \
-1 ./raw/BL6-TG-CUT-${i}_1.fq.gz \
-2 ./raw/BL6-TG-CUT-${i}_2.fq.gz \
-S ./bam/${i}_mm10_bowtie2.sam &> ./bowtie2_summary/${i}_mm10_bowtie2.txt &
done
