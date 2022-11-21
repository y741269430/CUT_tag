#!/bin/bash
## mapping ecoil (bowtie2) ##

spikeInRef="/home/yangjiajun/downloads/genome/ecoil_U00096.3/bowtie2_idx/ecoil"

cat filenames | while read i; 
do
nohup bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 \
-I 10 -X 700 -p 8 -x ${spikeInRef} \
-1 ./raw/BL6-TG-CUT-${i}_1.fq.gz \
-2 ./raw/BL6-TG-CUT-${i}_2.fq.gz \
-S ./bam/${i}_ecoil_bowtie2.sam &> ./bowtie2_summary/${i}_ecoil_bowtie2.txt &
done
