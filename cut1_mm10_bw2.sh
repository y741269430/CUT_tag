#!/bin/bash
## make id.config ##
## mapping (bowtie2) ##

ls raw/*1.fq.gz > 1
ls raw/*2.fq.gz > 2
ls raw/*1.fq.gz |cut -d "_" -f 1 |cut -d "/" -f 2 |cut -d "-" -f 4-5 > 0
paste 0 1 2 > id.config

cat id.config |while read id; 
do echo $id
arr=($id)
id2=${arr[2]}
id1=${arr[1]}
sample=${arr[0]}

mm10="/home/yangjiajun/downloads/genome/mm10_GRCm38/bowtie2_idx/mm10"

# paired end
nohup bowtie2 \
--end-to-end \
--very-sensitive \
--no-mixed \
--no-discordant \
--phred33 \
-I 10 \
-X 700 \
-p 8 \
-x ${mm10} \
-1 $id1 \
-2 $id2 \
-S ./bam/${sample}_mm10_bowtie2.sam &> ./info/${sample}_mm10_bowtie2.txt &
done

rm 0 1 2 id.config
