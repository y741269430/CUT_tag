# CUT_tag

## build envs  
  conda create -n cuttag  
  conda activate cuttag  

  conda install -c bioconda picard  
  conda install -c bioconda samtools
  conda install -c bioconda bedtools
  conda install -c bioconda deeptools
  conda install -c bioconda seacr
  conda install -c bioconda bowtie2
  conda install -c bioconda fastqc

## Build the bowtie2 reference genome index (mm10 and ecoil)  

Execute the following command once to generate a permanently used index!  

    nohup bowtie2-build /home/yangjiajun/downloads/genome/mm10_GRCm38/ucsc_fa/GRCm38.primary_assembly.genome.fa \
    /home/yangjiajun/downloads/genome/mm10_GRCm38/bowtie2_idx/mm10 &  
    nohup bowtie2-build /home/yangjiajun/downloads/genome/ecoil_U00096.3/ucsc_fa/GCF_000005845.2_ASM584v2_genomic.fa \
    /home/yangjiajun/downloads/genome/ecoil_U00096.3/bowtie2_idx/ecoil &

mkdir bam info

vim cut1_ecoil_bw2.sh

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

spikeInRef="/home/yangjiajun/downloads/genome/ecoil_U00096.3/bowtie2_idx/ecoil"

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
-x ${spikeInRef} \
-1 $id1 \
-2 $id2 \
-S ./bam/${sample}_ecoil_bowtie2.sam &> ./info/${sample}_ecoil_bowtie2.txt &
done

rm 0 1 2 id.config








