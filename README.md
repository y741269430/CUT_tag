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

ls raw/*1.fq.gz |cut -d "_" -f 1 |cut -d "/" -f 2 |cut -d "-" -f 4-5 > filenames

## Alignment to mm10
bash cut1_mm10_bw2.sh

## Alignment to ecoil
bash cut1_ecoil_bw2.sh

## Remove duplicates





