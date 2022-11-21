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

    cd /home/yangjiajun/downloads/genome/mm10_GRCm38/ucsc_fa/

    samtools faidx GRCm38.primary_assembly.genome.fa &
    cut -f1,2 GRCm38.primary_assembly.genome.fa.fai > mm10.chrom.sizes &
    
    nohup bowtie2-build GRCm38.primary_assembly.genome.fa \
    /home/yangjiajun/downloads/genome/mm10_GRCm38/bowtie2_idx/mm10 &  


    cd /home/yangjiajun/downloads/genome/ecoil_U00096.3/ucsc_fa/

    samtools faidx GCF_000005845.2_ASM584v2_genomic.fa &
    cut -f1,2 GCF_000005845.2_ASM584v2_genomic.fa.fai > ecoil.chrom.sizes &
    
    nohup bowtie2-build GCF_000005845.2_ASM584v2_genomic.fa \
    /home/yangjiajun/downloads/genome/ecoil_U00096.3/bowtie2_idx/ecoil &   

----

## 1. Activate the source and create the folder  
    
    conda activate cuttag  
    mkdir bam bowtie2_summary picard_summary bedgraph SEACR
    
## 2. Write the filenames  
    
    ls raw/*1.fq.gz |cut -d "_" -f 1 |cut -d "/" -f 2 |cut -d "-" -f 4-5 > filenames

## 3.1.1 Alignment to mm10  
    
    bash cut1_mm10_bw2.sh

## 3.1.2 Alignment to spike-in genome for spike-in calibration (ecoil)  
    
    bash cut1_ecoil_bw2.sh

## 3.3 Remove duplicates  
    
    bash cut2_picard.sh

## 4.2 File format conversion  
    
    bash cut3_bam2bed.sh

## 4.3 Assess replicate reproducibility  
    
    bash cut4.sh

## 5.1 Scaling factor   

    samtools view -F 0x04 ./bam/CFA3-1_ecoil_bowtie2.sam | wc -l > ./bowtie2_summary/CFA3-1_ecoil_bowtie2.seqDepthDouble &
    samtools view -F 0x04 ./bam/CFA3-2_ecoil_bowtie2.sam | wc -l > ./bowtie2_summary/CFA3-2_ecoil_bowtie2.seqDepthDouble &
    samtools view -F 0x04 ./bam/CFA3-B1_ecoil_bowtie2.sam | wc -l > ./bowtie2_summary/CFA3-B1_ecoil_bowtie2.seqDepthDouble &
    samtools view -F 0x04 ./bam/NADCFA3-1_ecoil_bowtie2.sam | wc -l > ./bowtie2_summary/NADCFA3-1_ecoil_bowtie2.seqDepthDouble &
    samtools view -F 0x04 ./bam/NADCFA3-2_ecoil_bowtie2.sam | wc -l > ./bowtie2_summary/NADCFA3-2_ecoil_bowtie2.seqDepthDouble &

    seqDepth=$((seqDepthDouble/2))
    scale_factor=`echo "10000 / $seqDepth" | bc -l`
    echo $scale_factor

    bedtools genomecov \
    -bg -scale $scale_factor \
    -i ./bam/${i}_mm10_bowtie2.fragments.bed \
    -g /home/yangjiajun/downloads/genome/mm10_GRCm38/ucsc_fa/mm10.chrom.sizes > ./bedgraph/${i}_mm10_bowtie2.fragments.normalized.bedgraph


## 6.1. SEACR  

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/CFA3-1_mm10_bowtie2.fragments.normalized.bedgraph \
         ./bedgraph/CFA3-B1_mm10_bowtie2.fragments.normalized.bedgraph \
         non stringent ./SEACR/CFA3-1_seacr_control.peaks &

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/CFA3-2_mm10_bowtie2.fragments.normalized.bedgraph \
         ./bedgraph/CFA3-B1_mm10_bowtie2.fragments.normalized.bedgraph \
         non stringent ./SEACR/CFA3-2_seacr_control.peaks &

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/NADCFA3-1_mm10_bowtie2.fragments.normalized.bedgraph \
         ./bedgraph/CFA3-B1_mm10_bowtie2.fragments.normalized.bedgraph \
         non stringent ./SEACR/NADCFA3-1_seacr_control.peaks &

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/NADCFA3-2_mm10_bowtie2.fragments.normalized.bedgraph \
         ./bedgraph/CFA3-B1_mm10_bowtie2.fragments.normalized.bedgraph \
         non stringent ./SEACR/NADCFA3-2_seacr_control.peaks &



