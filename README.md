# CUT&tag 
The pipeline were based on https://yezhengstat.github.io/CUTTag_tutorial/index.html  

## 0. Build source used for CUT&tag  
    
    conda create -n cuttag  
    conda activate cuttag  

    conda install -c bioconda picard  
    conda install -c bioconda samtools
    conda install -c bioconda bedtools
    conda install -c bioconda deeptools
    conda install -c bioconda seacr
    conda install -c bioconda bowtie2
    conda install -c bioconda fastqc

## 0. Build the bowtie2 reference genome index (mm10 and ecoil)  

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
    
    mkdir bam bed bowtie2_summary picard_summary bedgraph SEACR fragmentLen plot  
    
## 2. Write the filenames  
    
    ls raw/*1.fq.gz |cut -d "_" -f 1 |cut -d "/" -f 2 |cut -d "-" -f 4-5 > filenames

## 3.1.1 Alignment to mm10  
    
    vim cut1_mm10_bw2.sh
    
    #!/bin/bash
    ## Alignment to mm10 ##

    mm10="/home/yangjiajun/downloads/genome/mm10_GRCm38/bowtie2_idx/mm10"

    cat filenames | while read i; 
    do
    nohup bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 \
    -I 10 -X 700 -p 8 -x ${mm10} \
    -1 ./raw/BL6-TG-CUT-${i}_1.fq.gz \
    -2 ./raw/BL6-TG-CUT-${i}_2.fq.gz \
    -S ./bam/${i}_mm10_bowtie2.sam &> ./bowtie2_summary/${i}_mm10_bowtie2.txt &
    done

## 3.1.2 Alignment to spike-in genome for spike-in calibration (ecoil)   

E. coli DNA is carried along with bacterially-produced pA-Tn5 protein and gets tagmented non-specifically during the reaction. The fraction of total reads that map to the E.coli genome depends on the yield of epitope-targeted CUT&Tag, and so depends on the number of cells used and the abundance of that epitope in chromatin. Since a constant amount of pATn5 is added to CUT&Tag reactions and brings along a fixed amount of E. coli DNA, E. coli reads can be used to normalize epitope abundance in a set of experiments.  

    vim cut1_ecoil_bw2.sh
    
    #!/bin/bash
    ## Alignment to spike-in genome for spike-in calibration ##

    spikeInRef="/home/yangjiajun/downloads/genome/ecoil_U00096.3/bowtie2_idx/ecoil"

    cat filenames | while read i; 
    do
    nohup bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 \
    -I 10 -X 700 -p 8 -x ${spikeInRef} \
    -1 ./raw/BL6-TG-CUT-${i}_1.fq.gz \
    -2 ./raw/BL6-TG-CUT-${i}_2.fq.gz \
    -S ./bam/${i}_ecoil_bowtie2.sam &> ./bowtie2_summary/${i}_ecoil_bowtie2.txt &
    done

## 3.2 Report sequencing mapping summary  

Summarize the raw reads and uniquely mapping reads to report the efficiency of alignment. Alignment frequencies are expected to be >80% for high-quality data. CUT&Tag data typically has very low backgrounds, so as few as 1 million mapped fragments can give robust profiles for a histone modification in the human genome. Profiling of less-abundant transcription factors and chromatin proteins may require 10 times as many mapped fragments for downstream analysis.  

## 3.2.1 Sequencing depth  
R  

## 3.2.2 Spike-in alignment  
R  

## 3.2.3 Summarize the alignment to mm10 and E.coli  
R  

## 3.2.4 Visualizing the sequencing depth and alignment results.  
R  

In a typical CUT&Tag experiment targeting the abundant H3K27me3 histone modification in 65,000 K562 cells, the percentage of E. coli reads range from ~0.01% to 10%. With fewer cells or less abundant epitopes, E. coli reads can comprise as much as 70% or the total mapped reads. For IgG controls, the percentage of E. coli reads is typically much higher than that for an abundant histone modification.   

## 3.3 Remove duplicates   

CUT&Tag integrates adapters into DNA in the vicinity of the antibody-tethered pA-Tn5, and the exact sites of integration are affected by the accessibility of surrounding DNA. For this reason fragments that share exact starting and ending positions are expected to be common, and such ‘duplicates’ may not be due to duplication during PCR. In practice, we have found that the apparent duplication rate is low for high quality CUT&Tag datasets, and even the apparent ‘duplicate’ fragments are likely to be true fragments. Thus, we do not recommend removing the duplicates. In experiments with very small amounts of material or where PCR duplication is suspected, duplicates can be removed. The following commands show how to check the duplication rate using Picard.    

    vim cut2_picard.sh
    
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

R  

In these example datasets, the IgG control samples have relatively high duplication rates, since reads in this sample derive from non-specific tagmentation in the CUT&Tag reactions. Therefore, it is appropriate to remove the duplicates from the IgG datasets before downstream analysis.  

The estimated library size are the estimated number of unique molecules in the library based on PE duplication calculated by Picard.  

The estimated library sizes is proportional to the abundance of the targeted epitope and to the quality of the antibody used, while the estimated library sizes of IgG samples are expected to be very low.  

Unique fragment number is calculated by the MappedFragNum_hg38 * (1-DuplicationRate/100).  

## 3.4. Assess mapped fragment size distribution  

CUT&Tag inserts adapters on either side of chromatin particles in the vicinity of the tethered enzyme, although tagmentation within chromatin particles can also occur. So, CUT&Tag reactions targeting a histone modification predominantly results in fragments that are nucleosomal lengths (~180 bp), or multiples of that length. CUT&Tag targeting transcription factors predominantly produce nucleosome-sized fragments and variable amounts of shorter fragments, from neighboring nucleosomes and the factor-bound site, respectively. Tagmentation of DNA on the surface of nucleosomes also occurs, and plotting fragment lengths with single-basepair resolution reveal a 10-bp sawtooth periodicity, which is typical of successful CUT&Tag experiments.   

The smaller fragments (50-100 bp) can be due to that tethered Tn5 can tagment on the surface of a nucleosome as well as in linker regions, so the small fragments might not be background.  

    vim cut34_fragmentLen.sh
    
    #!/bin/bash

    cat filenames | while read i; 
    do
    ## Extract the 9th column from the alignment sam file which is the fragment length
    nohup samtools view \
    -@ 4 -F 0x04 ./bam/${i}_mm10_bowtie2.sam | awk \
    -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk \
    -v OFS="\t" '{print $2, $1/2}' > ./fragmentLen/${i}_fragmentLen.txt &
    done

R  

## 4.1 Filtering mapped reads by the mapping quality filtering  
Nothing

## 4.2 File format conversion  
    
    vim cut42_bam2bed.sh
    
    #!/bin/bash
    ## File format conversion ##

    cat filenames | while read i; 
    do
    ## Filter and keep the mapped read pairs
    nohup samtools view -@ 4 -bS -F 0x04 ./bam/${i}_mm10_bowtie2.sam > ./bam/${i}_mm10_bowtie2.mapped.bam &&

    ## Convert into bed file format
    bedtools bamtobed -i ./bam/${i}_mm10_bowtie2.mapped.bam -bedpe > ./bed/${i}_mm10_bowtie2.bed &&

    ## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    awk '$1==$4 && $6-$2 < 1000 {print $0}' ./bed/${i}_mm10_bowtie2.bed > ./bed/${i}_mm10_bowtie2.clean.bed &&

    ## Only extract the fragment related columns
    cut -f 1,2,6 ./bed/${i}_mm10_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n > ./bed/${i}_mm10_bowtie2.fragments.bed &
    done

## 4.3 Assess replicate reproducibility  
    
    vim cut43_bin500.sh
    
    #!/bin/bash
    ## Assess replicate reproducibility ##

    cat filenames | while read i; 
    do
    ## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.
    awk -v w=500 '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ./bed/${i}_mm10_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n > ./bed/${i}_mm10_bowtie2.fragmentsCount.bin500.bed &
    done

R  

## 5.1 Scaling factor   

    samtools view -@ 4 -F 0x04 ./bam/CFA3-1_ecoil_bowtie2.sam | wc -l > ./bowtie2_summary/CFA3-1_ecoil_bowtie2.seqDepthDouble &
    samtools view -@ 4 -F 0x04 ./bam/CFA3-2_ecoil_bowtie2.sam | wc -l > ./bowtie2_summary/CFA3-2_ecoil_bowtie2.seqDepthDouble &
    samtools view -@ 4 -F 0x04 ./bam/CFA3-B1_ecoil_bowtie2.sam | wc -l > ./bowtie2_summary/CFA3-B1_ecoil_bowtie2.seqDepthDouble &
    samtools view -@ 4 -F 0x04 ./bam/NADCFA3-1_ecoil_bowtie2.sam | wc -l > ./bowtie2_summary/NADCFA3-1_ecoil_bowtie2.seqDepthDouble &
    samtools view -@ 4 -F 0x04 ./bam/NADCFA3-2_ecoil_bowtie2.sam | wc -l > ./bowtie2_summary/NADCFA3-2_ecoil_bowtie2.seqDepthDouble &

    seqDepth=$((seqDepthDouble/2))
    scale_factor=`echo "10000 / $seqDepth" | bc -l`
    echo $scale_factor

    bedtools genomecov \
    -bg -scale $scale_factor \
    -i ./bed/${i}_mm10_bowtie2.fragments.bed \
    -g /home/yangjiajun/downloads/genome/mm10_GRCm38/ucsc_fa/mm10.chrom.sizes > ./bedgraph/${i}_mm10_bowtie2.fragments.normalized.bedgraph

    vim cut5_bedgraph.sh

    #!/bin/bash
    ## Spike-in calibration ##

    cat filenames | while read i; 
    do
    bedtools genomecov -bg -i ./bed/${i}_mm10_bowtie2.fragments.bed \
    -g /home/yangjiajun/downloads/genome/mm10_GRCm38/ucsc_fa/mm10.chrom.sizes > ./bedgraph/${i}_mm10_bowtie2.fragments.normalized.bedgraph &
    done

## 6.1. SEACR  

The Sparse Enrichment Analysis for CUT&RUN, SEACR, package is designed to call peaks and enriched regions from chromatin profiling data with very low backgrounds (i.e., regions with no read coverage) that are typical for CUT&Tag chromatin profiling experiments. SEACR requires bedGraph files from paired-end sequencing as input and defines peaks as contiguous blocks of basepair coverage that do not overlap with blocks of background signal delineated in the IgG control dataset. SEACR is effective for calling both narrow peaks from factor binding sites and broad domains characteristic of some histone modifications. The description of the method is published at Meers et al. 2019, and the user’s manual is available on github. Since we have normalized fragment counts with the E. coli read count, we set the normalization option of SEACR to “non”. Otherwise, the “norm” is recommended.  

Here, we used "norm".  

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/CFA3-1_mm10_bowtie2.fragments.normalized.bedgraph \
         ./bedgraph/CFA3-B1_mm10_bowtie2.fragments.normalized.bedgraph \
         norm stringent ./SEACR/CFA3-1_seacr_control.peaks &

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/CFA3-2_mm10_bowtie2.fragments.normalized.bedgraph \
         ./bedgraph/CFA3-B1_mm10_bowtie2.fragments.normalized.bedgraph \
         norm stringent ./SEACR/CFA3-2_seacr_control.peaks &

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/NADCFA3-1_mm10_bowtie2.fragments.normalized.bedgraph \
         ./bedgraph/CFA3-B1_mm10_bowtie2.fragments.normalized.bedgraph \
         norm stringent ./SEACR/NADCFA3-1_seacr_control.peaks &

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/NADCFA3-2_mm10_bowtie2.fragments.normalized.bedgraph \
         ./bedgraph/CFA3-B1_mm10_bowtie2.fragments.normalized.bedgraph \
         norm stringent ./SEACR/NADCFA3-2_seacr_control.peaks &

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/CT-H_mm10_bowtie2.fragments.normalized.bedgraph \
         ./bedgraph/bg-B2_mm10_bowtie2.fragments.normalized.bedgraph \
         norm stringent ./SEACR/CT-H_seacr_control.peaks &

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/CT-S1_mm10_bowtie2.fragments.normalized.bedgraph \
         ./bedgraph/bg-B2_mm10_bowtie2.fragments.normalized.bedgraph \
         norm stringent ./SEACR/CT-S1_seacr_control.peaks &

