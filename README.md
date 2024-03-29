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

    cd /home/yangjiajun/downloads/genome/mm39_GRCm39/ucsc_fa/

    samtools faidx GRCm38.primary_assembly.genome.fa &
    cut -f1,2 GRCm38.primary_assembly.genome.fa.fai > mm10.chrom.sizes &
    
    nohup bowtie2-build GRCm38.primary_assembly.genome.fa \
    /home/yangjiajun/downloads/genome/mm39_GRCm39/bowtie2_idx/mm10 &  


    cd /home/yangjiajun/downloads/genome/ecoil_U00096.3/ucsc_fa/

    samtools faidx GCF_000005845.2_ASM584v2_genomic.fa &
    cut -f1,2 GCF_000005845.2_ASM584v2_genomic.fa.fai > ecoil.chrom.sizes &
    
    nohup bowtie2-build GCF_000005845.2_ASM584v2_genomic.fa \
    /home/yangjiajun/downloads/genome/ecoil_U00096.3/bowtie2_idx/ecoil &   

----

## 1. Activate the source and create the folder  
    
    conda activate cuttag  
    
    mkdir raw clean bam mapbam bed bowtie2_summary picard_summary bedgraph SEACR fragmentLen plot trim  
    
## 2. Write the filenames  
    
    ls raw/*1.fq.gz |cut -d "_" -f 1 |cut -d "/" -f 2 |cut -d "-" -f 1-5 > filenames
    
## (Quick Start)  

    pre_trim.sh
    cut1_bw2.sh  
    cut42_bam2bed.sh  
    cut5_bedgraph.sh  
    cut6_seacr005.sh
    cut7_sort_idx.sh
    cut8_bw.sh  

## Remove adapter  
    vim pre_trim.sh  
    
    #!/bin/bash
    ## trim_galore ##
    
    cat filenames | while read i; 
    do
    # paired end
    nohup trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 1 --paired ./${i}*_1.fq.gz ./${i}*_2.fq.gz -o ./trim &
    
    done

## 3.1.1 Alignment to mm39  
注(报错)：Could not locate a Bowtie index corresponding to basename (下方的 ${mm39} 需要加上绝对路径(/home)而不是相对路径(~/))
    
    vim cut1_bw2.sh
    
    #!/bin/bash
    ## Alignment to mm39 ##

    mm39="/home/yangjiajun/downloads/genome/mm39_GRCm39/bowtie2_idx/mm39"

    cat filenames | while read i; 
    do
    nohup bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 \
    -I 10 -X 700 -p 4 -x ${mm39} \
    -1 ./trim/${i}_1_val_1.fq.gz \
    -2 ./trim/${i}_2_val_2.fq.gz \
    -S ./bam/${i}_mm39_bowtie2.sam &> ./bowtie2_summary/${i}_mm39_bowtie2.txt &
    done

## 3.1.2 Alignment to spike-in genome for spike-in calibration (ecoil)   

E. coli DNA is carried along with bacterially-produced pA-Tn5 protein and gets tagmented non-specifically during the reaction. The fraction of total reads that map to the E.coli genome depends on the yield of epitope-targeted CUT&Tag, and so depends on the number of cells used and the abundance of that epitope in chromatin. Since a constant amount of pATn5 is added to CUT&Tag reactions and brings along a fixed amount of E. coli DNA, E. coli reads can be used to normalize epitope abundance in a set of experiments.  

大肠杆菌的DNA与细菌产生的pA-Tn5蛋白一起携带，并在反应过程中受到非特异性标记。定位到大肠杆菌基因组的总reads的比例取决于表位靶向的CUT&Tag的产量，因此也取决于所使用的细胞数量和染色质中表位的丰度。由于在CUT&Tag反应中加入一定量的pA-Tn5并带来一定量的大肠杆菌DNA，因此大肠杆菌reads可用于一系列实验中表位丰度的标准化。  

    vim cut1_ecoil_bw2.sh
    
    #!/bin/bash
    ## Alignment to spike-in genome for spike-in calibration ##

    spikeInRef="/home/yangjiajun/downloads/genome/ecoil_U00096.3/bowtie2_idx/ecoil"

    cat filenames | while read i; 
    do
    nohup bowtie2 --end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 \
    -I 10 -X 700 -p 8 -x ${spikeInRef} \
    -1 ./raw/${i}_1.fq.gz \
    -2 ./raw/${i}_2.fq.gz \
    -S ./bam/${i}_ecoil_bowtie2.sam &> ./bowtie2_summary/${i}_ecoil_bowtie2.txt &
    done

## 3.2 Report sequencing mapping summary  

Summarize the raw reads and uniquely mapping reads to report the efficiency of alignment. Alignment frequencies are expected to be >80% for high-quality data. CUT&Tag data typically has very low backgrounds, so as few as 1 million mapped fragments can give robust profiles for a histone modification in the human genome. Profiling of less-abundant transcription factors and chromatin proteins may require 10 times as many mapped fragments for downstream analysis.  

总结原始reads和唯一比对reads，以报告比对的效率。对于高质量数据，比对频率预计为>80%。一般来说，CUT&Tag数据的背景非常低，因此，只需100万个比对上的片段就可以为人类基因组中的组蛋白修饰提供可靠的profiles。低丰度的转录因子和染色质蛋白的谱分析可能需要10倍于下游分析的图谱片段。  

## 3.2.1 Sequencing depth  
R  

## 3.2.2 Spike-in alignment  
R  

## 3.2.3 Summarize the alignment to mm10 and E.coli  
R  

## 3.2.4 Visualizing the sequencing depth and alignment results.  

R  

In a typical CUT&Tag experiment targeting the abundant H3K27me3 histone modification in 65,000 K562 cells, the percentage of E. coli reads range from ~0.01% to 10%. With fewer cells or less abundant epitopes, E. coli reads can comprise as much as 70% or the total mapped reads. For IgG controls, the percentage of E. coli reads is typically much higher than that for an abundant histone modification.   

在一项针对65,000个K562细胞中富集H3K27me3组蛋白修饰的CUT&Tag实验中，E.coli的reads比例在~0.01%至10%之间。如果细胞数量少或表位数量少，E.coli的reads可占总reads的70%。对于IgG对照，E.coli的reads比例通常比组蛋白修饰的要高得多。  

## 3.3 Remove duplicates   

CUT&Tag integrates adapters into DNA in the vicinity of the antibody-tethered pA-Tn5, and the exact sites of integration are affected by the accessibility of surrounding DNA. For this reason fragments that share exact starting and ending positions are expected to be common, and such ‘duplicates’ may not be due to duplication during PCR. In practice, we have found that the apparent duplication rate is low for high quality CUT&Tag datasets, and even the apparent ‘duplicate’ fragments are likely to be true fragments. Thus, we do not recommend removing the duplicates. In experiments with very small amounts of material or where PCR duplication is suspected, duplicates can be removed. The following commands show how to check the duplication rate using Picard.    

CUT&Tag将adaptors整合到抗体栓系pA-Tn5附近的DNA中，并且整合的准确位置受周围DNA可及性的影响。由于这个原因，共享精确起始和结束位置的片段是常见的，而这种“duplicates”可能不是由于PCR过程中的重复。在实践中，我们发现对于高质量的CUT&Tag数据集，重复率是很低的，甚至“重复”片段也很可能是真实的片段。因此，我们不建议删除duplicates。在用非常少量的材料或怀疑是PCR重复的时候，duplicates可以被去除。下面的命令展示了如何使用Picard检查重复率。  

    vim cut2_picard.sh
    
    #!/bin/bash
    ## Remove duplicates ##

    cat filenames | while read i; 
    do
    ## Sort by coordinate
    nohup picard SortSam -I ./bam/${i}_mm39_bowtie2.sam \
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

（1）In these example datasets, the IgG control samples have relatively high duplication rates, since reads in this sample derive from non-specific tagmentation in the CUT&Tag reactions. Therefore, it is appropriate to remove the duplicates from the IgG datasets before downstream analysis.  
（2）The estimated library size are the estimated number of unique molecules in the library based on PE duplication calculated by Picard.  
（3）The estimated library sizes is proportional to the abundance of the targeted epitope and to the quality of the antibody used, while the estimated library sizes of IgG samples are expected to be very low.  
（4）Unique fragment number is calculated by the MappedFragNum_hg38 * (1-DuplicationRate/100).  

（1）在这些样本数据集中，IgG对照样本有相对较高的重复率，因为该样本中的reads来自于CUT&Tag反应中的非特异性标记。因此，在进行下游分析之前，应该将重复的IgG数据集删除。  
（2）文库的大小是Picard根据PE重复计算出的文库中唯一分子的数量。  
（3）估计的文库大小与靶标表位的丰度和所用抗体的质量成正比，而IgG样本的文库估计大小非常低。  
（4）唯一的片段数由MappedFragNum_hg38 * (1-DuplicationRate/100)计算。  

## 3.4. Assess mapped fragment size distribution  

CUT&Tag inserts adapters on either side of chromatin particles in the vicinity of the tethered enzyme, although tagmentation within chromatin particles can also occur. So, CUT&Tag reactions targeting a histone modification predominantly results in fragments that are nucleosomal lengths (~180 bp), or multiples of that length. CUT&Tag targeting transcription factors predominantly produce nucleosome-sized fragments and variable amounts of shorter fragments, from neighboring nucleosomes and the factor-bound site, respectively. Tagmentation of DNA on the surface of nucleosomes also occurs, and plotting fragment lengths with single-basepair resolution reveal a 10-bp sawtooth periodicity, which is typical of successful CUT&Tag experiments.   

CUT&Tag在被酶栓附近的染色质颗粒的两侧插入adaptors，虽然染色质颗粒内的标记也可能发生。因此，针对组蛋白修饰的CUT&Tag反应主要导致核小体长度(~180 bp)或数倍于该长度的片段。CUT&Tag靶向转录因子主要产生核小体大小的片段和不同数量的短片段，分别来自邻近的核小体和因子结合位点。DNA的标记在核小体表面的也发生，用单碱基对分辨率绘制片段长度显示了10 bp的锯齿状周期性，这是成功的CUT&Tag实验的典型特征。  

The smaller fragments (50-100 bp) can be due to that tethered Tn5 can tagment on the surface of a nucleosome as well as in linker regions, so the small fragments might not be background.  

较小的片段(50-100 bp)可能是由于Tn5被栓系在核小体表面以及linker区域，因此小片段可能不是背景。  

    vim cut34_fragmentLen.sh (optional)
    
    #!/bin/bash

    cat filenames | while read i; 
    do
    ## Extract the 9th column from the alignment sam file which is the fragment length
    nohup samtools view \
    -@ 4 -F 0x04 ./bam/${i}_mm39_bowtie2.sam | awk \
    -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk \
    -v OFS="\t" '{print $2, $1/2}' > ./fragmentLen/${i}_fragmentLen.txt &
    done

R  

## 4.1 Filtering mapped reads by the mapping quality filtering  
Nothing

## 4.2 File format conversion  
    
No remove duplication  

    vim cut42_bam2bed.sh
    
    #!/bin/bash
    ## File format conversion ##

    cat filenames | while read i; 
    do
    ## Filter and keep the mapped read pairs
    nohup samtools view -@ 8 -bS -F 0x04 ./bam/${i}_mm39_bowtie2.sam > ./mapbam/${i}_bowtie2.mapped.bam &&

    ## Convert into bed file format
    bedtools bamtobed -i ./mapbam/${i}_bowtie2.mapped.bam -bedpe > ./bed/${i}_bowtie2.bed &&

    ## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    awk '$1==$4 && $6-$2 < 1000 {print $0}' ./bed/${i}_bowtie2.bed > ./bed/${i}_bowtie2.clean.bed &&

    ## Only extract the fragment related columns
    cut -f 1,2,6 ./bed/${i}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n > ./bed/${i}_bowtie2.fragments.bed &
    done
    
Remove duplication  

    vim cut42_bam2bed_rmd.sh
    
    #!/bin/bash
    ## File format conversion ##

    cat filenames | while read i; 
    do
    ## Filter and keep the mapped read pairs
    nohup samtools view -@ 8 -bS -F 0x04 ./bam/${i}_bowtie2.sorted.rmDup.sam > ./mapbam/${i}_bowtie2.mapped.rmDup.bam &&

    ## Convert into bed file format
    bedtools bamtobed -i ./mapbam/${i}_bowtie2.mapped.rmDup.bam -bedpe > ./bed/${i}_bowtie2.rmDup.bed &&

    ## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
    awk '$1==$4 && $6-$2 < 1000 {print $0}' ./bed/${i}_bowtie2.rmDup.bed > ./bed/${i}_bowtie2.clean.rmDup.bed &&

    ## Only extract the fragment related columns
    cut -f 1,2,6 ./bed/${i}_bowtie2.clean.rmDup.bed | sort -k1,1 -k2,2n -k3,3n > ./bed/${i}_bowtie2.fragments.rmDup.bed &
    done

## 4.3 Assess replicate reproducibility  
    
    vim cut43_bin500.sh (optional)
    
    #!/bin/bash
    ## Assess replicate reproducibility ##

    cat filenames | while read i; 
    do
    ## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.
    awk -v w=500 '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ./bed/${i}_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n > ./bed/${i}_bowtie2.fragmentsCount.bin500.bed &
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
    
No remove duplication  

    vim cut5_bedgraph.sh

    #!/bin/bash
    ## Spike-in calibration ##

    cat filenames | while read i; 
    do
    bedtools genomecov -bg -i ./bed/${i}_bowtie2.fragments.bed \
    -g /home/yangjiajun/downloads/genome/mm39_GRCm38/ucsc_fa/mm39.chrom.sizes > ./bedgraph/${i}_bowtie2.fragments.normalized.bedgraph &
    done
    
Remove duplication  

    vim cut5_bedgraph_rmd.sh
    
    #!/bin/bash
    ## Spike-in calibration ##

    cat filenames | while read i; 
    do
    bedtools genomecov -bg -i ./bed/${i}_bowtie2.fragments.rmDup.bed \
    -g /home/yangjiajun/downloads/genome/mm39_GRCm38/ucsc_fa/mm39.chrom.sizes > ./bedgraph/${i}_bowtie2.fragments.normalized.rmDup.bedgraph &
    done

## 6.1. SEACR  

The Sparse Enrichment Analysis for CUT&RUN, SEACR, package is designed to call peaks and enriched regions from chromatin profiling data with very low backgrounds (i.e., regions with no read coverage) that are typical for CUT&Tag chromatin profiling experiments. SEACR requires bedGraph files from paired-end sequencing as input and defines peaks as contiguous blocks of basepair coverage that do not overlap with blocks of background signal delineated in the IgG control dataset. SEACR is effective for calling both narrow peaks from factor binding sites and broad domains characteristic of some histone modifications. The description of the method is published at Meers et al. 2019, and the user’s manual is available on github. Since we have normalized fragment counts with the E. coli read count, we set the normalization option of SEACR to “non”. Otherwise, the “norm” is recommended.  

对CUT&RUN的稀疏富集分析SEACR包的设计用于从非常低的背景(即没有read覆盖的区域)的染色质中富集区域并call peaks，这是典型的CUT&Tag染色质分析实验。SEACR需要来自pairedend测序的bedGraph文件作为输入，并将峰定义为不与IgG对照数据集中描绘的背景信号块重叠的连续碱基对覆盖块。SEACR既能有效地从因子结合位点calling窄峰，也能calling出某些组蛋白修饰特征的宽域。由于我们已经使用大肠杆菌read count对片段计数进行了规范化，所以我们将SEACR的规范化选项设置为“non”。否则，建议使用“norm”。  

Here, we used "0.05".  

    vim cut6_seacr005.sh
    
    #!/bin/bash
    ## seacr 0.05 ##
    
    cat filenames | while read i; 
    do
    bash /home/jjyang/.conda/envs/cuttag/bin/SEACR_1.3.sh \
    ./bedgraph/${i}_bowtie2.fragments.normalized.bedgraph \
    0.05 norm stringent ./seacr005/${i}_005.peaks &
    done

Here, we used "norm".  

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/CFA3-1_mm10_bowtie2.fragments.normalized.bedgraph \
    ./bedgraph/BCFA3-1_mm10_bowtie2.fragments.normalized.bedgraph \
    norm stringent ./SEACR/CFA3-1_seacr_BCFA3-1.peaks &

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/CFA3-2_mm10_bowtie2.fragments.normalized.bedgraph \
    ./bedgraph/BCFA3-1_mm10_bowtie2.fragments.normalized.bedgraph \
    norm stringent ./SEACR/CFA3-2_seacr_BCFA3-1.peaks &

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/CFA3-3_mm10_bowtie2.fragments.normalized.bedgraph \
    ./bedgraph/BCFA3-1_mm10_bowtie2.fragments.normalized.bedgraph \
    norm stringent ./SEACR/CFA3-3_seacr_BCFA3-1.peaks &

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/CFA3-4_mm10_bowtie2.fragments.normalized.bedgraph \
    ./bedgraph/BCFA3-1_mm10_bowtie2.fragments.normalized.bedgraph \
    norm stringent ./SEACR/CFA3-4_seacr_BCFA3-1.peaks &

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/NADCFA3-1_mm10_bowtie2.fragments.normalized.bedgraph \
    ./bedgraph/BCFA3-1_mm10_bowtie2.fragments.normalized.bedgraph \
    norm stringent ./SEACR/NADCFA3-1_seacr_BCFA3-1.peaks &

    bash ~/miniconda3/envs/cuttag/bin/SEACR_1.3.sh ./bedgraph/NADCFA3-2_mm10_bowtie2.fragments.normalized.bedgraph \
    ./bedgraph/BCFA3-1_mm10_bowtie2.fragments.normalized.bedgraph \
    norm stringent ./SEACR/NADCFA3-2_seacr_BCFA3-1.peaks &

## 7. sort index

Sort bam files and index for diffbind in R.  

    vim cut7_sort_idx.sh

    #!/bin/bash
    ## sort bam and index for diffbind or igv ##

    cat filenames | while read i; 
    do
    nohup samtools sort -@ 8 ./mapbam/${i}_mm10_bowtie2.mapped.bam -o ./sortbam/${i}_mm10_bowtie2.sorted.bam &&
    samtools index -@ 8 ./sortbam/${i}_mm10_bowtie2.sorted.bam &
    done

## 8. bw file in IGV  

Transfer to BW in IGV.  

    vim cut8_bw.sh  
    
    #!/bin/bash
    ## make bigwig (deeptools) ##

    cat filenames | while read i; 
    do
    nohup bamCoverage --bam ./sortbam/${i}_mm10_bowtie2.sorted.bam -o bw_file/${i}.bw \
        --binSize 10 \
        --normalizeUsing RPKM \
        --numberOfProcessors 10 \
        --effectiveGenomeSize 2652783500 &
    done
    
    
