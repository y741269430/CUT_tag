# CUT&Tag 
参考：https://yezhengstat.github.io/CUTTag_tutorial/index.html    

实验大概流程：      
- 1.利用连接有刀豆蛋白 A 的磁珠结合细胞（刀豆蛋白 A 能与细胞膜或细胞核上的糖蛋白结合）；      
- 2.使用非离子去污剂洋地黄皂苷(digitonin)使细胞膜通透，加入目标蛋白的一抗抗体（假设是鼠源）进行孵育（抗体能够通过细胞膜、核孔进入细胞核与目标蛋白结合），
- 3.加入二抗（抗鼠）增强靶向结合能力；     
- 4.加入pG/pA（Tn5）孵育，其中pG上源自于链球菌G族的细胞表面蛋白，而pA源自于A型金黄色葡萄球菌的细胞壁蛋白。pG/pA（Tn5）带有Tn5，同时也带有接头。pG/pA与二抗结合，把Tn5带到目标蛋白附近，切割目标蛋白附近的序列，使目标蛋白和其结合的 DNA 序列从染色质上脱离下来并游离到细胞外；       
- 5.提取DNA，建库测序。     

生信分析流程：   
- 1.当我们数据下机之后，得到的fastq文件。使用`FastQC`软件对raw data进行质量评估。后续clean data同样需要评估。   
- 2.使用`Trimmomatic`软件对原始数据进行质控（这一步主要是去除3’端的接头污染、去除低质量序列（保留MPAQ >= 30））。得到clean data。
- 3.将clean data使用`bowtie2`软件与基因组进行比对，得到的sam文件使用`samtools`转换成bam。
- 4.得到的bam文件，获取其唯一比对以及去重复reads的结果bam文件。
- 5.使用`Deeptools`绘制TSS, Peak center 或GeneBody富集热图（依组学而定），展示数据在这些区域及前后3kb上的富集情况。
- 6.使用`SEACR`,`MACS2`或`MACS3`进行peak calling。
- 7.使用`IDR`软件进行样品间高可信度的峰筛选.
- 8.将bam文件转换成bigwig文件，使用`IGV`进行可视化。
- 9.使用r包`ChIPseeker`对peak进行注释。
- 10.使用`homer`或`MEME`进行motif预测。
- 11.使用`MAnorm`（无生物学重复）或`DiffBind`（有生物学重复）进行差异peak分析.


## 0. 创建环境         
    
    conda create -n cuttag  
    conda activate cuttag  

    conda install -c bioconda picard  
    conda install -c bioconda samtools
    conda install -c bioconda bedtools
    conda install -c bioconda deeptools
    conda install -c bioconda seacr
    conda install -c bioconda bowtie2
    conda install -c bioconda fastqc

## 0. 利用bowtie2构建小鼠基因组（mm39）索引以及ecoil基因组索引（构建一次以后都不用做了）        

    cd /home/yangjiajun/downloads/genome/mm39_GRCm39/ucsc_fa/

    samtools faidx GRCm38.primary_assembly.genome.fa &
    cut -f1,2 GRCm38.primary_assembly.genome.fa.fai > mm39.chrom.sizes &
    
    nohup bowtie2-build GRCm38.primary_assembly.genome.fa \
    /home/yangjiajun/downloads/genome/mm39_GRCm39/bowtie2_idx/mm10 &  


    cd /home/yangjiajun/downloads/genome/ecoil_U00096.3/ucsc_fa/

    samtools faidx GCF_000005845.2_ASM584v2_genomic.fa &
    cut -f1,2 GCF_000005845.2_ASM584v2_genomic.fa.fai > ecoil.chrom.sizes &
    
    nohup bowtie2-build GCF_000005845.2_ASM584v2_genomic.fa \
    /home/yangjiajun/downloads/genome/ecoil_U00096.3/bowtie2_idx/ecoil &   

----

## 1. 激活环境，创建所需文件夹           
    
    conda activate cuttag  
    
    mkdir raw clean bam mapbam bed bowtie2_summary picard_summary bedgraph SEACR fragmentLen plot trim  
    
## 2. 写入filenames     
    
    ls raw/*1.fq.gz |cut -d "_" -f 1 |cut -d "/" -f 2 |cut -d "-" -f 1-5 > filenames
    
## (Quick Start)  

    pre_trim.sh
    cut1_bw2.sh  
    cut42_bam2bed.sh  
    cut5_bedgraph.sh  
    cut6_seacr005.sh
    cut7_sort_idx.sh
    cut8_bw.sh  

## ~~Remove adapter~~  
    vim pre_trim.sh  
    
    #!/bin/bash
    ## trim_galore ##
    
    cat filenames | while read i; 
    do
    # paired end
    nohup trim_galore -q 25 --phred33 --length 20 -e 0.1 --stringency 1 --paired ./${i}*_1.fq.gz ./${i}*_2.fq.gz -o ./trim &
    
    done

## 利用trimmomatic去除接头(Illumina)   
```bash
#!/bin/bash
## trimmomatic ##

cat filenames | while read i; 
do
nohup trimmomatic PE -phred33 -threads 4 \
./RawData/${i}/${i}*_R1_001.fastq.gz \
./RawData/${i}/${i}*_R2_001.fastq.gz \
./trim/${i}_forward_paired.fq.gz \
./trim/${i}_forward_unpaired.fq.gz \
./trim/${i}_reverse_paired.fq.gz \
./trim/${i}_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &

done
```

## 3.1.1 比对到mm39  
注(报错)：Could not locate a Bowtie index corresponding to basename (下方的 ${mm39} 需要加上绝对路径(/home)而不是相对路径(~/))
    
    vim cut1_bw2.sh
    
    #!/bin/bash
    ## Alignment to mm39 ##

    mm39="/home/jjyang/downloads/genome/mm39_GRCm39/bowtie2_idx/mm39"

    cat filenames | while read i; 
    do
    nohup bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 \
    -I 10 -X 700 -p 4 -x ${mm39} \
    -1 trim/${i}_forward_paired.fq.gz \
    -2 trim/${i}_reverse_paired.fq.gz \
    -S ./bam/${i}_mm39_bowtie2.sam &> ./bowtie2_summary/${i}_mm39_bowtie2.txt &
    done

## 3.1.2 比对到 spike-in 基因组 (ecoil)   

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

## 3.2 比对并总结报告  

Summarize the raw reads and uniquely mapping reads to report the efficiency of alignment. Alignment frequencies are expected to be >80% for high-quality data. CUT&Tag data typically has very low backgrounds, so as few as 1 million mapped fragments can give robust profiles for a histone modification in the human genome. Profiling of less-abundant transcription factors and chromatin proteins may require 10 times as many mapped fragments for downstream analysis.  

总结原始reads和唯一比对reads，以报告比对的效率。对于高质量数据，比对频率预计为>80%。一般来说，CUT&Tag数据的背景非常低，因此，只需100万个比对上的片段就可以为人类基因组中的组蛋白修饰提供可靠的profiles。低丰度的转录因子和染色质蛋白的谱分析可能需要10倍于下游分析的图谱片段。  

## 3.2.1 测序深度        
R  

## 3.2.2 Spike-in 比对率    
R  

## 3.2.3 统计 mm39 和 E.coli 的比对率  
R  

## 3.2.4 可视化测序深度以及比对结果         
R  

In a typical CUT&Tag experiment targeting the abundant H3K27me3 histone modification in 65,000 K562 cells, the percentage of E. coli reads range from ~0.01% to 10%. With fewer cells or less abundant epitopes, E. coli reads can comprise as much as 70% or the total mapped reads. For IgG controls, the percentage of E. coli reads is typically much higher than that for an abundant histone modification.   

在一项针对65,000个K562细胞中富集H3K27me3组蛋白修饰的CUT&Tag实验中，E.coli的reads比例在~0.01%至10%之间。如果细胞数量少或表位数量少，E.coli的reads可占总reads的70%。对于IgG对照，E.coli的reads比例通常比组蛋白修饰的要高得多。  

## 3.3 是否去除重复序列        

CUT&Tag integrates adapters into DNA in the vicinity of the antibody-tethered pA-Tn5, and the exact sites of integration are affected by the accessibility of surrounding DNA. For this reason fragments that share exact starting and ending positions are expected to be common, and such ‘duplicates’ may not be due to duplication during PCR. In practice, we have found that the apparent duplication rate is low for high quality CUT&Tag datasets, and even the apparent ‘duplicate’ fragments are likely to be true fragments. Thus, we do not recommend removing the duplicates. In experiments with very small amounts of material or where PCR duplication is suspected, duplicates can be removed.     

因为CUT&Tag的特性，将会整合接头到pA-Tn5结合抗体周围的DNA区域中，且整合的位点受到DNA可及性的影响。所以在此因素的影响下，获取的片段有着一致的起始与终止位点是可能的，而不是因为PCR扩增的原因。在高质量的文库中，不需要进行CUT&Tag重复序列的去除，其重复率是很低的，“重复”片段也很可能是真实的片段。而当文库初始浓度、总量很少或者PCR扩增次数很多的时候，重复序列需要进行去除。所以一般推荐把重复序列去除。特别是一些检测转录因子结合信号的，一般重复序列的比例都很高，不大可能是由CUT&Tag的特性导致的，更可能的原因应该是PCR重复。    

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

## 3.4. 评估比对片段大小的分布         

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

## 4.1 比对结果的过滤与格式转换       
Nothing

## 4.2 文件格式转换       
    
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

## 4.3 评估重复样本的重复性       
    
    vim cut43_bin500.sh (optional)
    
    #!/bin/bash
    ## Assess replicate reproducibility ##

    cat filenames | while read i; 
    do
    ## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.
    awk -v w=500 '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ./bed/${i}_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n > ./bed/${i}_bowtie2.fragmentsCount.bin500.bed &
    done

R  

## 5.1 Spike-in校正      
能进行校正的基本假设是：在一系列使用相同数量细胞的样本中，比对到E.coli基因组上reads的比例是相同的。由于这个假设，分析流程没有进行样品间和批次间的标准化，这会使得残余的E.coli DNA数量差异很大。

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
    
    
