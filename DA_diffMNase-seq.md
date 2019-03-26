# Processing of cotton MNase-seq datasets
---

**Outline** 

1. [Dns-MNase-seq data preprocessing](#)
2. [Read mapping and quality control](##2.-Read-mapping-and-quality-control)
3. [Differential nuclease sensitivity profiling analysis]()
4. [Nucleosome positioning analysis](#4.-nucleosome-positioning-analysis)
5.[test](#5-here-is-test)

## 1. Dns-MNase-seq data preprocessing

### MNase-seq datasets
Short-read data will be deposited in the NCBI short read archive ([SRP??????](http://trace.ddbj.nig.ac.jp/DRASearch/study?acc=SRP??????)), also as Biobroject [PRJNA??????](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA??????).

Paired-end reads of 16 samples: 4 species (A2, D5, A2xD5, Maxxa) X 2 digestive conditions (Heavy and Light) X 2 technical reps.

    cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/
    cd rawfastq
    ln -s /lss/research/jfw-lab/RawData/HGJ_leafMNase-seq/WTNHHW163125/Feb2017/*gz .
    ln -s /lss/research/jfw-lab/RawData/HGJ_leafMNase-seq/WTNHHW163125/Feb2018/*gz .

#### Checking read quality with [FastQC](http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/)
    cd ..
    module load fastqc/0.11.3
    module load py-multiqc
    mkdir QCreport
    mkdir QCreport/raw
    mkdir QCreport/trimmed
    fastqc -o QCreport/raw rawfastq/*gz
    # results were collected from all samples into a single report for easy comparison with MultiQC (Ewels et al., 2016)
    cd QCreport/raw
    multiqc .
    
#### Quality trimming and adaptor removal
We usually use [Sickle](https://github.com/najoshi/sickle) to trim off sequences below quality threshold, and another popular tool is  [Fastx toolkit](http://hannonlab.cshl.edu/fastx_toolkit/). One more alternative to remove adpaters or primers is [cutadapt](https://cutadapt.readthedocs.io/). Based FastQC results above, illumina universal adaptor removal is necessay, and [Trim Galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) appears to be a really nice wrapper tool of FastQC and Cutadapt, quite easy to use! 

    module load python 
    # python needed for Cutadapt
    trim_galore_v0.4.2/trim_galore --paired -o trimmed/ rawfastq/M1H_1.fq.gz rawfastq/M1H_2.fq.gz
    # check results
    grep 'Total reads processed' trimmed/*report.txt >trimmed/summary.txt
    grep 'Reads with adapters' trimmed/*report.txt >>trimmed/summary.txt
    grep 'Total written' trimmed/*report.txt >>trimmed/summary.txt
    grep 'Number of sequence pairs' trimmed/*report.txt >>trimmed/summary.txt
    # QC again
    fastqc -o QCreport/trimmed/ trimmed/*val*
    cd QCreport/trimmed
    multiqc .

### Cotton reference genomes
[CottonGen](https://www.cottongen.org/data/download/genome#Ass) compiles all published cotton genomes, and I will need 4 different reference genomes for AD1, A2, D5 and A2xD5. An improved AD1 reference became available on [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Ghirsutum_er).

    mkdir refGenomes
    cd refGenomes
    module load bowtie2
    
    ### D5_JGI - Paterson et al. 2012 Nature
    ln -s ~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/Dgenome2_13.fasta
    
    ### AD1_458 - Saski et al, 2017 
    ln -s ~/jfw-lab/GenomicResources/archived_resources/AD1Saski/v1.1/assembly/Ghirsutum_458_v1.0.fa.gz
    zcat Ghirsutum_458_v1.0.fa.gz |grep -n '>scaffold'|head -10 
    # 27134894:>scaffold_27
    zcat Ghirsutum_458_v1.0.fa.gz |head -27134893 >TM1new_26.fasta 
  
    ### A2Du2018 - Du et al, 2018 
    ln -s ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Du2018/1.PacBio-Gar-Assembly-v1.0/G.arboreum.Chr.v1.0.fasta.gz
    zcat G.arboreum.Chr.v1.0.fasta.gz |grep -n '>tig'|head -10 
    # 40:>tig00000043
    zcat G.arboreum.Chr.v1.0.fasta.gz |head -39 >A2Du_13.fasta
    
    ### A2Du + D5 as ref for A2xD5
    sed 's/>Chr/>A/g' A2Du_13.fasta >At.fasta
    sed 's/>Chr/>D/g' Dgenome2_13.fasta >Dt.fasta
    cat At.fasta Dt.fasta >F1_26.fasta
    grep '>' F1_26.fasta 
    rm At.fasta
    rm Dt.fasta

    #-----------OLD---------------
    ### A2_BGI - Li et al. 2014 Nature Genetics
    ln -s ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Li/A2genome_13.fasta
    ### AD1_NBI - Zhang et al, 2016 Nature biotechnology
    ln -s ~/jfw-lab/GenomicResources/archived_resources/gmapdb/AD1TM1/TM1.fasta
    grep -n '>scaffold' TM1.fasta |head -10   #32244283:>scaffold27_A01
    head -32244282 TM1.fasta >TM1_26.fasta 
    ### make my own ref for A2xD5
    cat A2genome_13.fasta Dgenome2_13.fasta >F1_26t.fasta
    grep '>' F1_26t.fasta 
    sed 's/>Chr/>D5_chr/g' F1_26t.fasta >F1_26.fasta
    grep '>' F1_26.fasta
    rm F1_26t.fasta
    
Now build bowtie reference

    bowtie2-build TM1_26.fasta TM1
    bowtie2-build TM1new_26.fasta TM1new
    bowtie2-build F1_26.fasta F1
    bowtie2-build A2genome_13.fasta A2
    bowtie2-build Dgenome2_13.fasta D5
    bowtie2-build A2Du_26.fasta A2Du

## 2.Read mapping and quality control

### 2.1 Bowtie2 mapping
Default setting `-k 1` report 1 alignment for each read/pair) should work, while some might need to be modified as required by downstream tools.

    mkdir mapping
    bowtie2 -q -p 6 -t --no-mixed --no-discordant --no-unal --dovetail -x refGenomes/D5 -1 <(zcat trimmed/D1H_1_val_1.fq.gz) -2 <(zcat trimmed/D1H_2_val_2.fq.gz) -S mapping/D1H.sam 2>mapping/D1H.log
       
    # save chromosome size
    head -50 mapping/D1H.sam|grep '@SQ'|awk -v OFS='\t' '{print $2,$3}'|sed 's/.N://g' > mappingF/chr.size.txt
     
* `-x refGenomes/D5`: use ref genome
* `-1 trimmed/D1H_1_val_1.fq`: paired end read 1
* `-2 trimmed/D1H_2_val_2.fq`: paired end read 2
* `-q`: takes fastq files
* `-p 6`: use 6 thread
* `-t`: Print the amount of wall-clock time taken by each phase.
* `--no-mixed --no-discordant`: discard discordant mapping
* `--no-unal`: Suppress SAM records for reads that failed to align.
* `--dovetail`: allow pair to overlap and over extend

Bowtie2 mapping results in `SAM` format were converted to `BAM` with a quality filter >=20.

    samtools view -q 20 -bS mapping/D1H.sam | samtools sort - -o mapping/D1H_q20.bam ; samtools index mapping/D1H_q20.bam

### 2.2 Mapping QC

[Deeptools](http://deeptools.readthedocs.io/en/latest/index.html) provide basic QC utilities:

* Assess sequencing depth - `plotCoverage` generates two plots. The first one simply represents the frequencies of the found read coverages. The second plot shows what is the fraction of the genome that has a certain depth.
* Check fragment size - `bamPEFragmentSize`.
* Assess the reproducibility of replicates by computing sample correlations - `multiBamSummary` followed by `plotCorrelation`, `plotPCA`.

08/07/18: Command can be ran on biocrunch2 and speedy, but not biocrunch:

    module load py-deeptools
    cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingD
    # input bams need to be sorted and index
    plotCoverage -b *q20.bam -o plotCoverage.pdf --ignoreDuplicates --minMappingQuality 20
    bamPEFragmentSize -b *q20.bam -hist plotPEFragmentSize.png
    
    multiBamSummary bins -b *q20.bam -out mapping.results.npz --outRawCounts mapping.results.tab
    head mapping.results.tab
    plotCorrelation -in mapping.results.npz -o plotHeatmap_pearson.pdf --corMethod pearson --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix plotHeatmap_pearson.tab.txt
    plotPCA -in mapping.results.npz -o plotPCA.pdf -T "PCA of mapping profiles"
    mkdir mappingQC
    mv plot* mappingQC/
    mv mapping* mappingQC/

If good correlations are seen between replicates, we will combind replicates in following analysis.

### 2.3 Convert BAM to BED and pool replicates

Quality filtered `BAM` results will be converted to `BED` format, e.g. 
   
    samtools sort -n D1H_q20.bam | bedtools bamtobed -i - -bedpe  | cut -f 1,2,6,7,8,9 | sort -T . -k1,1 -k2,2n -S 1G  > D1L_q20.bed

For both heavy and light digestions (e.g. `D1H_q20.bed`, `D2H_q20.bed`), a combo `BED` file (`DcH_q20.bed`) will be made by merging equal amount of mapped fragments from replicates. 

See wrapper r script [repBam2Bed.r](Scripts/repBam2Bed.r).

## 3.Differential nuclease sensitivity profiling analysis

### 3.0 Differential analysis of light and heavy nucleosomal profiles

Originally, each `BED` file (individual replicate, not pooled) was parsed into two fragment size ranges, 0-130bp and 130-260bp. Genome coverage was next calculated in RPM for Light, heavy and difference of Light-Heavy digestive conditions at each range, the resulted `BedGraph` files were compressed to `BigWig` files. The idea is that different sized fragments respond to digestion conditions differently, and it was suggested that small fragments by light digestion is comparable to DNase-seq profiles. I may further explore this direction later, based upon [dns.old.r](Scripts/dns.old.r).


### 3.1 Segmentation of sensitive and resistant fragments

The algorithm [iSeg](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2140-3) ([web server](http://lebesgue.fgcu.edu/iSeg)) identifies candidate segments from normalized differential genome coverage - [latest source v190207, v1.3.2](http://lebesgue.fgcu.edu/iSeg/downloads.html).

First, input files for iSeg need to be prepared from BED files, including individual and pooled runs, using script [iseg.pre.r](Scripts/iseg.pre.r). The coverage of Heavy and Light digested reads were estimated and quantile normalized across 20 bp windows tilling the genome, to generate DNS (D = L - H) profiles. Once BedGrapgh files (e.g. `A6D_q20_chr.size_w20_qnorm.bg`) are ready, inspect their MAD and SD values across genome. The BigWig files were inspected with IGV, and an issue of [mitochondria DNA insertion into D5 chr01](Discussion/mtInsertionD5.docx) was discussed.

Noting the differences in MAD between genomes (see "sumMDS.txt"), quantile normalization was performed for pooled D profiles ("sumMDS_quantile.txt"), and then subjected to iSeg analysis, using script [iseg.run.r](Scripts/iseg.run.r). The example commands of iSeg are:

    ~/iSegv190207/iSeg -cz -ctp -bc 1.0-2.0-3.0 -minwl 10 -maxwl 100 -d iseg/A6Dn_q20_chr.size_w20_qnorm_qnorm.bg -of A6Dn >isegv1.3.2/A6Dn.manifest.txt
    ~/iSegv190207/iSeg -cz -ctp -bc 1.0-2.0-3.0 -minwl 10 -maxwl 100 -d iseg/DcD_q20_chr.size_w20_qnorm_qnorm.bg -of DcD >isegv1.3.2/DcD.manifest.txt
    ~/iSegv190207/iSeg -cz -ctp -bc 1.0-2.0-3.0 -minwl 10 -maxwl 100 -d iseg/FcD_q20_chr.size_w20_qnorm_qnorm.bg -of FcD >isegv1.3.2/FcD.manifest.txt
    ~/iSegv190207/iSeg -cz -ctp -bc 1.0-2.0-3.0 -minwl 10 -maxwl 100 -d iseg/McD_q20_chr.size_w20_qnorm_qnorm.bg -of McD >isegv1.3.2/McD.manifest.txt

The raw segmentation output file contains 6 columns: chr, start, end, height, t-stat, p-value. Use script [iseg.summary.r](Scripts/iseg.summary.r) to generate some summary plots, including the number, size and distribution of detected segments, and output BED format. This can also be done using the iSeg python wrapper version.

**Summary**: IGV visualization of BED outputs revealed consistent segmentation profiles between v1.3.2 to v180809b; [iSegSummary021119.xlsx](iSegSummary021119.xlsx) showed that BC6 produced 1% of MSFs and MRFs each, will be used for following analysis. 
 
### 3.2 Genomic annotation of segments

Use ChIPseeker for peak annotation with [segAnno.r](Scripts/segAnno.r): 

1. visualizing segment locations on chromosomes
2. plot coverage profile and heatmap over TSS
3. annotate genomic location of segments with nearby genes.

### 3.3 Motif discovery and fuctional annotation of segements


## 4. Nucleosome positioning analysis

### 4.1 Nucleosome calling and classification

Bowtie2 mapping results in `BED` with a quality filter >=20, were first filtered to remove extra long reads over 260bp, and then trimmed fragments to 50 bp around nucleosome dyad. After read pre-processing, genome coverage was calculated in RPM for each digestive condition and each rep. Fast Fourier Transform (FFT) was applied to filter noise for coverage profile prior to peak detection of nucleosome positioning. Well and loosely positioned nucleosomes were characterized. [nucleR.r](Scripts/nucleR.r)

Nucleosome repeat length was calculated and compared between genomes. [NRL.r](Scripts/NRL.r)

What percentage of genome are covered by well-positioned and fuzzy nucleosomes? [nucleR.post.r](Scripts/nucleR.post.r)

### 4.2 Other potential analysis

[NucTools](https://homeveg.github.io/nuctools/) used replicates instead of peak shape by [nucleR](http://bioconductor.org/packages/release/bioc/html/nucleR.html) to distinguish stable and fuzzy nucleosomes. [DANPOS2](https://sites.google.com/site/danposdoc/) provides tool for calculating nucleosome occupancy change scores between conditions. [NUCwave](http://nucleosome.usal.es/nucwave/)is a wavelet-based bioinformatic tool that generates nucleosome occupation maps.

Analyses of dinucleotide frequency followed previously published methods (Locke et al., 2010; Valouev et al., 2011).

### 5. here is test
