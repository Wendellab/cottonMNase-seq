# Data analysis of cotton MNase-seq datasets
---

Outline: 

* Process reads: quality filter and trimming
* Map reads: mapping rates, check consistency acorss libraries
* Generate and compare genomic coverage profiles: nucleosome occupancy (heavy & light), nuclease sensitivity (Light-heavy), nucleosome dycad positioning
* Peak calling from nuclease sensitivity profile to infer potential regulatory regions
* Investigate the genomic/genetic location of peaks, correlation with expression, etc. 

## 1. Preprocessing of DNA sequencing datasets

### MNase-seq datasets
Short-read data will be deposited in the NCBI short read archive ([SRP??????](http://trace.ddbj.nig.ac.jp/DRASearch/study?acc=SRP??????)), also as Biobroject [PRJNA??????](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA??????).

Paired-end reads of 16 samples: 4 species (A2, D5, A2xD5, Maxxa) X 2 digestive conditions (Heavy and Light) X 2 technical reps.

    cd cottonLeaf/rawfastq
    ln -s ~/jfw-lab/RawData/HGJ_leafMNase-seq/WTNHHW163125/Feb2017/*gz .
    ln -s ~/jfw-lab/RawData/HGJ_leafMNase-seq/WTNHHW163125/Feb2018/*gz .

#### Checking read quality with [FastQC](http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/)
    cd ..
    module load fastqc/0.11.3
    mkdir QCreport
    mkdir QCreport/raw
    mkdir QCreport/trimmed
    fastqc -o QCreport/raw rawfastq/*gz
    
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

### Cotton reference genomes
[CottonGen](https://www.cottongen.org/data/download/genome#Ass) compiles all published cotton genomes, and I will need 4 different reference genomes for AD1, A2, D5 and A2xD5. An improved AD1 reference became available on [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Ghirsutum_er).

    mkdir refGenomes
    cd refGenomes
    module load bowtie2
    ### D5_JGI - Paterson et al. 2012 Nature
    ln -s ~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/Dgenome2_13.fasta
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
    
    ### AD1_458 - Saski et al, 2017 
    ln -s ~/jfw-lab/GenomicResources/archived_resources/AD1Saski/v1.1/assembly/Ghirsutum_458_v1.0.fa.gz
    zcat Ghirsutum_458_v1.0.fa.gz |grep -n '>scaffold'|head -10 
    # 27134894:>scaffold_27
    zcat Ghirsutum_458_v1.0.fa.gz |head -27134893 >TM1new_26.fasta 
  
    ### A2Du2018 - Du et al, 2018 
    ln -s ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Du2018/1.PacBio-Gar-Assembly-v1.0/G.arboreum.Chr.v1.0.fasta.gz
    zcat G.arboreum.Chr.v1.0.fasta.gz |grep -n '>tig'|head -10 
    # 40:>tig00000043
    zcat G.arboreum.Chr.v1.0.fasta.gz |head -39 >A2Du_26.fasta 
    
    # build bowtie2 ref
    bowtie2-build TM1_26.fasta TM1
    bowtie2-build TM1new_26.fasta TM1new
    bowtie2-build F1_26.fasta F1
    bowtie2-build A2genome_13.fasta A2
    bowtie2-build Dgenome2_13.fasta D5
    bowtie2-build A2Du_26.fasta A2Du

## 2.Read mapping and calling of hypersensive sites

### Bowtie2 mapping
Default setting `-k 1` report 1 alignment for each read/pair) should work, while some might need to be modified as required by downstream tools.

    mkdir mapping
    bowtie2 -q -p 6 -t --no-mixed --no-discordant --no-unal --dovetail -x refGenomes/D5 -1 <(zcat trimmed/D1H_1_val_1.fq.gz) -2 <(zcat trimmed/D1H_2_val_2.fq.gz) -S mapping/D1H.sam 2>mapping/D1H.log
    samtools view -bS D1H.sam | samtools sort - -o D1H.sort.bam ; samtools index D1H.sort.bam

* `-x refGenomes/D5`: use ref genome
* `-1 trimmed/D1H_1_val_1.fq`: paired end read 1
* `-2 trimmed/D1H_2_val_2.fq`: paired end read 2
* `-q`: takes fastq files
* `-p 6`: use 6 thread
* `-t`: Print the amount of wall-clock time taken by each phase.
* `--no-mixed --no-discordant`: discard discordant mapping
* `--no-unal`: Suppress SAM records for reads that failed to align.
* `--dovetail`: allow pair to overlap and over extend

### Differential nuclease sensitivity profiling analysis
Bowtie2 mapping results in `SAM` format were converted to `BAM` with a quality filter >=20, and then further converted to `BED` format. Each `BED` file was parsed into two fragment size ranges, 0-130bp and 130-260bp. Genome coverage was next calculated in RPM for Light, heavy and difference of Light-Heavy digestive conditions at each range, the resulted `BedGraph` files were compressed to `BigWig` files.

See pipeline scripts in [dns.r](dns.r).

Taking sample **A6** for example, by inputing mapping results of its heavy and light MNase-seq data - `A6H.sam` & `A6L.sam`, output results include:

* Reference chromosome sizes: `chr.size.txt`
* Density plot of mapped fragment sizes: `sizeDistribution.pdf`
* BAM files: `A6H_q20.bam` & `A6L_q20.bam`
* BED files: `A6H_q20.bed` & `A6L_q20.bed`
* BED files parsed by size range: `A6H_q20_000-130.bed`, `A6H_q20_130-260.bed`, `A6L_q20_000-130.bed`, `A6L_q20_130-260.bed`
* BedGrapgh files of genome coverage: `A6H_q20_unified.bg`,`A6L_q20_unified.bg`,`A6D_q20_unified.bg`; `A6H_q20_000-130_unified.bg`, `A6L_q20_000-130_unified.bg`, `A6D_q20_000-130_unified.bg`;  `A6H_q20_130-260_unified.bg`,  `A6L_q20_130-260_unified.bg`, `A6D_q20_130-260_unified.bg`. 
* BigWig files:  `A6H_q20_unified.bw`, `A6L_q20_unified.bw`, `A6D_q20.bw_unified`; `A6H_q20_000-130_unified.bw`, `A6L_q20_000-130_unified.bw`, `A6D_q20_000-130_unified.bw`;  `A6H_q20_130-260_unified.bw`,  `A6L_q20_130-260_unified.bw`, `A6D_q20_130-260_unified.bw`. 


### Segmentation of sensitive and resistent fragments

The algorithm [iSeg](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2140-3) identifies candidate segments from normlized differential genome coverage.

See pipeline scripts in [iseg.r](iseg.r). Output results include:

* BedGrapgh files of 20bp sliding window coverage: `A6H_q20_chr.size_w20.bg`, `A6L_q20_chr.size_w20.bg`.
* Quantile normalized BedGraph: `A6H_q20_chr.size_w20_qnorm.bg`, `A6L_q20_chr.size_w20_qnorm.bg`,
`A6D_q20_chr.size_w20_qnorm.bg`.

Convert BedGraph to gff file
    
    awk '{print $1,$2,$3,".",".",$4}' OFS='\t' A6D_q20_chr.size_w20_qnorm.bg > A6D_q20_chr.size_w20_qnorm.gff

Run iSeg:

    ~/iseg -p MeanCut.csv -fdr 0.05 -bc 3.0 -minwl 10 -maxwl 100 -sig 1e-3 -d D2D_Chr13.bg -of D2D_Chr13_x0_bc3.0.txt -mad 0.10657938 -sd 0.170028794833 -cz -ctp

    ## Options:
    # -p: file for test statistic cutoff. default: MeanCut.csv.
    # -maxwl: Maximum Window length, the maximum size of a detected segment.
    # -minwl: Minimum Window length, the minimum size of a detected segment.
    # -bc: biological cutoff. For instance, a segment is significant only if the segment mean exceeds twice the standard deviation of the segment when set to 2.
    # -sig: An intial cut-off to determine significance. Usually, 0.001 or 0.0001 is sufficient.
    # -fdr: false discovery rate (FDR); determines the cut-off to call a significant region using Benjamini-Hochberg procedure.
    # -mad:
    # -sd:
    # -cz -ctp: true z-stat and t-stat; must given, otherwise found NO significant regions
    # -d: data file. no default.
    # -tf: file containing actual significant regions. no default.
    # -of: output file name. dedault: out.txt.


### Nucleosome calling and classification
Bowtie2 mapping results in `BED` with a quality filter >=20, were first filtered to remove extra long reads over 260bp, and then trimmed fragments to 50 bp around nucleosome dyad. After read pre-processing, genome coverage was calculated in RPM for each digestive condition and each rep. Fast Fourier Transform (FFT) was applied to filter noise for coverage profile prior to peak detection of nucleosome positioning. Well and loosely positioned nucleosomes were characterized.

See pipeline scripts in [nucleR.r](nucleR.r).

Taking sample **A6H** for example, by inputing quality filtered mapping results  - `A6H_q20.bed`, output results include:
* Plot mapping depth by chromosomeL: `A6H_q20.checkRawDepth.pdf`
* BigWig file: `A6H_q20.coverage.bw`
* GFF file of nucleosome position: `A6H_q20.nucleosome.gff`

## Downsteam analysis

Visualization and exploration of sequencing tracks (`BigWig`) around genomic features (e.g. TSS) is the key of downstream analysis. Many tools are availble and reviewed in a biostars thread [here](https://www.biostars.org/p/180314/):

* [Deeptools](http://deeptools.readthedocs.io/en/latest/index.html) - faster than R packages, waiting for bioIT to install
* [SeqPlots](http://przemol.github.io/seqplots/) - corresponding genomic packages are required from bioconductor, useless for cotton, shame
* [Genomation](https://bioconductor.org/packages/release/bioc/html/genomation.html) - R package
* [ngsplot](https://github.com/shenlab-sinai/ngsplot)

### Assess the reproducibility of replicates
Using coverages before size partition, check if Replicates correlate better than non-replicates. Using [Deeptools](http://deeptools.readthedocs.io/en/latest/index.html) for visualization of sample correlations based on the output of `multiBamSummary` or `multiBigwigSummary`. Pearson or Spearman methods are available to compute correlation coefficients. 

    multiBigwigSummary bins -b *q20_unified.bw -out results.npz
    plotCorrelation -in results.npz -o plotHeatmap_pearson.pdf --corMethod pearson --skipZeros --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix plotHeatmap_pearson.tab.txt

If good correlations are seen between replicates, we can consider pooling reps in following analysis.

### Check sequencing coverage, fragment size, etc.
Two plots were generated. The first one simply represents the frequencies of the found read coverages, which helps you judge how relevant the mean coverage value (printed next to the sample name) is. If the distribution of read coverages is more or less homoskedatic and, ideally, normally distributed (most likely it won’t be), then the mean is a very appropriate proxy for sequencing depth. The second plot shows what is the fraction of the genome that has a certain depth. 

    plotCoverage -b *sort.bam -o plotCoverage.pdf --ignoreDuplicates --minMappingQuality 20

For paired-end samples, we often additionally check whether the fragment sizes are more or less what we would expected based on the library preparation.

    bamPEFragmentSize -b *sort.bam -hist plotPEFragmentSize.png
    
### Visualize nucleosome occupancy coverage on genomic features

First, need to prepare genomic feature `bed` file:

    module load bedops
    
    cd mappingD
    head ~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/Dgenome2_13.gene.gff
    gff2bed < ~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/Dgenome2_13.gene.gff >D5.gene.bed
    cd ..
    
    cd mappingM
    head ~/jfw-lab/GenomicResources/archived_resources/gmapdb/AD1TM1/Gossypium_hirsutum_v1.1.gene.gff3
    gff2bed < <(grep 'gene' ~/jfw-lab/GenomicResources/archived_resources/gmapdb/AD1TM1/Gossypium_hirsutum_v1.1.gene.gff3) > TM1.gene.bed
    grep '^A' TM1.gene.bed >TM1.gene.A.bed
    grep '^D' TM1.gene.bed >TM1.gene.D.bed
    cd ..
    
    cd mappingA
    head ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Li/A2Li.exons.gff
    gff2bed < <(grep "mRNA" ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Li/A2Li.gene.gff) > A2.gene.bed
    cd ..
    
    cd mappingF
    cat mappingA/A2.gene.bed <(sed 's/^Chr/D5_chr/g' mappingDref/mappingD/D5.gene.bed) >F.gene.bed
    grep '^A' F.gene.bed >F.gene.A.bed
    grep '^D' F.gene.bed >F.gene.D.bed
    cd ..
    
Make plots for visualization. Generate heatmap of read coverages, for	visualizing the read coverages for genomic regions. The default setting plots profile on top of heatmaps.

    computeMatrix reference-point -S *q20_unified.bw -R D5.gene.bed -o TSS.gz --referencePoint TSS -b 1500 -a 1500 --skipZeros 
    plotHeatmap -m TSS.gz -out plotHeatmap_TSS.png
    
    computeMatrix reference-point -S *q20_unified.bw -R D5.gene.bed -o TES.gz --referencePoint TES -b 1500 -a 1500 --skipZeros 
    plotHeatmap -m TES.gz -out plotHeatmap_TES.png
    
    
    computeMatrix scale-regions -S *q20_unified.bw -R D5.gene.bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1500 -a 1500 --skipZeros
    plotHeatmap -m scaleGeneBody.gz -out plotHeatmap_scaleGeneBody.png
    
In addition, generate summary plot of "meta-profile", for visualizing the average read coverages over a group of genomic regions, and output `.tab` file for R analysis.
    plotProfile -m scaleGeneBody.gz --perGroup -out plotProfileGroup.png --outFileNameData plotProfile.tab

### bookmark: below to be sorted 
---

### Differential nucleosome occupancy analysis - [DANPOS2](https://sites.google.com/site/danposdoc/)

### Profiling Nucleopostioning - [nucleR](http://bioconductor.org/packages/release/bioc/html/nucleR.html) & [NUCwave](http://nucleosome.usal.es/nucwave/)
Briefly, we need to examine the wave-length patterns of mapping coverage on the reference genome, and then specifically locate peaks that fit the description of nucleosomes - proper length, non-overlapping, etc.

(Zhang et al. 2015): "Well-positioned and loosely positioned nucleosomes were identified using nucleR (Flores and Orozco, 2011). ... We used the filterFFT function of nucleR to remove noise and smooth the read count score of each position along chromosomes with the parameter pcKeepComp = 0.01. After noise removal, nucleosome peaks and centers/dyads were determined using the peakDetection function (threshold = 25%, score = true, width = 140). Overlapped peaks were merged into longer regions, which were defined as loosely positioned nucleosomes, and distinct individual peaks were defined as well-positioned nucleosomes. If the length of merged peaks is longer than 150 bp, this region is considered to contain more than two nucleosome dyads and thus, contains loosely positioned nucleosomes. If the length of merged peaks is shorter than 150 bp, this region is considered to contain a well-positioned nucleosome."

The phasogram and average distance between two adjacent nucleosomes were calculated using our previously reported methods (Zhang et al., 2013). The nucleosome occupancy change scores were calculated by DANPOS (Chen et al., 2013). Analyses of dinucleotide frequency followed previously published methods (Locke et al., 2010; Valouev et al., 2011).

## Visualization and other result presentation
