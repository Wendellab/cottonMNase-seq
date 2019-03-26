# Question driven integrative analysis
---

**Outline Questions** 

1. [What are the peak distributions of variable chromatin structure profiles over different genomic features?](#)




## 1. Aggregation and visualization of chromatin structure profile data

Check [aggreNvisual.r](scripts/aggreNvisual.r) for R scripts covering below tasks.

### 5.1 Prepare chromatin structure profile data

Multiple chromatin structure profiles were generated regarding nucleosome occupancy and open chromatin, including:

1. nucleosome occupancy profiles from replicated heavy and light digestion, full-ranged fragment size
2. DNS profiles from normalized light minus heavy nucleosome occupancy profiles, full-ranged fragment size
3. nucleosome positioning profiles from heavy digestion only, <200 bp fragments trimmed to 50 bp around dyad
4. Size partitioned nucleosome occupancy profiles of heavy, light and DNS, range 0-130bp as sub-ucleosomal fragments and 130-260 bp.

Although these profiles were generated in DNS and nucleosome positioning steps, their different normalization methods bother me. I would rather re-make these BigWig files from BAM using Deeptools bamCoverage with RPGC (reads per 1x sequenced genomic content) normalization, as explained [here](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html ).


CPM = Counts Per Million mapped reads, same as CPM in RNA-seq; 
BPM = Bins Per Million mapped reads, same as TPM in RNA-seq; 
RPGC = reads per genomic content (1x normalization); Mapped reads are considered after blacklist filtering (if applied). 


RPKM (per bin) = number of reads per bin / (number of mapped reads (in millions) * bin length (kb)). 
CPM (per bin) = number of reads per bin / number of mapped reads (in millions). 
BPM (per bin) = number of reads per bin / sum of all reads per bin (in millions). 
RPGC (per bin) = number of reads per bin / scaling factor for 1x average coverage. 
None = the default and equivalent to not setting this option at all. 


scaling factor = (total number of mapped reads * fragment length) / effective genome size. 

The scaling factor used is the inverse of the sequencing depth computed for the sample to match the 1x coverage. 

This option requires –effectiveGenomeSize. Each read is considered independently, if you want to only count one mate from a pair in paired-end data, then use the –samFlagInclude/–samFlagExclude options.


### 5.2 Prepare genomic annotation
 
Visualization and exploration of sequencing tracks (`BigWig`) around genomic features (e.g. TSS) is the key of downstream analysis. Annotation of genomic features need to be prepared in BED format. BEDOPS tool [gff2bed](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gff2bed.html) can be used to convert gff to bed for deeptools.

    # 08/07/18 module bedops no longer supported, run gff2bed locally
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
    
    cd mappingMnew
    zcat ~/jfw-lab/GenomicResources/archived_resources/AD1Saski/annotation/Ghirsutum_458_v1.1.gene_exons.gff3.gz |head
    
    cd mappingA
    head ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Li/A2Li.exons.gff
    gff2bed < <(grep "mRNA" ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Li/A2Li.gene.gff) > A2.gene.bed
    cd ..
 
    cd mappingAnew
    head ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Du2018/1.PacBio-Gar-Assembly-v1.0/G.arboreum.Chr.v1.0.gff
    gff2bed < <(grep "mRNA" G.arboreum.Chr.v1.0.gff) > A2.gene.bed
    cd ..
    
    cd mappingF
    cat mappingA/A2.gene.bed <(sed 's/^Chr/D5_chr/g' mappingDref/mappingD/D5.gene.bed) >F.gene.bed
    grep '^A' F.gene.bed >F.gene.A.bed
    grep '^D' F.gene.bed >F.gene.D.bed
    cd ..

Alternatively, load txDB in R.

    library(GenomicFeatures)
    library(genomation)
    txdb <- loadDb("refGenomes/txdb.TM1saski.sqlite")
    TSS = promoters(txdb, upstream=1000, downstream=1000)
    bw.files<-list.files("mappingMnew", pattern="..D.*w20_qnorm.bw",full.names=TRUE)
    sml = ScoreMatrixList(target = bw.files, windows = TSS, type="bigWig",strand.aware=TRUE)


### 5.3 Aggregation and visualization
    
Many tools are availble and reviewed in a biostars thread [here](https://www.biostars.org/p/180314/):

* [Deeptools](http://deeptools.readthedocs.io/en/latest/index.html) - faster than R packages, waiting for bioIT to install
* [SeqPlots](http://przemol.github.io/seqplots/) - corresponding genomic packages are required from bioconductor, useless for cotton, shame
* [Genomation](https://bioconductor.org/packages/release/bioc/html/genomation.html) - R package
* [ngsplot](https://github.com/shenlab-sinai/ngsplot) - unix stripts calling R and python usage

---

Make plots for visualization. Generate heatmap of read coverages, for visualizing the read coverages for genomic regions. The default setting plots profile on top of heatmaps.

    computeMatrix reference-point -S *q20_unified.bw -R D5.gene.bed -o TSS.gz --referencePoint TSS -b 1500 -a 1500 --skipZeros 
    plotHeatmap -m TSS.gz -out plotHeatmap_TSS.png
    
    computeMatrix reference-point -S *q20_unified.bw -R D5.gene.bed -o TES.gz --referencePoint TES -b 1500 -a 1500 --skipZeros 
    plotHeatmap -m TES.gz -out plotHeatmap_TES.png
    
    
    computeMatrix scale-regions -S *q20_unified.bw -R D5.gene.bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1500 -a 1500 --skipZeros
    plotHeatmap -m scaleGeneBody.gz -out plotHeatmap_scaleGeneBody.png
    
In addition, generate summary plot of "meta-profile", for visualizing the average read coverages over a group of genomic regions, and output `.tab` file for R analysis.

    plotProfile -m scaleGeneBody.gz --perGroup -out plotProfileGroup.png --outFileNameData plotProfile.tab


    

