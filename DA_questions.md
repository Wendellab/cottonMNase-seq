# Question driven integrative analysis
---

**Outline Questions** 

1. [How does allopolyploidy alter chromatin architecture?]()
  - [Characterize and compare chromatin structure profiles over]
  - dd
2. [How does the alternation of chromatin structure drive gene expression evolution?]()


## 1. Aggregation and visualization of chromatin structure profiles on genome

Check [aggreNvisual.r](scripts/aggreNvisual.r) for R scripts covering below tasks.

### 1.1 Prepare chromatin structure profile data

Chromatin structure profiles of nucleosome occupancy and open chromatin were generated in previous analysis of [Processing of differential MNase-seq datasets](DA_diffMNase-seq.md), including:

1. nucleosome occupancy profiles from replicated heavy and light digestion, full-ranged fragment size
2. DNS profiles from normalized light minus heavy nucleosome occupancy profiles, full-ranged fragment size
3. nucleosome positioning profiles from heavy digestion only, <200 bp fragments trimmed to 50 bp around dyad
4. Size partitioned nucleosome occupancy profiles of heavy, light and DNS, range **0-130bp as sub-nucleosomal fragments** and 130-260 bp.
5. 

The RPM(reads per million) and quantile normalization in DNS pipeline bothered me, given the 3-fold genome size difference between species surveyed. For comparative aggragation over genomic features, I generated new BigWig files from BAM using Deeptools bamCoverage first with RPGC (reads per 1x sequenced genomic content) normalization, and then switched to RPKM, as explained [here](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html ).

RPGC = reads per genomic content (1x normalization); scaling factor = (total number of mapped reads * fragment length) / effective genome size. 

Effective genomic size is the part of genome that is uniquely mappable, according to [General-deepTools-FAQs](https://github.com/deeptools/deepTools/wiki/General-deepTools-FAQs). I first tried `faCount` to report non-'N' size, but little Ns was found in A2Du, not right.

    ~/kent/faCount /work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/A2Du_13.fasta >/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/A2Du_13.facount.txt
    ~/kent/faCount /work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/Dgenome2_13.fasta >/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/Dgenome2_13.facount.txt
    ~/kent/faCount /work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/TM1new_26.fasta >/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/TM1new_26.facount.txt

The BEDtool2 `genomecov -d` can be used to calculate the coverage per site in the genome, and regions with mapped reads have depth >0.

    module load bedtools2
    bedtools genomecov -ibam D1H_q20.bam -d |awk '$3>0 {print $3}' |wc -l
    # 704408202
    
Using the obtained effective size to convert BAM to BigWig

    bamCoverage -b D1H_q20.bam -o F2L_q20_142-200.bw --binSize 20 --normalizeTo1x 704408202
    bamCoverage -b D1H_q20.bam -o F2L_q20_142-200.bw --binSize 20 --minFragmentLength 130 --maxFragmentLength 260 --normalizeTo1x 1878118078bamCoverage -b D1H_q20.bam -o F2L_q20_142-200.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 130 --normalizeTo1x 704408202
    bamCoverage -b D1H_q20.bam -o F2L_q20_142-200.bw --binSize 20 --minFragmentLength 130 --maxFragmentLength 260 --normalizeTo1x 704408202


### 1.2 Prepare genomic annotation
 
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


### 1.3 Aggregation and visualization
    
Many tools are available and reviewed in a biostars thread [here](https://www.biostars.org/p/180314/):

* [Genomation](https://bioconductor.org/packages/release/bioc/html/genomation.html) - **R package used here**.
* [Deeptools](http://deeptools.readthedocs.io/en/latest/index.html) - faster than R packages
* [SeqPlots](http://przemol.github.io/seqplots/) - corresponding genomic packages are required from bioconductor, useless for cotton, shame
* [ngsplot](https://github.com/shenlab-sinai/ngsplot) - unix stripts calling R and python usage


## 2. Comparison of chromatin profiles in association with duplicated gene expression patterns

[expNuc.r](scripts/expNuc.r)

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


    

