# Question-driven Integrative Analysis
---

**Regulatory landscape questions**

1. [How are various chromatin features compared with each other and in associated with gene expression?](#1-Inspect various chromatin features)
2. [Do different approaches - DNS, SPO, and ATAC provide a conserved view of chromatin assessiblity?](#2-Compare different profiles of chromatin accessibility)


**Polyploidy Questions** 

3. [How does allopolyploidy alter chromatin architecture, in terms of nucleosome origanization and chromatin accessibility?](#3-Compare chromatin profiles between diploid and polyploid cotton)
4. [How does the alternation of chromatin structure drive gene expression evolution?](#4-Associate chromatin profiles with duplicated gene expression patterns)

## 1. Inspect various chromatin features

Use [Qregulation1.r](scripts/Qregulation1.r) to conduct following analyses:

* Assemble various chromatin feature profiles that were preprocessed from [differential MNase-seq](DA_diffMNase-seq.md) (heavy and light nucleosome occupancy, DNS, SPO, nucleosome positioning) and [ATAC-seq](DA_otherDatasets.md)
* Plot Pearson/Spearman correlations and plot PCA
* Aggregate profiles on gene regions (scaled gene body, TSS, TES) relative to expression levels
* Aggregate profiles on Gypsy regions
* _TO-DO: Aggregate NP over 3 ACR categories for DNS and SPO each_.

## 2. Compare different profiles of chromatin accessibility

Use [Qregulation2.r](scripts/Qregulation2.r) to conduct following analyses in D5 genome:

* Compare characteristics (structure and annotation) of accesible chromatin regions (ACRs) resulted from DNS, SPO, and ATAC. 
* Aggregate profiles on each other ACR centers

## 3. Compare chromatin profiles between diploid and polyploid cotton

Examine hybridization and allopolyploidization effects: diploid A2/D5 vs F1 vs AD1.

0. TE composition [TE.r](scripts/TE.r): compare TE compositions, TE insertion relative to genes, Gypsy relative to more detailed genic structure。What percentages of the genome are covered by gene and TEs? What are the frequency of different TE categories in gene vicinity? What are the frequency of genes in TE vicinity? 
1. Nucleosome organization [nuc.post.r](scripts/nuc.post.r):  constrast the genome coverages, nucleosome repeat link lengths for observed and predicted nucleosome position (observed Well/Fuzzy categorization)
2. Accessible chromatin regions [ACR.r](scripts/ACR.r): marginal comparison of ACR categorization (DNS, SPO, DNase-seq);  ChIPseeker genomic annotation for pACR and gACR (DNS); Peak overlaps for A2 vs F1.At and D5 vs F1.Dt
3. Differential accessibility analysis [diffACRs_csaw_AD1ref.r](scripts/diffACRs_csaw_AD1ref.r) and [diffACRs_csaw_F1ref.r](scripts/diffACRs_csaw_F1ref.r) : diffbind-like vs csaw methods
4. 3D genome organization [hic.r](scripts/hic.r): 



## 4. Associate chromatin profiles with duplicated gene expression patterns

Use [Qregulation4.r](scripts/Qregulation4.r) to conduct following analyses:

* Examine promoter DNS profiles in association with B, Bp, Hr, Wr, and Pr

---



Aggregate profiles on gene regions

* Use [aggreNvisual.prep.r](scripts/aggreNvisual.prep.r) to address the following question:

  * What is the genomic composition of annotated gene and TE regions?
    * What percentages of the genome are covered by gene and TEs? What are the frequency of different TE categories in gene vicinity? What are the frequency of genes in TE vicinity? 
  * What are the patterns of nucleosome organization, nuclease sensitivity and chromatin openness associated with different genomic features?
    * H1: higher level of chromatin openess in gene promoter than other genomic regions.
    * H2: higher openess in gene promters depleted of TEs.
  * How do above patterns differ between A and D diploids, and between the subgenomes of A2xD5 F1 and the tetraploid AD?    

## 0. Prep and Tools

### 0.1 Prepare chromatin structure profile data

Use [aggreNvisual.prep.r](scripts/aggreNvisual.prep.r) to prepare 

Chromatin profiles that were preprocessed from [differential MNase-seq](DA_diffMNase-seq.md), [ATAC-seq](DA_otherDatasets.md) , including:

1. nucleosome occupancy profiles from replicated heavy and light digestion, full-ranged fragment size
2. DNS profiles from normalized light minus heavy nucleosome occupancy profiles, full-ranged fragment size
3. nucleosome positioning profiles from heavy digestion only, <200 bp fragments trimmed to 50 bp around dyad
4. Size partitioned nucleosome occupancy profiles of heavy, light and DNS, range **0-130bp as sub-nucleosomal fragments** and 130-260 bp.
5. Frenters of MOA-seq (MNase Open and Accessible chromatin mapping) by ~50-130 bp fragments from light digestion, 

The RPM (reads per million) and quantile normalization in original DNS pipeline bothered me, given the 3-fold genome size difference between species surveyed. For comparative aggragation over genomic features, I generated new BigWig files from BAM using Deeptools bamCoverage first using RPGC (reads per 1x sequenced genomic content) normalization, and then switched to RPKM for simplicity, as explained [here](https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html ).

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


### 0.3 Tools

Many tools are available and reviewed in a biostars thread [here](https://www.biostars.org/p/180314/):

* [Genomation](https://bioconductor.org/packages/release/bioc/html/genomation.html) - **R package used here**.
* [Deeptools](http://deeptools.readthedocs.io/en/latest/index.html) - generally faster than R.
* [SeqPlots](http://przemol.github.io/seqplots/) - corresponding genomic packages are required from bioconductor, useless for cotton, shame
* [ngsplot](https://github.com/shenlab-sinai/ngsplot) - unix stripts calling R and python usage; pretty figure but need to build genome

## 1. Aggregation and visualization of chromatin structure profiles

Use [aggreNvisual.prep.r](scripts/aggreNvisual.prep.r) to address the following question:

* What is the genomic composition of annotated gene and TE regions?
    * What percentages of the genome are covered by gene and TEs? What are the frequency of different TE categories in gene vicinity? What are the frequency of genes in TE vicinity? 
* What are the patterns of nucleosome organization, nuclease sensitivity and chromatin openness associated with different genomic features?
    * H1: higher level of chromatin openess in gene promoter than other genomic regions.
    * H2: higher openess in gene promters depleted of TEs.
* How do above patterns differ between A and D diploids, and between the subgenomes of A2xD5 F1 and the tetraploid AD?    

## 2. Comparison of chromatin profiles in association with duplicated gene expression patterns

[expNuc.r](scripts/expNuc.r)

---

Make plots for visualization. Generate heatmap of read coverages, for visualizing the read coverages for genomic regions. The default setting plots profile on top of heatmaps.

    computeMatrix reference-point -S *q20_unified.bw -R D5.gene.bed -o TSS.gz --referencePoint TSS -b 1500 -a 1500 --skipZeros 
    plotHeatmap -m TSS.gz -out plotHeatmap_TSS.png
    
    computeMatrix reference-point -S *q20_unified.bw -R D5.gene.bed -o TES.gz --referencePoint TES -b 1500 -a 1500 --skipZeros 
    plotHeatmap -m TES.gz -out plotHeatmap_TES.png


​    
​    computeMatrix scale-regions -S *q20_unified.bw -R D5.gene.bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1500 -a 1500 --skipZeros
​    plotHeatmap -m scaleGeneBody.gz -out plotHeatmap_scaleGeneBody.png

In addition, generate summary plot of "meta-profile", for visualizing the average read coverages over a group of genomic regions, and output `.tab` file for R analysis.

    plotProfile -m scaleGeneBody.gz --perGroup -out plotProfileGroup.png --outFileNameData plotProfile.tab


​    

