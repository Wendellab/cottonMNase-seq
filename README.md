# Differential Nuclease Sensitivity Profiling of Chromatin in Diploid and Allopolyploid Cotton

Integated MNase-seq and RNA-seq analysis of duplicated gene expression and regulation in allopolyploid cotton.

## Project Description

To investigate differential regulatory control of duplicated genes (*aka* homoeologs) in allopolyploid cotton, we initiated a collaborative project in 2016 with Dr. Daniel L. Vera and Professor Hank W. Bass (Florida State University, USA) to explore the role of chromatin accessibility and nucleosome position on duplicated gene expression in allopolyploid cotton. Nucleosome formation is known to directly regulate the access of regulatory proteins to DNA sequences, and strongly associated with gene expression and other features of epigenetic modifications. Using a technique based on chromatin digestion by micrococcal nuclease followed by illumine sequencing (MNase-seq), we prepared mono-nucleosomal DNAs from four *Gossypium* species (both diploid parents, their F1 hybrid, and the natural allopolyploid cotton) using two different levels of MNase digestion.

### Timeline ###

* **2016**: A modified protocol was estblished to extract high-quality nucleus from mature cotton tissues including leaf and petal, which proves to be the most critical step of chromatin accessibility assays. 
* **2017**: Library preparation and sequencing was finished for 14 out of 16 samples. 
* **2018**: All sequencing datasets were completed. Data analysis in ongoing.
* **2019-2021**: Data analysis and manuscript writing.

## Data analysis workflow

Working directory `/work/LAS/jfw-lab/hugj2006/cottonLeaf`

Long-term storage `/lss/research/jfw-lab/Projects/MNase-seq/cottonLeaf/`

1. [Processing of differential MNase-seq datasets](DA_diffMNase-seq.md) resulted into the genome-wide characterization of nucleosome positioning, differential nuclease sensitivity (DNS), and sub-nucleosomal particle occupancy (SPO). Chromatin accessible regions were annotated by MNase sensitive and resistant footprints (MSFs/MRFs) and SPO fragment centers. 
2. [Processing of other datasets](DA_otherDatasets.md) including ATAC-seq, DNase-seq, Hi-C, ChIP-seq, etc.
3. [Analyze RNA-seq data](DA_RNAseq.md) to examine gene expression patterns regarding bias, dominance, *cis*/*trans* regulation, and partition of hybridization and genome doubling effects.
4. [Integrative analyses of above datasets driven by research questions](DA_questions.md)



**Misc** Other proposed tasks include:

- Assess the reproducibility of replicates: 

- homoeolog/ortholog group

- Bioinformatic prediction of cis elements and regions

- Other layering datasets, e.g. Hi-C, DNA fragility, etc.

  

**Discussion**

* [mitochondria DNA insertion into D5 chr01](Discussion/mtInsertionD5.docx)
* [Higher variability of genome-wide DNS in D/Dt vs A/Dt](Discussion/HankDiscussion2018.docx)

**Scripts archived but no longer included in the workflow**

- [dns.old.r](Scripts/dns.old.r) for original differential nuclease sensitivity profiling analysis on 0-130bp and 130-260bp size ranges. 

