# Differential Nuclease Sensitivity Profiling of Chromatin in Diploid and Allopolyploid construction

MNase-seq and RNA-seq analysis of duplicated gene expression and regulation in allopolyploid cotton

## Project description and progress

To investigate differential regulatory control of duplicated genes (*aka* homoeologs) in allopolyploid cotton, we initiated a collaborative project in 2016 with Dr. Daniel L. Vera and Professor Hank W. Bass (Florida State University, USA) to explore the role of chromatin accessibility and nucleosome position on duplicated gene expression in allopolyploid cotton. Nucleosome formation is known to directly regulate the access of regulatory proteins to DNA sequences, and strongly associated with gene expression and other features of epigenetic modifications. Using a technique based on chromatin digestion by micrococcal nuclease followed by illumine sequencing (MNase-seq), we prepared mono-nucleosomal DNAs from four *Gossypium* species (both diploid parents, their F1 hybrid, and the natural allopolyploid cotton) using two different levels of MNase digestion.

**2016**: A modified protocol was estblished to extract high-quality nucleus from mature cotton tissues including leaf and petal, which proves to be the most critical step of chromatin accessibility assays. 

**2017**: Library preparation and sequencing was finished for 14 out of 16 samples. 

**2018**: All sequencing datasets were completed. Data analysis in ongoing.

## Data analysis:

- [Summary](dataAnalysis.md)
- [MNase-seq processing](DA_diffMNase-seq.md)
- [RNA-seq processing](DA_RNAseq.md)
- [External datasets: ChIP-seq, DNase-seq, etc.](DA_otherDatasets.md)
- homoeolog/ortholog group
- Prediction of nucleosome positioning
- Bioinformatic prediction of cis elements and regions
- Other layering datasets, e.g. Hi-C, DNA fragility, etc.

## Scripts

- Quality trimming and sequence adaptor removal: [runTrimGalore.sh](runTrimGalore.sh)
- Bowtie2 mapping: [runBowtie2.D.sh](runBowtie2.D.sh) is an example for D-genome reads mapped against D-genome reference.
- Differential nuclease sensitivity profiling analysis: [dns.r](dns.r)
- Nucleosome calling and classification: [nucleR.r](nucleR.r) & [nucleR.source.r](nucleR.source.r)
- Assess the reproducibility of replicates: 


## Manuscript preparation

- [Introduction](introduction.md): to be written
- [Materials and Methods](methods.md): ongoing
- Results: [outline](results.md); [preliminary results](results_preliminary.md)
- Discussions:

