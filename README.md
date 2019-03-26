# Differential Nuclease Sensitivity Profiling of Chromatin in Diploid and Allopolyploid Cotton

MNase-seq and RNA-seq analysis of duplicated gene expression and regulation in allopolyploid cotton.

## Project description and progress

To investigate differential regulatory control of duplicated genes (*aka* homoeologs) in allopolyploid cotton, we initiated a collaborative project in 2016 with Dr. Daniel L. Vera and Professor Hank W. Bass (Florida State University, USA) to explore the role of chromatin accessibility and nucleosome position on duplicated gene expression in allopolyploid cotton. Nucleosome formation is known to directly regulate the access of regulatory proteins to DNA sequences, and strongly associated with gene expression and other features of epigenetic modifications. Using a technique based on chromatin digestion by micrococcal nuclease followed by illumine sequencing (MNase-seq), we prepared mono-nucleosomal DNAs from four *Gossypium* species (both diploid parents, their F1 hybrid, and the natural allopolyploid cotton) using two different levels of MNase digestion.

### Timeline ###

* **2016**: A modified protocol was estblished to extract high-quality nucleus from mature cotton tissues including leaf and petal, which proves to be the most critical step of chromatin accessibility assays. 
* **2017**: Library preparation and sequencing was finished for 14 out of 16 samples. 
* **2018**: All sequencing datasets were completed. Data analysis in ongoing.
* **2019**: Data analysis and manuscript writing.

## Data analysis pipelines

**First**, [Processing of differential MNase-seq datasets](DA_diffMNase-seq.md) resulted into the genome-wide detection of nucleosome position, MNase sensitive and resistant footprints (MSFs/MRFs). Scripts includes:

- [runTrimGalore.sh](Scripts/runTrimGalore.sh) for quality trimming and sequence adaptors removal.
- [runBowtie2.D.sh](Scripts/runBowtie2.D.sh) is an example for Bowtie2 mapping against the D-genome reference, which needs to be moditified for different reference genome used.
- [repBam2Bed.r](Scripts/repBam2Bed.r) converted mapping results BAM to BED files and pooled replicates.
- [iseg.pre.r](Scripts/iseg.pre.r), [iseg.run.r](Scripts/iseg.run.r) and [iseg.summary.r](Scripts/iseg.summary.r) performed L-H and iSeg analysis and summary.
- [makeTxDb.r](Scripts/makeTxDb.r) prepared `txdb.<refGenome>.sqlite` for handling gene anonation with R
- [segAnno.r](Scripts/segAnno.r) conducted genomic annotation for MSFs and MRFs with ChIPseeker.
- [nucleR.r](Scripts/nucleR.r), [NRL.r](Scripts/NRL.r), and [nucleR.post.r](Scripts/nucleR.post.r) performed nucleosome positioning analyses.


**Second**, [Processing of RNA-seq data](DA_RNAseq.md) provides expression profiles coupling with the nucleosome positioning and open chromatin profiles.

**Third**, [Processing of other data]()

**Third**, [Question driven integrative analysis](DA_questions.md)


- [expNuc.r](expNuc.r)

Third, [Integrated analysis]() 

1. inspect H, L and D profiles around TSS and TTS in association with expression quantiles in each species
2. Inspect 
inspect the presence of MSFs in 1kb promoter region, 
[]()




- Nucleosome calling and classification: [nucleR.r](nucleR.r) & [nucleR.source.r](nucleR.source.r)
- Assess the reproducibility of replicates: 


S
- [RNA-seq processing](DA_RNAseq.md)
- [External datasets: ChIP-seq, DNase-seq, etc.](DA_otherDatasets.md)
- homoeolog/ortholog group
- Prediction of nucleosome positioning
- Bioinformatic prediction of cis elements and regions
- Other layering datasets, e.g. Hi-C, DNA fragility, etc.

Following issues were discussed:

* [mitochondria DNA insertion into D5 chr01](Discussion/mtInsertionD5.docx)
* [Higher variability of genome-wide DNS in D/Dt vs A/Dt](Discussion/HankDiscussion2018.docx)

**Scripts archived but no longer included in pipeline**

- [dns.old.r](Scripts/dns.old.r) for original differential nuclease sensitivity profiling analysis on 0-130bp and 130-260bp size ranges. 


Working directory - `/work/LAS/jfw-lab/hugj2006/cottonLeaf`
Long-term storage - `/lss/research/jfw-lab/Projects/MNase-seq/cottonLeaf/'


## Manuscript preparation

- [Introduction](introduction.md): to be written
- [Materials and Methods](methods.md): ongoing
- Results: [outline](results.md); [preliminary results](results_preliminary.md)
- Discussions:

