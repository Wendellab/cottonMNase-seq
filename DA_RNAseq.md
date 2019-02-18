# RNA-seq analysis corresponding to cotton MNase-seq datasets
---

## Preprocessing and mapping of RNA-seq datasets

### Locare FASTQ files
For both MNase-seq and RNA-seq experiments, mature leaf tissue was harvested from flowering branches at 5pm, and immediately flash frozen in liquid nitrogen and stored at −80°C. Four *Gossypium* accessions were used including a natural allopolyploid, *G. hirsutum* (AD1) cultivar Acala Maxxa, and models of its A- and D-genome diploid progenitors, *G. arboreum* (A2) and *G. raimondii* (D5), as well as the corresponding interspecific diploid F1 hybrid (A2 × D5).
    
    cd /home/hugj2006/jfw-lab/Projects/MNase-seq/cottonLeaf/RNAseq/rawfastq
    ln -s ~/jfw-lab/RawData/Duplicated_Networks/transcriptomic/SD5-Maxxa* .
    ln -s ~/jfw-lab/RawData/Duplicated_Networks/transcriptomic/SD5-D5* .
    ln -s ~/jfw-lab/RawData/Duplicated_Networks/transcriptomic/SD5-A2* .
    ls

### Quality trimming

    bash runTrimGalore.sh >runTrimGalore.071117.txt 2>&1
    grep -E "Total reads|Reads written" trimmed/*txt

### Cotton reference transcriptomes
[CottonGen](https://www.cottongen.org/data/download/genome#Ass) compiles all published cotton genomes, and I will need 4 different reference genomes for AD1, A2, D5 and A2xD5. An improved AD1 reference became available on [Phytozome](https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Ghirsutum_er).

    mkdir refTranscriptomes
    cd refTranscriptomes
    module load bowtie2
    
    ### D5_JGI - Paterson et al. 2012 Nature
    cp /lss/research/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/Dgenome2_13.cds.fasta >D5.transcripts.fasta
    grep -c '>' D5.transcripts.fasta
    # 37223
    # or see: https://github.com/huguanjing/homoeologGeneExpression-Coexpression/blob/master/bowtie2hylite.sh
        
    ### AD1_458 - Saski et al, 2017
    zcat /lss/research/jfw-lab/GenomicResources/archived_resources/AD1Saski/annotation/Ghirsutum_458_v1.1.transcript_primaryTranscriptOnly.fa.gz >TM1new.transcripts.fasta
    grep -c '>' TM1new.transcripts.fasta
    # 66577
    
    ### A2Du2018 - Du et al, 2018 
    zcat  /lss/research/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Du2018/1.PacBio-Gar-Assembly-v1.0/G.arboreum.Chr.v1.0.cds.v1.0.fasta.gz >A2du.transcripts.fasta
    grep -c '>' A2du.transcripts.fasta
    # 40960
    
    ### A2Du + D5 as ref for A2xD5
    cat A2du.transcripts.fasta D5.transcripts.transcripts.fa >F1.transcripts.fasta

    # build rsem-bowtie2 ref
    rsem-prepare-reference D5.transcripts.fasta --bowtie2 D5.transcripts
    rsem-prepare-reference TM1new.transcripts.fasta --bowtie2 TM1new.transcripts
    rsem-prepare-reference F1.transcripts.fasta --bowtie2 F1.transcripts
    rsem-prepare-reference A2du.transcripts.fasta --bowtie2 A2du.transcripts

In addition to mapping and expression analysis based on genome-specific transcript, I used the D5 reference transcriptomes based on SNPs 4.0 and 4.1.

    ls /lss/research/jfw-lab/GenomicResources/pseudogenomes/
    more /lss/research/jfw-lab/GenomicResources/pseudogenomes/README.txt 
    # for diploids and F1
    grep -c ">" /lss/research/jfw-lab/GenomicResources/pseudogenomes/A2D5.transcripts.fa
    # 74446
    grep -c ">" /lss/research/jfw-lab/GenomicResources/pseudogenomes/AtDt.transcripts.fa
    # 74446

     
### Bowtie2-RSEM mapping

    bash runBowtie2.D.sh >runBowtie2.D.101318.txt 2>&1
    bash runBowtie2.A.sh >runBowtie2.A.101318.txt 2>&1
    bash runBowtie2.F1.sh >runBowtie2.F1.101318.txt 2>&1
    bash runBowtie2.AD1.sh >runBowtie2.AD1.101218.txt 2>&1
    bash runBowtie2.Dref.sh >runBowtie2.Dref.101918.txt 2>&1

Mapping rate using D5 reference transcriptomes based on SNPs 4.0 and 4.1.

    grep 'overall alignment rate' /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingDref/*log
    `/work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingDref/SD5-A2-S1.log:65.09% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingDref/SD5-A2-S4.log:67.49% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingDref/SD5-A2-S5.log:67.59% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingDref/SD5-A2xD5-S1.log:71.99% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingDref/SD5-A2xD5-S2.log:71.35% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingDref/SD5-A2xD5-S3.log:66.41% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingDref/SD5-D5-S1.log:73.73% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingDref/SD5-D5-S2.log:70.11% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingDref/SD5-D5-S3.log:76.26% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingDref/SD5-Maxxa-S1.log:69.34% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingDref/SD5-Maxxa-S2.log:69.72% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingDref/SD5-Maxxa-S3.log:66.70% overall alignment rate`
   
Mapping rate using individual reference transcriptomes.

    grep 'overall alignment rate' /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingIndiv/*log
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mapping/SD5-A2-S1.log:50.49% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mapping/SD5-A2-S4.log:52.25% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mapping/SD5-A2-S5.log:53.11% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mapping/SD5-A2xD5-S1.log:74.01% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mapping/SD5-A2xD5-S2.log:73.17% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mapping/SD5-A2xD5-S3.log:67.97% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mapping/SD5-D5-S1.log:47.99% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mapping/SD5-D5-S2.log:51.41% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mapping/SD5-D5-S3.log:53.56% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mapping/SD5-Maxxa-S1.log:78.32% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mapping/SD5-Maxxa-S2.log:78.86% overall alignment rate
    /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mapping/SD5-Maxxa-S3.log:75.40% overall alignment rate


## Expression profile of gene sets

The r code [expTests.r](expTests.r) performs following:

* collect rsem count and rpkm results
* differential expression analysis
* cis-trans analysis
* impact of genome evolution: Hr, Pr, Wr
* expression dominance analyisis
