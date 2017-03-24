# Data analysis of cotton MNase-seq datasets
---

## Preprocessing of DNA sequencing datasets

### MNase-seq datasets
Short-read data will be deposited in the NCBI short read archive ([SRP??????](http://trace.ddbj.nig.ac.jp/DRASearch/study?acc=SRP??????)), also as Biobroject [PRJNA??????](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA??????).

#### Access to fastq files
Paired-end reads of 16 samples: 4 species (A2, D5, A2xD5, Maxxa) X 2 digestive conditions (Heavy and Light) X 2 technical reps.

    cd cottonLeaf/rawfastq
    ln -s ~/jfw-lab/RawData/HGJ_leafMNase-seq/WTNHHW163125/data_release/raw_data/*gz .

#### Checking read quality with [FastQC](http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/)
    cd ..
    module load fastqc/0.11.3
    fastqc rawfastq/*
    # move result files to QCreport folder
    mkdir QCreport
    mv rawfastq/*zip QCreport/
    mv rawfastq/*html QCreport/
    
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
[CottonGen](https://www.cottongen.org/data/download/genome#Ass) compiles all published cotton genomes, and I will need 4 different reference genomes for AD1, A2, D5 and A2xD5.

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
    
    # build bowtie2 ref
    bowtie2-build TM1_26.fasta TM1
    bowtie2-build F1_26.fasta F1
    bowtie2-build A2genome_13.fasta A2
    bowtie2-build Dgenome2_13.fasta D5
    ####bookmark

## Read mapping and calling of hypersensive sites
(Rodgers-Melnick et al. PNAS 2016): "After the computational trimming of adaptor sequences using CutAdapt (40), paired-end reads were mapped to the maize B73 AGPv3 reference genome, using Bowtie2 with options “no-mixed,” “no-discordant,” “no-unal,” and “dovetail” (41) for each replicate digest and for the genomic DNA. BED files were made from the resulting BAM files, using bedtools bamtobed, filtered for minimal alignment quality (≥10), and read coverage in 10-bp intervals was calculated using coverageBed (42). The DNS values were obtained by subtracting the mean normalized depth (in reads per million) of the heavy digest replicates from those of the light digest replicates. In this way, positive DNS values correspond to MNase hypersensitive footprints (as defined by ref. 8; and referred to here as MNase HS regions), whereas negative DNS values correspond to nuclease hyper-resistant footprints (MRF, as per ref. 8). A Bayes factor criterion was used to classify as significantly hypersensitive."

### Bowtie2 mapping
Default setting `-k 1` report 1 alignment for each read/pair) should work, while some might need to be modified as required by downstream tools.

    mkdir mapping
    bowtie2 -q -p 6 -t --no-mixed --no-discordant --no-unal --dovetail -x refGenomes/D5 -1 <(zcat trimmed/D1H_1_val_1.fq.gz) -2 <(zcat trimmed/D1H_2_val_2.fq.gz) -S mapping/D1H.sam 2>mapping/D1H.log
    samtools view -bS D1H.sam | samtools sort -o D1H.sort.bam ; samtools index D1H.sort.bam

* `-x maize`: use ref maize genome
* `-1 trimmed/SRR2542701_1_val_1.fq`: paired end read 1
* `-2 trimmed/SRR2542701_2_val_2.fq`: paired end read 2
* `-q`: takes fastq files
* `-p 6`: use 6 thread
* `-t`: Print the amount of wall-clock time taken by each phase.
* `--no-mixed --no-discordant`: discard discordant mapping
* `--no-unal`: Suppress SAM records for reads that failed to align.
* `--dovetail`: allow pair to overlap and over extend

### Differential nucleosome occupancy analysis - [DANPOS2](https://sites.google.com/site/danposdoc/)

### Profiling Nucleopostioning - [nucleR](http://bioconductor.org/packages/release/bioc/html/nucleR.html) & [NUCwave](http://nucleosome.usal.es/nucwave/)
Briefly, we need to examine the wave-length patterns of mapping coverage on the reference genome, and then specifically locate peaks that fit the description of nucleosomes - proper length, non-overlapping, etc.

(Zhang et al. 2015): "Well-positioned and loosely positioned nucleosomes were identified using nucleR (Flores and Orozco, 2011). ... We used the filterFFT function of nucleR to remove noise and smooth the read count score of each position along chromosomes with the parameter pcKeepComp = 0.01. After noise removal, nucleosome peaks and centers/dyads were determined using the peakDetection function (threshold = 25%, score = true, width = 140). Overlapped peaks were merged into longer regions, which were defined as loosely positioned nucleosomes, and distinct individual peaks were defined as well-positioned nucleosomes. If the length of merged peaks is longer than 150 bp, this region is considered to contain more than two nucleosome dyads and thus, contains loosely positioned nucleosomes. If the length of merged peaks is shorter than 150 bp, this region is considered to contain a well-positioned nucleosome."

The phasogram and average distance between two adjacent nucleosomes were calculated using our previously reported methods (Zhang et al., 2013). The nucleosome occupancy change scores were calculated by DANPOS (Chen et al., 2013). Analyses of dinucleotide frequency followed previously published methods (Locke et al., 2010; Valouev et al., 2011).

## Visualization and other result presentation
