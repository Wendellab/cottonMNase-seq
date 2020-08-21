#!/bin/bash

###########################
## 0. Set up Working Dir ##
###########################
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase
mkdir RNA_A2
mkdir RNA_D5
mkdir RNA_AD1
mkdir MappingIndiv
mkdir MappingDref

module load trimgalore kallisto

# D5
echo ''
echo "==Running D5"
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase/RNA_D5
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/008/SRR5886148/SRR5886148_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/008/SRR5886148/SRR5886148_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/009/SRR5886149/SRR5886149_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/009/SRR5886149/SRR5886149_2.fastq.gz
ref='/work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/refTranscriptomes/D5.transcripts.kallisto'
ref_D='/work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/refTranscriptomes/A2D5.transcripts.kallisto'
SUFFIX=_1.fastq.gz
for j in $( ls | grep _1.fastq.gz$ );    do
echo ''
echo "==========Running Trim_galore for"
echo ${j%$SUFFIX}
trim_galore --cores 4 --paired -o . $j ${j%$SUFFIX}_2.fastq.gz
echo "==========Kallisto"
kallisto quant -t $ncores -i $ref -o ../MappingIndiv/${j%%_*} ${j%%_*}_1_val_1.fq.gz ${j%%_*}_2_val_2.fq.gz 2>../MappingIndiv/${j%%_*}.log
kallisto quant -t $ncores -i $ref_D -o ../MappingDref/${j%%_*} ${j%%_*}_1_val_1.fq.gz ${j%%_*}_2_val_2.fq.gz 2>../MappingDref/${j%%_*}.log
done
echo "==========Quality control"
module load fastqc
module load py-multiqc
fastqc *fq.gz
multiqc .

# A2
echo ''
echo "==Running A2"
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase/RNA_A2
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/005/SRR5886155/SRR5886155_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/005/SRR5886155/SRR5886155_2.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/006/SRR5886156/SRR5886156_1.fastq.gz
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/006/SRR5886156/SRR5886156_2.fastq.gz
ref='/work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/refTranscriptomes/A2_WHU.transcripts.kallisto'
ref_D='/work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/refTranscriptomes/A2D5.transcripts.kallisto'
SUFFIX=_1.fastq.gz
for j in $( ls | grep _1.fastq.gz$ );    do
echo ''
echo "==========Running Trim_galore for"
echo ${j%$SUFFIX}
trim_galore --cores 4 --paired -o . $j ${j%$SUFFIX}_2.fastq.gz
echo "==========Kallisto"
kallisto quant -t $ncores -i $ref -o ../MappingIndiv/${j%%_*} ${j%%_*}_1_val_1.fq.gz ${j%%_*}_2_val_2.fq.gz 2>../MappingIndiv/${j%%_*}.log
kallisto quant -t $ncores -i $ref_D -o ../MappingDref/${j%%_*} ${j%%_*}_1_val_1.fq.gz ${j%%_*}_2_val_2.fq.gz 2>../MappingDref/${j%%_*}.log
done
echo "==========Quality control"
module load fastqc
module load py-multiqc
fastqc *fq.gz
multiqc .


# AD1
echo ''
echo "==Running AD1"
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase/RNA_AD1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/007/SRR5886147/SRR5886147_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/007/SRR5886147/SRR5886147_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/006/SRR5886146/SRR5886146_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/006/SRR5886146/SRR5886146_2.fastq.gz
ref='/work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/refTranscriptomes/AD1_UTX.transcripts.kallisto'
ref_D='/work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/refTranscriptomes/AtDt.transcripts.kallisto'
SUFFIX=_1.fastq.gz
for j in $( ls | grep _1.fastq.gz$ );    do
echo ''
echo "==========Running Trim_galore for"
echo ${j%$SUFFIX}
trim_galore --cores 4 --paired -o . $j ${j%$SUFFIX}_2.fastq.gz
echo "==========Kallisto"
kallisto quant -t $ncores -i $ref -o ../MappingIndiv/${j%%_*} ${j%%_*}_1_val_1.fq.gz ${j%%_*}_2_val_2.fq.gz 2>../MappingIndiv/${j%%_*}.log
kallisto quant -t $ncores -i $ref_D -o ../MappingDref/${j%%_*} ${j%%_*}_1_val_1.fq.gz ${j%%_*}_2_val_2.fq.gz 2>../MappingDref/${j%%_*}.log
done
echo "==========Quality control"
module load fastqc
module load py-multiqc
fastqc *fq.gz
multiqc .

# Summary
echo ''
echo "==Make RNA-seq summary table"
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase/
grep 'Total reads processed' RNA*/*trimming_report.txt >RNA.summary.text
grep 'Reads written' RNA*/*trimming_report.txt >>RNA.summary.text
grep "pseudoaligned" MappingIndiv/*log |cut -f3,5 >>RNA.summary.text
grep "pseudoaligned" MappingDref/*log |cut -f3,5 >>RNA.summary.text

# DE analysis, etc.
Rscript expTests.YL.r
