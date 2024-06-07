#!/bin/bash

###########################
## 0. Set up Working Dir ##
###########################
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase
mkdir A2
mkdir D5
mkdir AD1
mkdir F1

## download
mkdir rawFastq

# 2023-2-16 system update need this line to call old modules
#module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core
#module load wget

# D5
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase/D5
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590556/Gr_DHS_rep1_1.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590556/Gr_DHS_rep1_2.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590557/Gr_DHS_rep2_1.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590557/Gr_DHS_rep2_2.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590558/Gr_DHS_rep3_1.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590558/Gr_DHS_rep3_2.fq.gz
ref='../../refGenomes/D5'
ref_chr_size='../../mappingD/chr.size.txt'
ref_fasta='../../refGenomes/Dgenome2_13.fasta'
gsize="7.5e8"
blacklist="../../ATAC/blacklist.bed"


# A2
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase/A2
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590508/Ga_DHS_rep1_1.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590508/Ga_DHS_rep1_2.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590509/Ga_DHS_rep2_1.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590509/Ga_DHS_rep2_2.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590510/Ga_DHS_rep3_1.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590510/Ga_DHS_rep3_2.fq.gz
ref='../../refGenomes/A2WHU'
ref_chr_size='chr.size.txt'
ref_fasta='../../refGenomes/A2WHU_13.fasta'

# AD1
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase/AD1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR499/006/SRR4996276/SRR4996276_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR499/006/SRR4996276/SRR4996276_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR499/009/SRR4996279/SRR4996279_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR499/009/SRR4996279/SRR4996279_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR499/001/SRR4996281/SRR4996281_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR499/001/SRR4996281/SRR4996281_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR499/002/SRR4996282/SRR4996282_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR499/002/SRR4996282/SRR4996282_2.fastq.gz
ref='../../refGenomes/AD1'
ref_chr_size='../../mappingM_UTX/chr.size.txt'
ref_fasta='../../refGenomes/TM1utx_26.fasta'

# F1
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase/F1
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590518/F1_DHS_rep1_1.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590518/F1_DHS_rep1_2.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590519/F1_DHS_rep2_1.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590519/F1_DHS_rep2_2.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590520/F1_DHS_rep3_1.fq.gz
wget ftp.sra.ebi.ac.uk/vol1/run/ERR659/ERR6590520/F1_DHS_rep3_2.fq.gz
# AD1 ref
ref='../../refGenomes/AD1'
ref_chr_size='../../mappingM_UTX/chr.size.txt'
ref_fasta='../../refGenomes/TM1utx_26.fasta'
# A2+D5 ref
ref='../../refGenomes/F2020'
ref_chr_size='chr.size.txt'
ref_fasta='../../refGenomes/F2020_26.fasta'

##########################
## 1. Download and trim ##
##########################

module load trimgalore
SUFFIX=_1.fastq.gz
for j in $( ls | grep _1.fq.gz$ );    do
echo ''
echo "==========Running Trim_galore for"
echo ${j%$SUFFIX}
time trim_galore --cores 4 --paired -o . $j ${j%$SUFFIX}_2.fq.gz
done

### Quality Control
module load fastqc
module load py-multiqc
fastqc *fq.gz
multiqc .

#################################################
## 2. General workflow of DNase-seq with MACS2 ##
#################################################
module load r/3.5.0-py2-ufvuwmm
module load bowtie2 samtools
module load picard

## From trimmed FASTQ mapping to refence genome, filtering and quality control: Bowtie2, Samtools, phantompeakqualtools

echo ''
echo "======Map sequences to the genome"
for j in $( ls | grep _1_val_1.fq.gz$ );    do
echo ''
echo ${j%%_*}
bowtie2 -p $ncores -t --no-mixed --no-discordant --no-unal --dovetail -x $ref -1 $j -2 ${j%%_*}_2_val_2.fq.gz > ${j%%_*}.sam 2>${j%%_*}.log
samtools view -Sb ${j%%_*}.sam | samtools sort - -o ${j%%_*}.bam ; samtools index ${j%%_*}.bam
picard MarkDuplicates I=${j%%_*}.bam O=${j%%_*}.nodup.bam M=${j%%_*}.dups.txt REMOVE_DUPLICATES=true
samtools view -q 20 -b ${j%%_*}.nodup.bam > ${j%%_*}.nodup.q20.bam; samtools index ${j%%_*}.nodup.q20.bam
Rscript /work/LAS/jfw-lab/hugj2006/tools/phantompeakqualtools/run_spp.R -c=${j%%_*}.nodup.q20.bam -savp -out=${j%%_*}.nodup.q20.phantomPeak.txt >${j%%_*}.nodup.q20.phantomPeak.log
done

echo ''
echo "======Compare replicate BAMs"
module purge
module load py-deeptools
multiBamSummary bins --bamfiles *nodup.q20.bam -o MBresults.npz
plotCorrelation -in MBresults.npz --corMethod spearman --skipZeros --plotTitle 'Spearman Correlation of Read Counts' --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o heatmap_SpearmanCorr_readCounts.png --outFileCorMatrix SpearmanCorr_readCounts.tab
plotCorrelation -in MBresults.npz --corMethod pearson --skipZeros --plotTitle 'Pearson Correlation of Read Counts' --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o heatmap_PearsonCorr_readCounts.png --outFileCorMatrix PearsonnCorr_readCounts.tab
bamPEFragmentSize -hist fragmentSize.png -T 'Fragment size of PE DNase-seq data' --maxFragmentLength 1000 -b *nodup.q20.bam
# Generate BigWig, IGV visualization to identify blacklist region, exclude blastlist region
module load bedtools2 samtools
for j in $( ls *nodup.q20.bam );    do
echo ${j%%[.]*}
bamCoverage -b $j -o ${j%%[.]*}.bw
bedtools intersect -v -abam $j -b $blacklist >  ${j%%[.]*}.nodup.q20.filter.bam; samtools index ${j%%[.]*}.nodup.q20.filter.bam
bamCoverage -b ${j%%[.]*}.nodup.q20.filter.bam -o ${j%%[.]*}.filter.bw
done

echo ''
echo "======MACS2 peak calling"
module purge
module load py-macs2
module load bedtools2
for j in $( ls *nodup.q20.filter.bam );    do
macs2 callpeak -t $j -f BAMPE -g $gsize -n ${j%%[.]*}.macs2 --keep-dup auto --outdir .
#bedtools intersect -v -a ${j%%[.]*}.macs2_peaks.narrowPeak -b $blacklist >${j%%[.]*}.macs2_peaks.narrowPeak.filtered.bed
done
#intersect and combine
bedtools intersect -a $(ls *macs2_peaks.narrowPeak | head -1) -b $(ls *macs2_peaks.narrowPeak | tail -1) > macs2.intersect.bed
cat *macs2_peaks.narrowPeak | bedtools sort | bedtools merge > macs2.combine.bed

echo ''
echo "======Genrich peak calling"
###------ Genrich peak calling
module load samtools
for j in $( ls *sam );    do
echo ''
echo $j
samtools sort -n $j -o ${j%.sam}.sortN.bam
done
/work/LAS/jfw-lab/hugj2006/tools/Genrich/Genrich -t SRR5886144.sortN.bam,SRR5886145.sortN.bam \
  -o genrich.narrowPeak  \
  -f genrich.PQvalue.txt -k genrich.pileup.txt \
  -r  -x  -q 0.05  -a 20.0  -v  \
  -E $blacklist
  >genrich.log 2>&1
# split pileup files to bedgrapgh each rep
#sed s/chr/#chr/g genrich.pileup.txt | cut -f1-4 | csplit -z - /'# experimental file'/ {*}
# bedgrapgh to bigwig
#/work/LAS/jfw-lab/hugj2006/tools/kent/bedGraphToBigWig xx00 $ref_chr_size genrich.pileup.SRR5886142.bw
#/work/LAS/jfw-lab/hugj2006/tools/kent/bedGraphToBigWig xx01 $ref_chr_size genrich.pileup.SRR5886143.bw
bedtools intersect -a macs2.intersect.bed -b genrich.narrowPeak > MvG.intersect.bed
bedtools intersect -a macs2.combine.bed -b genrich.narrowPeak > MvG.combine.bed
# Genrich not suitable for DNase-seq: tested for D genome, macs2 called >25k peaks, genrich called ~2k peaks with ~400 overlapped with macs2 result

echo ''
echo "======Generate sequence and mapping statistics"
module load r/3.5.0-py2-ufvuwmm
Rscript dnase.summary.r

```{dnase.summary.r}
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase/")
fqL <- c(grep("_1.fq.gz$",list.files("A2", full=T),value=T),grep("_1.fq.gz$",list.files("D5", full=T),value=T),grep("_1.fq.gz$",list.files("AD1", full=T),value=T))
###------Generate sequence and mapping statistics:
sumStat<-data.frame(sample=fqL, Type="DNase-seq", raw.reads=NA, filtered.reads=NA, mapping.input=NA, mapped=NA, mapped.Perc=NA, mapped.unique=NA, mapped.unique.Perc=NA, corr_estFragLen=NA, NSC=NA, RSC=NA,Qtag=NA)
rownames(sumStat) = gsub("_.*","",fqL)
for(fq in fqL){
    each = gsub("_.*","",fq)
    # from trimming report
    p1 <- readLines(paste0(each,"_1.fastq.gz_trimming_report.txt"))
    p1.raw.reads = gsub(".* ","",grep("Total reads",p1,value=TRUE))
    p1.filtered.reads = gsub(".*  | \\(.*","",grep("Reads written",p1,value=TRUE))
    p2 <- readLines(paste0(each,"_2.fastq.gz_trimming_report.txt"))
    p2.raw.reads = gsub(".* ","",grep("Total reads",p2,value=TRUE))
    p2.filtered.reads= gsub(".*  | \\(.*","",grep("Reads written",p2,value=TRUE))
    sumStat[each,"raw.reads"]=ifelse(p1.raw.reads==p2.raw.reads,p1.raw.reads,"Error: p1!=p2")
    sumStat[each,"filtered.reads"]=ifelse(p1.filtered.reads==p2.filtered.reads,p1.filtered.reads,"error: p1!=p2")
    # from mapping
    m <- readLines(paste0(each,".log"))
    input.reads <- as.numeric(gsub(" reads; of these:","",grep("reads",m,value=TRUE)))
    unmapped <- as.numeric(gsub(" |\\(.*","",grep("0 times",m,value=TRUE)))
    mapped <-input.reads-unmapped
    mappedUni <- as.numeric(gsub(" |\\(.*","",grep("exactly 1 time",m,value=TRUE)))
    sumStat[each,"mapping.input"]=input.reads
    sumStat[each,"mapped"]=mapped
    sumStat[each,"mapped.unique"]=mappedUni
    # from PCARD output
    p <- readLines(paste0(each,".dups.txt"))[8]
    ps<-strsplit(p,split="\t")
    sumStat[each,"duplicates"] =ps[[1]][7]
    sumStat[each,"optical.duplicates"] =ps[[1]][8]
    sumStat[each,"percent.duplication"] =ps[[1]][9]
    sumStat[each,"after.dup.removal"] =ps[[1]][10]
    # from phantompeakqual
    p <- read.table(paste0(each,".nodup.q20.phantomPeak.txt"),sep="\t")
    sumStat[each,"corr_estFragLen"] = p$V3
    sumStat[each,"NSC"] = p$V9
    sumStat[each,"RSC"] =p$V10
    sumStat[each,"Qtag"] =p$V11
}
sumStat$mapped.Perc=round(sumStat$mapped/sumStat$mapping.input*100,2)
sumStat$mapped.unique.Perc=round(sumStat$mapped.unique/sumStat$mapping.input*100,2)
sumStat
write.table(sumStat, "seqSumStats.txt", quote=FALSE, sep="\t")

###-------Report peak calling summary
library(data.table)
bedL <- list.files("D5",pattern="narrowPeak$|intersect.bed$|combine.bed$", full=T)
peakStat = data.frame(Bed=bedL)
rownames(peakStat)=bedL
for(each in bedL){
    df<-fread(each)
    peakStat[each,"peakN"] = nrow(df)
    peakStat[each,"peakBp"] = sum(df$V3-df$V2)
}
peakStat$bpPerGenome = peakStat$peakBp/749228090
peakStat
t=peakStat
#
bedL <- list.files("A2",pattern="narrowPeak$|intersect.bed$|combine.bed$", full=T)
peakStat = data.frame(Bed=bedL)
rownames(peakStat)=bedL
for(each in bedL){
    df<-fread(each)
    peakStat[each,"peakN"] = nrow(df)
    peakStat[each,"peakBp"] = sum(df$V3-df$V2)
}
peakStat$bpPerGenome = peakStat$peakBp/1509187826
peakStat
t= rbind(t,peakStat)
#
bedL <- list.files("AD1",pattern="narrowPeak$|intersect.bed$|combine.bed$", full=T)
peakStat = data.frame(Bed=bedL)
rownames(peakStat)=bedL
for(each in bedL){
    df<-fread(each)
    peakStat[each,"peakN"] = nrow(df)
    peakStat[each,"peakBp"] = sum(df$V3-df$V2)
}
peakStat$bpPerGenome = peakStat$peakBp/2281618630
t= rbind(t,peakStat)
write.table(t, "peakSumStats.txt", quote=FALSE, sep="\t",row.names=F)



```




