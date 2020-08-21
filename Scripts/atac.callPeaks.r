###------Generate working dir structure
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/ATAC")
system("bash makeWorkenvir.sh")

###------Locate FASTQ file
system("ln -s /work/LAS/jfw-lab/ATAC/G.raimondii/ATAC/*gz rawFastq/")

###------Quality control of fastq: trim, filter and reports: Trim_galore, FastQC
cat("module load python\nmodule load trimgalore\nmodule load py-setuptools\nmodule load py-multiqc", file="runQC.sh",sep="\n")
fqL <- list.files("rawFastq",pattern="_R1",full=T)
for(each in fqL){
    cmd <- paste0("trim_galore --paired -o trimmedFastq/ ",each," ",gsub("_R1","_R2",each))
    cat(cmd, file="runQC.sh",sep="\n", append=TRUE)
}
# QC raw and trimmed fastq
cat("module load fastqc", file="runQC.sh",sep="\n", append=TRUE)
cat("fastqc -o QCreport/raw/ rawFastq/*",file="runQC.sh",sep="\n", append=TRUE)
cat("fastqc -o QCreport/trimmed/ trimmedFastq/*",file="runQC.sh",sep="\n", append=TRUE)
cat("cd QCreport/raw/\nmultiqc .\ncd ../trimmed/\nmultiqc .",file="runQC.sh",sep="\n", append=TRUE)
system("cat runQC.sh")
system("bash runQC.sh")
# inspect multiqc results

###------Read mapping to refence genome, filtering and quality control: Bowtie2, Samtools, phantompeakqualtools
#sh# bowtie2 -q -p 8 -t --no-mixed --no-discordant --no-unal --dovetail -x ../refGenomes/TM1new -1 <(zcat trimmedFastq/SRR4996230_1_val_1.fq.gz) -2 <(zcat trimmedFastq/SRR4996230_2_val_2.fq.gz) -S mapping/SRR4996230.sam 2>mapping/SRR4996230.log
#sh# samtools view -Sb mapping/SRR4996230.sam | samtools sort - -o mapping/SRR4996230.bam ; samtools index mapping/SRR4996230.bam
#sh# samtools view -q 20 -b mapping/SRR4996230.bam > mapping/SRR4996230.q20.bam
#sh# samtools view mapping/SRR4996230.q20.bam | grep 'XS:' |head
#sh# samtools view mapping/SRR4996230.q20.bam | grep 'XT:A:U' | samtools view -bS -T referenceSequence.fa - > reads.uniqueMap.bam
# Calculate fragment length, NSC and RSC by phantompeakqualtools
#sh# Rscript /work/LAS/jfw-lab/hugj2006/tools/phantompeakqualtools/run_spp.R -c=mapping/SRR4996230.q20.bam -savp -out=mapping/SRR4996230.q20.phantomPeak.txt >mapping/SRR4996230.q20.phantomPeak.log
# sequence alignment with Bowtie2, set parameters first
ref<-"../refGenomes/D5"
nThread <- 8
# prepare bash
cat("module load bowtie2 samtools bedtools2\nmodule load picard\n", file="runBowtie2.sh",sep="\n")
fqL <- gsub("_R1_val_1.fq.gz","",list.files("trimmedFastq",pattern="_R1_val_1"))
for(each in fqL){
    cmd.bowtie2 <- paste0("bowtie2 -q -p ",nThread," -t --no-mixed --no-discordant --no-unal --dovetail -x ", ref," -1 <(zcat trimmedFastq/",each,"_R1_val_1.fq.gz) -2 <(zcat trimmedFastq/", each,"_R2_val_2.fq.gz) -S mapping/", each,".sam 2>mapping/", each,".log")
cmd.sam2bam <- paste0("samtools view -Sb mapping/", each, ".sam | samtools sort - -o mapping/", each, ".bam ; samtools index mapping/", each,".bam")
cmd.picard <- paste0("picard MarkDuplicates I=mapping/",each,".bam O=mapping/",each,".nodup.bam M=mapping/",each,".dups.txt REMOVE_DUPLICATES=true")
cmd.filter <- paste0("samtools view -q 20 -b mapping/",each,".nodup.bam | bedtools intersect -v -abam - -b blacklist.bed >  mapping/",each,".nodup.q20.bam; samtools index mapping/",each,".nodup.q20.bam")
cmd.qc <- paste0("Rscript /work/LAS/jfw-lab/hugj2006/tools/phantompeakqualtools/run_spp.R -c=mapping/",each,".nodup.q20.bam -savp -out=mapping/",each,".nodup.q20.phantomPeak.txt >mapping/",each,".nodup.q20.phantomPeak.log\n")
cat(cmd.bowtie2,cmd.sam2bam,cmd.picard, cmd.filter,cmd.qc, file="runBowtie2.sh",sep="\n", append=TRUE)
}
system("cat runBowtie2.sh")
system("bash runBowtie2.sh >runBowtie2.xxxx18.log")
# how to calculate library complexity matric? ENCODE prefered values are NRF>0.9, PBC1>0.9, and PBC2>10.


###------Generate sequence and mapping statistics:
sumStat<-data.frame(sample=fqL, Type="ATAC-seq", raw.reads=NA, filtered.reads=NA, mapping.input=NA, mapped=NA, mapped.Perc=NA, mapped.unique=NA, mapped.unique.Perc=NA, corr_estFragLen=NA, NSC=NA, RSC=NA,Qtag=NA)
rownames(sumStat) = fqL
for(each in fqL){
    # from trimming report
    p1 <- readLines(paste0("trimmedFastq/",each,"_R1.fastq.gz_trimming_report.txt"))
    p1.raw.reads = gsub(".* ","",grep("Total reads",p1,value=TRUE))
    p1.filtered.reads = gsub(".*  | \\(.*","",grep("Reads written",p1,value=TRUE))
    p2 <- readLines(paste0("trimmedFastq/",each,"_R2.fastq.gz_trimming_report.txt"))
    p2.raw.reads = gsub(".* ","",grep("Total reads",p2,value=TRUE))
    p2.filtered.reads= gsub(".*  | \\(.*","",grep("Reads written",p2,value=TRUE))
    sumStat[each,"raw.reads"]=ifelse(p1.raw.reads==p2.raw.reads,p1.raw.reads,"Error: p1!=p2")
    sumStat[each,"filtered.reads"]=ifelse(p1.filtered.reads==p2.filtered.reads,p1.filtered.reads,"error: p1!=p2")
    # from mapping
    m <- readLines(paste0("mapping/",each,".log"))
    input.reads <- as.numeric(gsub(" reads; of these:","",grep("reads",m,value=TRUE)))
    unmapped <- as.numeric(gsub(" |\\(.*","",grep("0 times",m,value=TRUE)))
    mapped <-input.reads-unmapped
    mappedUni <- as.numeric(gsub(" |\\(.*","",grep("exactly 1 time",m,value=TRUE)))
    sumStat[each,"mapping.input"]=input.reads
    sumStat[each,"mapped"]=mapped
    sumStat[each,"mapped.unique"]=mappedUni
    # from PCARD output
    p <- readLines(paste0("mapping/",each,".dups.txt"))[8]
    ps<-strsplit(p,split="\t")
    sumStat[each,"duplicates"] =ps[[1]][7]
    sumStat[each,"optical.duplicates"] =ps[[1]][8]
    sumStat[each,"percent.duplication"] =ps[[1]][9]
    sumStat[each,"after.dup.removal"] =ps[[1]][10]
    # from phantompeakqual
    p <- read.table(paste0("mapping/",each,".nodup.q20.phantomPeak.txt"),sep="\t")
    sumStat[each,"corr_estFragLen"] = p$V3
    sumStat[each,"NSC"] = p$V9
    sumStat[each,"RSC"] =p$V10
    sumStat[each,"Qtag"] =p$V11
}
sumStat$mapped.Perc=round(sumStat$mapped/sumStat$mapping.input*100,2)
sumStat$mapped.unique.Perc=round(sumStat$mapped.unique/sumStat$mapping.input*100,2)
sumStat
write.table(sumStat, "seqSumStats.txt", quote=FALSE, sep="\t")

###-------Compare replicates
cat("## Compare replicates\nmodule load py-deeptools", file="compareReps.sh",sep="\n")
cmd <- paste0("multiBamSummary bins --bamfiles ",paste(list.files("mapping",pattern="q20.bam$",full=T), collapse=" "), " -o mapping/MBresults.npz; plotCorrelation -in mapping/MBresults.npz --corMethod spearman --skipZeros --plotTitle 'Spearman Correlation of Read Counts' --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o mapping/heatmap_SpearmanCorr_readCounts.png --outFileCorMatrix mapping/SpearmanCorr_readCounts.tab; plotCorrelation -in mapping/MBresults.npz --corMethod pearson --skipZeros --plotTitle 'Pearson Correlation of Read Counts' --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o mapping/heatmap_PearsonCorr_readCounts.png --outFileCorMatrix mapping/PearsonnCorr_readCounts.tab")
cat(cmd, file="compareReps.sh",sep="\n", append=TRUE)
cmd<-paste0("bamPEFragmentSize -hist mapping/fragmentSize.png -T 'Fragment size of PE ATAC-seq data' --maxFragmentLength 1000 -b ",paste(list.files("mapping",pattern="q20.bam$",full=T), collapse=" ")," --samplesLabel rep1 rep2")
cat(cmd, file="compareReps.sh",sep="\n", append=TRUE)
system("cat compareReps.sh")
system("bash compareReps.sh")


###-------HOMER peak calling
gsize="7.5e8"
# http://plant-plasticity.github.io/resources/3_ATAC-seq%20data%20processing.pdf
cat("## call peaks with HOMER\nexport PATH=$PATH:/work/LAS/jfw-lab/hugj2006/tools/homer/bin/\nmodule load bedtools2\nmkdir peaks", file="runHomer.sh",sep="\n")
bamL <- list.files("mapping",pattern="*q20.bam$", full=T)
for(each in bamL){
    flag=gsub("mapping/|.nodup.q20.bam","",each)
    cmd.tag <- paste0("makeTagDirectory peaks/homerTag_",flag, " ",each)
    cmd.call <- paste0("findPeaks peaks/homerTag_",flag, " -o peaks/",flag,".peaks.txt -gsize ",gsize," minDist 150 -region")
    cmd.2bed <- paste0("pos2bed.pl peaks/",flag,".peaks.txt| bedtools sort | bedtools merge >peaks/",flag,".merged.bed")
    # filter out the chr1 mt region
    #cmd.filter <-paste0("bedtools intersect -v -a peaks/",flag,".merged.bed -b blacklist.bed >peaks/",flag,".merged.filtered.bed")
    cat(cmd.tag,cmd.call,cmd.2bed, file="runHomer.sh",sep="\n", append=TRUE)
}
bedL <- paste0("peaks/",gsub("mapping/|.nodup.q20.bam","",bamL),".merged.bed")
cmd1 = paste0("bedtools intersect -a ",bedL[1]," -b ", bedL[2]," > peaks/homer.intersect.bed")
cmd2 = paste0("cat ",bedL[1]," ", bedL[2],"| bedtools sort | bedtools merge > peaks/homer.combine.bed")
cat(cmd1,cmd2, file="runHomer.sh",sep="\n", append=TRUE)
system("cat runHomer.sh")
system("bash runHomer.sh")

###-------MACS2 peak calling
cat("## call peaks with MACS2\nmodule load py-macs2\nmodule load bedtools2\n", file="runMacs2.sh",sep="\n")
bamL <- list.files("mapping",pattern="*q20.bam$", full=T)
for(each in bamL){
    flag=gsub("mapping/|.nodup.q20.bam","",each)
    cmd.call <- paste0("macs2 callpeak -t ",each ," -f BAMPE -g ",gsize," -n ",flag,".macs2 --keep-dup auto --outdir peaks/")
    #cmd.filter <-paste0("bedtools intersect -v -a peaks/",flag,".macs2_peaks.narrowPeak -b blacklist.bed >peaks/",flag,".macs2_peaks.narrowPeak.filtered.bed")
    cat(cmd.call, file="runMacs2.sh",sep="\n", append=TRUE)
}
npL <- paste0("peaks/",gsub("mapping/|.nodup.q20.bam","",bamL),".macs2_peaks.narrowPeak")
cmd1 = paste0("bedtools intersect -a ",npL[1]," -b ", npL[2]," > peaks/macs2.intersect.bed")
cmd2 = paste0("cat ",npL[1]," ", npL[2],"| bedtools sort | bedtools merge > peaks/macs2.combine.bed")
cmd3=paste0("for j in $( ls mapping/*nodup.q20.bam; do bamCoverage -b $j -o ${j%%[.]*}.filter.bw; done")
cat(cmd1, cmd2, cmd3, file="runMacs2.sh",sep="\n", append=TRUE)
system("cat runMacs2.sh")
system("bash runMacs2.sh")

###-------MACS2 vs Homer peak calling
system("module load bedtools2; bedtools intersect -a peaks/macs2.intersect.bed -b peaks/homer.intersect.bed > peaks/MvH.eachIntersect.intersect.bed; bedtools intersect -a peaks/macs2.combine.bed -b peaks/homer.combine.bed > peaks/MvH.intersect.bed")

###------ Genrich peak calling
system("bash atac.genrich_command.sh")
system("module load bedtools2; bedtools intersect -a peaks/macs2.intersect.bed -b peaks/genrich.narrowPeak > peaks/MvG.intersect.bed; bedtools intersect -a peaks/homer.combine.bed -b peaks/genrich.narrowPeak > peaks/HvG.intersect.bed")


###-------Report peak calling summary
library(data.table)
bedL <- list.files("peaks",pattern="narrowPeak$|merged.bed$|intersect.bed$|combine.bed$", full=T)
peakStat = data.frame(Bed=bedL)
rownames(peakStat)=bedL
for(each in bedL){
    df<-fread(each)
    peakStat[each,"peakN"] = nrow(df)
    peakStat[each,"peakBp"] = sum(df$V3-df$V2)
}
peakStat$bpPerGenome = peakStat$peakBp/749228090
peakStat$method="intersection"
peakStat$method[grep("macs",peakStat$Bed)]="macs"
peakStat$method[grep("homer|merged",peakStat$Bed)]="homer"
peakStat$method[grep("genrich",peakStat$Bed)]="genrich"
peakStat$method[grep("macs",peakStat$Bed)]="macs"
peakStat$rep="both"
peakStat$rep[grep("R25",peakStat$Bed)]="R25"
peakStat$rep[grep("R26",peakStat$Bed)]="R26"
peakStat
write.table(peakStat, "peakSumStats.txt", quote=FALSE, sep="\t",row.names=F)
