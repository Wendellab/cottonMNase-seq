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
cat("module load bowtie2 samtools\n", file="runBowtie2.sh",sep="\n")
fqL <- gsub("_R1_val_1.fq.gz","",list.files("trimmedFastq",pattern="_R1_val_1"))
for(each in fqL){
    cmd.bowtie2 <- paste0("bowtie2 -q -p ",nThread," -t --no-mixed --no-discordant --no-unal --dovetail -x ", ref," -1 <(zcat trimmedFastq/",each,"_R1_val_1.fq.gz) -2 <(zcat trimmedFastq/", each,"_R2_val_2.fq.gz) -S mapping/", each,".sam 2>mapping/", each,".log")
cmd.sam2bam <- paste0("samtools view -Sb mapping/", each, ".sam | samtools sort - -o mapping/", each, ".bam ; samtools index mapping/", each,".bam")
cmd.filter <- paste0("samtools view -q 20 -b mapping/",each,".bam > mapping/",each,".q20.bam; samtools index mapping/",each,".q20.bam")
cmd.qc <- paste0("Rscript /work/LAS/jfw-lab/hugj2006/tools/phantompeakqualtools/run_spp.R -c=mapping/",each,".q20.bam -savp -out=mapping/",each,".q20.phantomPeak.txt >mapping/",each,".q20.phantomPeak.log\n")
cat(cmd.bowtie2,cmd.sam2bam,cmd.filter,cmd.qc, file="runBowtie2.sh",sep="\n", append=TRUE)
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
    # from phantompeakqual
    p <- read.table(paste0("mapping/",each,".q20.phantomPeak.txt"),sep="\t")
    sumStat[each,"corr_estFragLen"] = p$V3
    sumStat[each,"NSC"] = p$V9
    sumStat[each,"RSC"] =p$V10
    sumStat[each,"Qtag"] =p$V11
}
sumStat$mapped.Perc=round(sumStat$mapped/sumStat$mapping.input*100,2)
sumStat$mapped.unique.Perc=round(sumStat$mapped.unique/sumStat$mapping.input*100,2)
sumStat
write.table(sumStat, "seqSumStats.txt", quote=FALSE, sep="\t")

###-------HOMER peak calling
# http://plant-plasticity.github.io/resources/3_ATAC-seq%20data%20processing.pdf
cat("## call peaks with HOMER\nexport PATH=$PATH:/work/LAS/jfw-lab/hugj2006/tools/homer/bin/\nmodule load bedtools2\nmkdir peaks", file="runHomer.sh",sep="\n")
bamL <- list.files("mapping",pattern="*q20.bam$", full=T)
for(each in bamL){
    flag=gsub("mapping/|.q20.bam","",each)
    cmd.tag <- paste0("makeTagDirectory peaks/homerTag_",flag, " ",each)
    cmd.call <- paste0("findPeaks peaks/homerTag_",flag, " -o peaks/",flag,".peaks.txt -gsize 7.1e8 minDist 150 -region")
    cmd.2bed <- paste0("pos2bed.pl peaks/",flag,".peaks.txt| bedtools sort | bedtools merge >peaks/",flag,".merged.bed")
    # filter out the chr1 mt region
    cmd.filter <-paste0("bedtools intersect -v -a peaks/",flag,".merged.bed -b blacklist.bed >peaks/",flag,".merged.filtered.bed")
    cat(cmd.tag,cmd.call,cmd.2bed,cmd.filter, file="runHomer.sh",sep="\n", append=TRUE)
}
system("cat runHomer.sh")
system("bash runHomer.sh")
system("bedtools intersect -a peaks/raimondii_R25_ATAC.merged.bed -b peaks/raimondii_R26_ATAC.merged.bed > peaks/homer.res.bed")
library(data.table)
df<-fread("peaks/raimondii_R25_ATAC.merged.bed"); sum(df$V3-df$V2) #6886316 bp
df<-fread("peaks/raimondii_R26_ATAC.merged.bed"); sum(df$V3-df$V2) #3560336
df<-fread("peaks/homer.res.bed"); sum(df$V3-df$V2) #2987522


###-------MACS2 peak calling
cat("## call peaks with MACS2\nmodule load py-macs2\nmodule load bedtools2\n", file="runMacs2.sh",sep="\n")
bamL <- list.files("mapping",pattern="*q20.bam$", full=T)
for(each in bamL){
    flag=gsub("mapping/|.q20.bam","",each)
    cmd.call <- paste0("macs2 callpeak -t ",each ," -f BAMPE -g 7.1e8 -n ",flag,".macs2 --keep-dup auto --outdir peaks/")
    cmd.filter <-paste0("bedtools intersect -v -a peaks/",flag,".macs2_peaks.narrowPeak -b blacklist.bed >peaks/",flag,".macs2_peaks.narrowPeak.filtered.bed")
    cat(cmd.call,cmd.filter, file="runMacs2.sh",sep="\n", append=TRUE)
}
system("cat runMacs2.sh")
system("bash runMacs2.sh")
system("bedtools intersect -a peaks/raimondii_R25_ATAC.macs2_peaks.narrowPeak.filtered.bed -b peaks/raimondii_R26_ATAC.macs2_peaks.narrowPeak.filtered.bed > peaks/macs2.res.bed")
library(data.table)
df<-fread("peaks/raimondii_R25_ATAC.macs2_peaks.narrowPeak.filtered.bed"); sum(df$V3-df$V2) #5699116
df<-fread("peaks/raimondii_R26_ATAC.macs2_peaks.narrowPeak.filtered.bed"); sum(df$V3-df$V2) #4229540
df<-fread("peaks/macs2.res.bed"); sum(df$V3-df$V2) #2827122


###------iSeg peak calling
cat("## call peaks with iSeg v1.3.2 in bigram\nmodule load samtools\nmodule load bedtools2\nmodule load py-deeptools\n", file="runiSeg.pre.sh",sep="\n")
bams <- list.files("mapping",pattern="*q20.bam$", full=T)
# bam to bed
for(i in bams){
    cmd <- paste0("samtools sort -n ",i," | bedtools bamtobed -i - -bedpe  | cut -f 1,2,6,7,8,9 | sort -T . -k1,1 -k2,2n -S 1G  > ",gsub("bam","bed",i))
    cat(cmd, file="runiSeg.pre.sh",sep="\n", append=TRUE)
}
system("cat runiSeg.sh")
system("bash runiSeg.sh")
# bed -> bedgraph ->qnorm
cat("#use library Travis in R for bed -> bedgraph ->qnorm\n", file="runiSeg.pre.sh",sep="\n",append=TRUE)
library(travis)
beds = files("mapping/*q20.bed")
w= "../mappingD/chr.size_w20.bed"
bgs=bedtoolsCoverage(beds,w)
bgs  # bgs=list.files(pattern="w20.bg")
qbgs=bgQuantileNorm(bgs)
qbgs # qbgs=list.files(pattern="qnorm.bg")
# chromosome length and non-zero region
sl=list()
nonzero = list()
for(bg in bgs)
{
    flag = gsub("raimondii_|_ATAC.*","",bg)
    dt=fread(bg,sep="\t")
    sl[[flag]] = tapply(dt$V3,dt$V1,max)
    dt$width = dt$V3-dt$V2
    select= (dt$V4!=0 )
    x = aggregate(dt$width[select],by=list(dt$V1[select]),sum)
    nonzero[[flag]] =  x$x
    names(nonzero[[flag]]) =x$Group.1
}
sl # chr length
nonzero # nonzero region length
nonzero[[1]]/sl[[1]] # mean 46% coverge
nonzero[[2]]/sl[[2]] # mean 17%
# check MAD
library(data.table)
getMADnSD<-function(df_bg)
{
    val<-rep(df_bg$cov, df_bg$end-df_bg$start)
    MAD = median(abs(val-median(val)))
    SD = 1.4826 * MAD
    res<-c(MAD,SD)
    names(res)=c("MAD","SD")
    return(res)
}
for(input in qbgs){
    print(input)
    bg<-fread(input)
    names(bg)<-c("chr","start","end","cov")
    chrs<-unique(bg$chr)
    MS<-c("MAD","SD")
    for(i in chrs){
        MS<-rbind(MS,getMADnSD(bg[chr==i]))
    }
    MS<-rbind(MS,getMADnSD(bg))
    MS<-cbind(c("chr",chrs,"all"),MS)
    #write.table(MS, file= gsub("_.*","_MAD&SD.txt",input), row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
    print(MS)
}  # most chromosomes and overall MAD=0
# run iseg
system("mkdir peaks/iseg")
maxwl =100
minwl =10
for(input in qbgs){
    outputDir = paste0("peaks/iseg/", gsub("raimondii_|_ATAC.*","",input))
    # run iseg for BC 1, 2, 3
    cmd<-paste0("~/iSegv190207/iSeg -cz -ctp -bc 1.0-2.0-3.0-4.0-5.0-6.0 -minwl ",minwl," -maxwl ",maxwl," -d ",input," -of ",outputDir, " >",outputDir,".manifest.txt")
    message(cmd)
    system(cmd)
}
# iSEg report error, not sure how to proceed with iSeg


###-------Deeptools visualzation: ATAC BW on TSS
# Convert .bam file to a normalized bigwig using the bamCoverage tool in deepTools:
# check aggregation over TSS
cat("module load py-deeptools", file="runDeeptools.sh",sep="\n")
bamL <- list.files("mapping",pattern="*q20.bam$", full=T)
for(each in bamL){
    cmd <- paste0("bamCoverage -b ",each," -o ",gsub("bam","bw",each)," -bs=1 --normalizeUsingRPKM -p=6")
    cat(cmd, file="runDeeptools.sh",sep="\n", append=TRUE)
}
cmd="computeMatrix reference-point -S mapping/*bw -R ../refGenomes/D5.gene.bed -o mapping/TSS.gz --referencePoint TSS -b 1000 -a 1000 --skipZeros; plotHeatmap -m mapping/TSS.gz -out mapping/atac_TSS_Heatmap.png"
cat(cmd, file="runDeeptools.sh",sep="\n", append=TRUE)
system("cat runDeeptools.sh")
system("bash runDeeptools.sh")

###-------Deeptools visualzation: ATAC BW on Peak centers detected by homer and macs2
#####visualize reads on detected Homer peaks
cat("module load py-deeptools", file="runDeeptools.homer.sh",sep="\n")
cmd1="computeMatrix reference-point -S mapping/*bw -R peaks/*merged.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_homer.png --colorMap=Blues --refPointLabel=center"
cmd2="computeMatrix reference-point -S mapping/*bw -R peaks/*merged.filtered.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_homer.filter.png --colorMap=Blues --refPointLabel=center"
cat(cmd1,cmd2, file="runDeeptools.homer.sh",sep="\n", append=TRUE)
system("cat runDeeptools.homer.sh")
system("bash runDeeptools.homer.sh")
#####visualize reads on detected MACS2 peaks
cat("module load py-deeptools", file="runDeeptools.macs2.sh",sep="\n")
cmd1="computeMatrix reference-point -S mapping/*bw -R peaks/*narrowPeak -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_macs2.png --colorMap=Blues --refPointLabel=center"
cmd2="computeMatrix reference-point -S mapping/*bw -R peaks/*narrowPeak.filtered.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_macs2.filter.png --colorMap=Blues --refPointLabel=center"
cat(cmd1,cmd2, file="runDeeptools.macs2.sh",sep="\n", append=TRUE)
system("cat runDeeptools.macs2.sh")
system("bash runDeeptools.macs2.sh")

###-------Deeptools visualzation: ATAC BW on MSF and subnucleosmal centers of DNS-MNase-seq
#####visualize reads on detected Homer peaks
system("grep "20,20,255" ../isegRes/isegv1.3.2_032519/DcD_bc6.0.Fus.bed  >DcD_bc6.0.Fus.MSF.bed")
system("grep "255,20,20" ../isegRes/isegv1.3.2_032519/DcD_bc6.0.Fus.bed  >DcD_bc6.0.Fus.MRF.bed")
cat("module load py-deeptools", file="runDeeptools.MSF&L0-130.sh",sep="\n")
cmd1="computeMatrix reference-point -S mapping/*bw -R ../isegRes/MOA/iseg_v1.3.2_041219/DcL_bc4.0.Fus.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_L130.png --colorMap=Blues --refPointLabel=center"
cmd2="computeMatrix reference-point -S mapping/*bw -R DcD_bc6.0.Fus.MSF.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_msf.png --colorMap=Blues --refPointLabel=center"
cmd3="computeMatrix reference-point -S mapping/*bw -R DcD_bc6.0.Fus.MRF.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_mrf.png --colorMap=Blues --refPointLabel=center"
cat(cmd1,cmd2, file="runDeeptools.MSF&L0-130.sh",sep="\n", append=TRUE)
system("cat runDeeptools.MSF\&L0-130.sh")
system("bash runDeeptools.MSF\&L0-130.sh")

###-------Deeptools visualzation: MNase-seq BW on TSS peak centers (intersect between reps)
######check DNS profiles
cat("module load py-deeptools", file="runDeeptools.DNS.sh",sep="\n")
cmd1="computeMatrix reference-point -S ../aggregation/combo/D*full.bw -R peaks/*res.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.dnsFull_atacPeaks.png --colorMap=Blues --refPointLabel=center"
cmd2="computeMatrix reference-point -S ../aggregation/combo/D*0-130*.bw -R peaks/*res.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.dnsSizeS_atacPeaks.png --colorMap=Blues --refPointLabel=center"
cmd3="computeMatrix reference-point -S ../aggregation/combo/D*130-260*.bw -R peaks/*res.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.dnsSizeL_atacPeaks.png --colorMap=Blues --refPointLabel=center"
cmd4="computeMatrix reference-point -S ../aggregation/combo/D*mnase.bw -R peaks/*res.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.mnase_atacPeaks.png --colorMap=Blues --refPointLabel=center"
cmd5="computeMatrix reference-point -S mapping/*bw -R peaks/*res.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_atacPeaks.png --colorMap=Blues --refPointLabel=center"
cat(cmd1,cmd2,cmd3,cmd4,cmd5, file="runDeeptools.DNS.sh",sep="\n", append=TRUE)
system("cat runDeeptools.DNS.sh")
system("bash runDeeptools.DNS.sh")

--book
### Correlaion of ATAC-seq and DNS coverages
multiBigwigSummary bins -b mapping/*bw ../aggregation/D*0-130.bw ../aggregation/D*130-260.bw ../aggregation/D*full.bw ../isegRes/MOA/D*bw -out mapping.results.npz --outRawCounts mapping.results.tab
head mapping.results.tab.txt
plotCorrelation -in mapping.results.npz -o plotHeatmap_pearson.pdf --corMethod pearson --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix plotHeatmap_pearson.tab.txt


multiBigwigSummary bins -b mapping/*bw ../isegRes/MOA/D*bw ../aggregation_RPGC/D*full.bw -out mapping.results.npz --outRawCounts mapping.results.tab.txt
head mapping.results.tab.txt
plotCorrelation -in mapping.results.npz -o plotHeatmap_pearson.pdf --corMethod pearson --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix plotHeatmap_pearson.tab.txt


################
## CONCLUSION ##
################
# 1. the quality of ATAC-seq datasets is questionable: two replicates contain 21 and 16 miliion raw PE reads (93-95% mapped), corresponding to 46% and 16% of genome region being covered by non zero reads
# 2. Homer and MACS2 peak detection were performed for both replicates, and the intersection between replicates resulted in 2.8-3.0MB peak regions, about 0.4% of the 780MB genome
# 3. iSeg peak detection failed due to overall MAD=0
# 4. Aggregation profile of DNS-seq (H, L, D, L_0-130, mnase) showed tiny bumps over ATAC peak centers; hard to tell resemablance of L_0-130 with ATAC.
#
