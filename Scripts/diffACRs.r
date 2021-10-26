# module load py-deeptools r/3.5.0-py2-ufvuwmm bedtools2
# cd /work/LAS/jfw-lab/hugj2006/cottonLeaf
# mkdir Qregulation_diffACRs
# R
options(scipen=999)
sbatch=function(cfile, name="JOB", time = "7-00:00:00", N=1, n=8, partition="whatever", mem="200G",email="hugj2006@iastate.edu"){
    slurmF=paste0(cfile,".slurm")
    cat("#!/bin/bash\n#SBATCH --constraint=AVX2\n", file=slurmF,sep="\n")
    cat(paste0("#SBATCH -t ",time,"   # walltime\n#SBATCH -N ",N,"   # number of nodes in this job\n#SBATCH -n ",n,"   # total number of processor cores in this job; each node has 272 cores, Nova has 36\n#SBATCH --mem=",mem," # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G\n#SBATCH --partition=",partition," # use sinfo to get queue names\n#SBATCH -J ",name,"   # job name\n#SBATCH --output=job%j.txt\n#SBATCH --mail-type=FAIL,END\n#SBATCH --mail-user=",email,"   # email address\n"), file=slurmF, sep="\n",append=TRUE)
    cat(paste0("## LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE\necho 'Start'\ndate\n"), file=slurmF,sep="\n",append=TRUE)
    system(paste0("cat ",slurmF," ",cfile, " >>",slurmF))
    cat(paste0("\necho 'END'\ndate\n"), file=slurmF,sep="\n",append=TRUE)
    system(paste0("cat ",slurmF))
    system(paste0("sbatch ",slurmF))
}

load("Qregulation3/ACRs.rdata")->l;l
# "chr.size.A2"   "chr.size.AD1"  "chr.size.AD1a" "chr.size.AD1d"
# "chr.size.D5"   "chr.size.F1"   "chr.size.F1a"  "chr.size.F1d"
# "MSF.A2"        "MSF.AD1"       "MSF.AD1a"      "MSF.AD1d"
# "MSF.D5"        "MSF.F1"        "MSF.F1a"       "MSF.F1d"
# "MSFsum"        "SPO.A2"        "SPO.AD1"       "SPO.D5"
# "SPO.F1"        "SPOsum"
load("Qregulation3/ACRsC.rdata")->l;l
#  "ACR.A2"          "ACR.AD1a"        "ACR.AD1d"        "ACR.D5"
#  "ACR.F1a"         "ACR.F1d"         "ACRsum"          "annotateACRs"
#  "annotateACRs_TE" "chr.size.A2"     "chr.size.AD1"    "chr.size.AD1a"
#  "chr.size.AD1d"   "chr.size.D5"     "chr.size.F1"     "chr.size.F1a"
#  "chr.size.F1d"

###################################################
## Differential ACRs: map MNase-seq reads to ref ##
###################################################
library(rtracklayer)
bl =GRanges("D01", IRanges(23162000, 23820000))  # blacklist region

## TM1 ref
# export ACR Granges to GFF3 for featurecount
peaks = ACR.AD1a;peaks # 260278
a=peaks
peaks = ACR.AD1d;peaks # 189068
d=peaks
ad=c(a,d)
#
export(ad,"Qregulation_diffACRs/AD1acr.gtf",format="GTF")
ad<-read.table("Qregulation_diffACRs/AD1acr.gtf",sep="\t")
ad$V9 = paste0('acr_id "',ad$V1,'.',ad$V4,'"; acr_name "',ad$V1,'.',ad$V4,'"; ')
ad$V3="acr"
write.table(ad,file="Qregulation_diffACRs/AD1acr.gtf",sep="\t",row.names=FALSE,col.names=FALSE,quote=F)
system("head Qregulation_diffACRs/AD1acr.gtf")
#
export(ad,"Qregulation_diffACRs/AD1acr.gff3",format="GFF")
ad<-read.table("Qregulation_diffACRs/AD1acr.gff3",sep="\t")
ad$V9 = paste0("ID=",ad$V1,".",ad$V4,";Name=",ad$V1,".",ad$V4)
ad$V3="ACR"
write.table(ad,file="Qregulation_diffACRs/AD1acr.gff3",sep="\t",row.names=FALSE,col.names=FALSE,quote=F)
system("head Qregulation_diffACRs/AD1acr.gff3")

## F2020 ref
# convert MSF BED to GFF3 for featurecount
peaks = ACR.F1a;peaks # 357171
a=peaks
peaks = ACR.F1d;peaks # 224483
d=peaks
f=c(a,d)
h=findOverlaps(f,bl); length(h) # 1351 in black region
f=f[-queryHits(h)]

## A2 and D5
# convert MSF BED to GFF3 for featurecount
peaks = ACR.A2;peaks # 296312
seqlevels(peaks) = gsub("^Chr","A",seqlevels(peaks))
a=peaks
peaks = ACR.D5;peaks # 190795
seqlevels(peaks) = gsub("^Chr","D",seqlevels(peaks))
d=peaks
di=c(a,d)
h=findOverlaps(di,bl); length(h) # 0 in black region
#di=di[-queryHits(h)]

# out of 581654 F1 ACRs and 487107 A2+D5 ACRs
all<- reduce(c(f,di)) #1014201
length(h<-findOverlaps(f,di)) #52979 <10%
length(h<-findOverlaps(di,f)) #52979
#
export(all,"Qregulation_diffACRs/ADFacr.gtf",format="GTF")
ad<-read.table("Qregulation_diffACRs/ADFacr.gtf",sep="\t")
ad$V9 = paste0('acr_id "',ad$V1,'.',ad$V4,'"; acr_name "',ad$V1,'.',ad$V4,'"; ')
ad$V3="acr"
write.table(ad,file="Qregulation_diffACRs/ADFacr.gtf",sep="\t",row.names=FALSE,col.names=FALSE,quote=F)
system("head Qregulation_diffACRs/ADFacr.gtf")
#
export(all,"Qregulation_diffACRs/ADFacr.gff3",format="GFF")
all<-read.table("Qregulation_diffACRs/ADFacr.gff3",sep="\t")
tag=paste0(all$V1,".s",all$V4,".w",all$V5-all$V4+1)
all$V9 = paste0("ID=",tag,";Name=",tag)
all$V3="ACR"
write.table(all,file="Qregulation_diffACRs/ADFacr.gff3",sep="\t",row.names=FALSE,col.names=FALSE,quote=F)
system("head Qregulation_diffACRs/ADFacr.gff3")


# 1. map all reads against TM1 ref, count reads aligned to TM1 ACRs
# 2. map A2, D5, F1 reads against F1 ref, count reads aligned to combined ACRs
cfile="Qregulation_diffACRs/runHisat2_AD1.sh"
sbatch(cfile, name="diffACR_AD1", time = "4-06:00:00", N=1, n=8)
cfile="Qregulation_diffACRs/runHisat2_F1.sh"
sbatch(cfile, name="diffACR_F1", time = "4-06:00:00", N=1, n=8)


########################################
## Differential ACRs: DE using DEseq2 ##
########################################
## collect htseq-count genome table: A, F1-At, AD1-AT x 2 rep x H and L; PCA
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/Qregulation_diffACRs/AD1ref")
## examine featureCounts
count = read.table("featureCounts.txt",sep="\t",header=T)
names(count)=gsub("Qregulation_diffACRs.AD1ref.|_q20.sort.bam","",names(count))
head(count)
# Geneid Chr Start   End Strand Length A6Hn A6H A6Ln A6L D1H D1L D2H D2L F2H
# 1   A01.499 A01   499   781      +    283   40   2   95 158   0   0   0   0   9
# 2  A01.1032 A01  1032  1059      +     28   16   5   38  21   0   0   0   0   9
# 3 A01.14340 A01 14340 14380      +     41   13   1   35  32   0   0   0   0  22
# 4 A01.34920 A01 34920 34948      +     29   88   7  150 135   0   0   0   0  84
# 5 A01.40800 A01 40800 40840      +     41   27   4   29  34   0   0   0   0  13
# 6 A01.42280 A01 42280 42440      +    161    4   2   32  40   0   0   0   0   1
# F2L F3H F3L M1H M1L M2H M2L
# 1  33  14  37  60  66  41 118
# 2  14  12  14  18  13  25  35
# 3  19   4  24  15   9  18  19
# 4  76  64  78  24  15  18  34
# 5  17   9  25   8   5  16  20
# 6   6   7  10   7  16   6  75
colSums(count[,-c(1:6)])/10^6
#     A6Hn       A6H      A6Ln       A6L       D1H       D1L       D2H       D2L
# 2.620373  0.319568  5.362069  5.749111  0.869014  0.678230  0.533522  1.127835
#      F2H       F2L       F3H       F3L       M1H       M1L       M2H       M2L
# 2.994498  3.140645  2.075804  3.871355  6.539694  7.061200  5.492616 14.153271
colSums(count[grep("^A",count$Chr),-c(1:6)])/colSums(count[grep("^D",count$Chr),-c(1:6)])
#        A6Hn         A6H        A6Ln         A6L         D1H         D1L
# 33.80419450 30.04108791 22.20026739 20.81726449  0.09822997  0.08001331
#        D2H         D2L         F2H         F2L         F3H         F3L
# 0.11106668  0.08011446  1.78347704  1.87981445  1.74300771  1.84671643
#        M1H         M1L         M2H         M2L
# 1.23996797  1.38528428  1.21924955  1.34962411
# diploid reads were mapped to AD1 ref genome, clearly A2 reads mapped mainly to At, vice versa; for F1 and AD1, twice mapped to At as to Dt, because more ACRs detected in At.

## table preprocessing, normalization by overall library size, PCA plot
x=count[,-c(1:6)]; rownames(x)=count$Geneid
dim(x)   # 449346     16
rpm<-sweep(x,2,colSums(x),"/")*10^6
# my custom PCA function
library(ggplot2)
library(scales)
library(ape)
plotGrouping <- function(norm_log, color, shape, text, tip, save = "plotGrouping.pdf"){
    # norm<-sweep(total,2,info$lib_size,"/")*10^6
    # norm_log <- log2(norm+1)
    pca=prcomp(t(norm_log))
    dat = as.data.frame(pca$x)
    proportion<-summary(pca)$importance[2,1:2]*100
    proportion<-paste0(names(proportion)," (", proportion, "%)")
    p<-ggplot(aes(PC1, PC2, color=color, shape=shape,label=tip),data=dat) + geom_point() +xlab(proportion[1]) + ylab(proportion[2])
    pdf(save)
    print( p + geom_text(color="grey", hjust = 0, nudge_x = 0.09) )
    
    hc<-hclust( dist(t(norm_log)) )
    tre <- as.phylo(hc)
    tre$tip.label <- as.character(tip)
    # try to match ggplot color: library(scale);
    # show_col(col4<-hue_pal()(4))
    tipCol <- as.factor(color)
    levels(tipCol) <- hue_pal()(nlevels(tipCol))
    plot(tre,tip.col =as.character(tipCol),  type="unrooted",cex=0.6, no.margin=TRUE)
    plot(tre,tip.col =as.character(tipCol), type="fan",cex=0.6, no.margin=TRUE)
    dev.off()
}
# plots
pdf("checkACRcounts.pdf")
boxplot(count[,-c(1:6)])
boxplot(log2(count[,-c(1:6)]+1))
boxplot(log2(rpm+1))
dev.off()
## PCA and cladogram plots
info = data.frame(sample=names(x), genome = substring(names(x),1,1), digestion = substring(names(x),3,3),acr_size=colSums(x))
info$lib_size =c(94034291,11686216,114138711,107018821,21133900,11427993,13290606,20095084,93246674,59574178,61991464,76845486,156023849,94961844,138165016,190402124)  #q20 mapped
info$acr_size/info$lib_size
plotGrouping(log2(rpm+1), color=info$genome, shape=info$digestion, tip=info$sample, text=info$genome, save = "plotGrouping.log2rpm.pdf")

## Examine the changes across genome M vs F
# unify the total reads mapped to ref
library(DESeq2)
y=x[,grep("M|F",names(x))]
coldata=info[grep("M|F",info$sample),]
dds <- DESeqDataSetFromMatrix( countData = y, colData = coldata, design =~genome+digestion+digestion:genome )
sizeFactors(dds) <-coldata$lib_size
print(sizeFactors(dds))
levels(dds$genome)=c("M","F") # M as ref
levels(dds$digestion)=c("H","L") # H as ref
ddsMF=DESeq(dds)
resultsNames(ddsMF) #  "Intercept"          "genome_F_vs_M"      "digestion_L_vs_H"   "genomeF.digestionL"
# the condition effect for ref genotype M (the main effect)
print( summary(results(ddsMF, contrast=c("digestion","L","H")),alpha=.05) )  #0+4
print( summary(results(ddsMF, name="digestion_L_vs_H"),alpha=.05) ) # 0+4; clearly wrong because all 445822 were detected as ACRs with L>H
# the condition effect for genotype L: main effect *plus* the interaction term， (the extra condition effect in genotype III compared to genotype I).
print( summary(results(ddsMF, contrast=list( c("digestion_L_vs_H","genomeF.digestionL") )),alpha=.05) ) #129+9
# the interaction term for condition effect in genotype F vs genotype M，tests if the condition effect is different in III compared to I
print( summary(results(ddsMF, name="genomeF.digestionL"),alpha=.05) )  #0+0
# LTR model testing, again showing no interaction
ddsLRT <- DESeq(ddsMF, test="LRT", reduced= ~ genome+digestion)
summary(results(ddsLRT),alpha=.05)  #0+0

## Conclusion ##
# I already used total mapped reads of each sample as lib_size
# This approach is flawed as only counts mapped to ACRs were extracted for comparison, while all ACRs should have significant L>H; normalization was impossible


##------------------OLD METHOD BELOW, ABANDON ------------------------------##


########################################
## Differential Peak Analysis by csaw ##
########################################
# 1. map all reads against F2020=A2+D5 ref
# cfile="runHisat2.sh"
# sbatch(cfile, name="differentialACR", time = "4-06:00:00", N=1, n=8)
# source("diffACRs_csaw_Fref.r")
# 2. map all reads against AD1 ref
# cfile="runHisat2.sh"
# sbatch(cfile, name="differentialACR", time = "4-06:00:00", N=1, n=8)
# source("diffACRs_csaw_AD1ref.r")

#----------------------------------------------------------------------------------
BiocManager::install("csaw")
library(csaw)

1. Loading in data from BAM files.
param <- readParam(minq=20)
data <- windowCounts(bam.files, ext=110, width=10, param=param)
2. Filtering out uninteresting regions.
library(edgeR)
keep <- aveLogCPM(asDGEList(data)) >= -1
data <- data[keep,]
3. Calculating normalization factors.
binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
data <- normFactors(binned, se.out=data)
4. Identifying DB windows.
y <- asDGEList(data)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
results <- glmQLFTest(fit)
5. Correcting for multiple testing.
merged <- mergeResults(data, results$table, tol=1000L)
#----------------------------------------------------------------------------------

for j in $( ls *bam ); do
echo ''
echo $j
samtools sort $j -o ${j%%_*}_q20.f.sort.bam; samtools index ${j%%_*}_q20.f.sort.bam
done

setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/Qregulation_diffACRs/Fref")
library(csaw)
library(edgeR)

bam.files=list.files(patter="sort.bam$")
df=data.frame(BAM=bam.files, genome=substring(bam.files,1,1),digestion=substring(bam.files,3,3))
# separate A and D chromosomes
for(g in c("A","D"))
{
   ## 1. input sorted and indexed bam files
   use = df[df$genome=="F"|df$genome==g,]
   use.chr <- c(paste0(g,"0", 1:9), paste0(g,10:13))
   print(use)
   print(use.chr)
   # minimal MAPQ quality 20, paired-end
   param <- readParam(pe="both",restrict=use.chr);param
   # Counting reads into windows
   bams=as.character(use$BAM)
   win.data <- windowCounts(bams, ext=150, width=150, param=param)
   win.data
   head(assay(win.data))
   
   ## 2. Filtering out uninteresting regions.
   # simple restriction to filter window abundances >-1
   abundances <- aveLogCPM(asDGEList(win.data))
   summary(abundances)
   keep.simple <- abundances > -1
   filtered.data <- win.data[keep.simple,]
   summary(keep.simple)
   
   ## 3. Calculating normalization factors.
   binned <- windowCounts(bams, bin=TRUE, width=10000, param=param)
   #data <- normFactors(binned, se.out=data)
   filtered.data <- normFactors(binned, se.out=filtered.data)
   filtered.data$norm.factors
   # MA plots
   pdf(paste0("normalizationOutcome_",g,".pdf"))
   par(mfrow=c(3, 3), mar=c(5, 4, 2, 1.5))
   adj.counts <- cpm(asDGEList(binned), log=TRUE)
   normfacs <- filtered.data$norm.factors
   for (i in seq_len(length(bams)-1)) {
       cur.x <- adj.counts[,1]
       cur.y <- adj.counts[,1+i]
       smoothScatter(x=(cur.x+cur.y)/2+6*log2(10), y=cur.x-cur.y, xlab="A", ylab="M", main=paste("1 vs", i+1))
       all.dist <- diff(log2(normfacs[c(i+1, 1)]))
       abline(h=all.dist, col="red")
   }
   dev.off()

   ## 4. Identifying DB windows.
   # one-way layout with 4 groups
   group=factor(paste0(use$genome,use$digestion))
   design <- model.matrix(~0+group)
   #genome = factor(use$genome,levels=c("F",g))
   #digestion = factor(use$digestion, levels=c("H","L"))
   #design=model.matrix(~genome+digestion+genome:digestion)
   y <- asDGEList(filtered.data,group=group)
   # examine replicats through multi-dimensional scaling (MDS) plots.
   pdf(paste0("MDSplot_",g,".pdf"))
   par(mfrow=c(2,2), mar=c(5,4,2,2))
   adj.counts <- cpm(y, log=TRUE)
   for (top in c(100, 500, 1000, 5000)) {
   out <- plotMDS(adj.counts, main=top, labels=group, top=top)
   }
   dev.off()
   # Estimate dispersion, plot e coefficient of biological variation (BCV) and quasi-likelihood (QL)dispersion
   y <- estimateDisp(y, design, robust=TRUE)
   fit <- glmQLFit(y, design, robust=TRUE)
      #pdf(paste0("Dispersion_",g,".pdf"))
      #plotBCV(y)
      #plotQLDisp(fit)
      #dev.off()
   # Differential Binding analyses
   con <- makeContrasts(groupAL-groupAH,levels=design)
   resDi <- glmQLFTest(fit, contrast=con)
   con <- makeContrasts(groupFL-groupFH,levels=design)
   resF <- glmQLFTest(fit, contrast=con)
   con <- makeContrasts((groupFL-groupFH)-(groupAL-groupAH),levels=design)
   resH <- glmQLFTest(fit, contrast=con)
   summary(decideTests(resA)) # FDR <0.05
   summary(decideTests(resF)) # FDR <0.05
   summary(decideTests(resH)) # FDR <0.05
   rowData(filtered.data) <- cbind(rowData(filtered.data), resA$table,resF$table,resH$table)
   
   save(resDi,resF,resH, filtered.data, file=paste0("DBres_",g,".rdata"))
}


   ## 5. Correcting for multiple testing.
   mergeWindows followed by combineTests and getBestTests.
   merge.res <- mergeResults(filtered.data, res$table, tol=100, merge.args=list(max.width=5000))
   names(merge.res)
   
   
   merged <- mergeResults(filtered.data, results$table, tol=1000L)

   
rowRanges(filtered.data)
   
   
   
   
}


## 2. Filtering out uninteresting regions.
# 2.1 for histone modification: the abundance of the window is greater than the background abundance estimate by log2(3) or more, not suitable for widespread MNase-seq peaks
bins <- windowCounts(bams, bin=TRUE, width=2000, param=param)
filter.stat <- filterWindows(win.data, bins, type="global")
quantile(filter.stat$filter)  # check this!!!!
min.fc <- 3
keep <- filter.stat$filter > log2(min.fc)
summary(keep) # most will be filtered out, not right
pdf("checkBackground.pdf")
plot(hist(filter.stat$back.abundances, main="", breaks=50, xlab="Background abundance (log2-CPM)"))
threshold <- filter.stat$abundances[1] - filter.stat$filter[1] + log2(min.fc)
abline(v=threshold, col="red")
dev.off()
filtered.data <- win.data[keep,]
# 2.2 Normalizing for library-specific trended biases
filtered.data <- win.data
win.ab <- filter.stat$abundances
offsets <- normOffsets(filtered.data, type="loess")
head(offsets)
# before offset
adjc <- log2(assay(filtered.data)+0.5)
logfc <- adjc[,1] - adjc[,2]
smoothScatter(win.ab, logfc, ylim=c(-6, 6), xlim=c(0, 5), xlab="Average abundance", ylab="Log-fold change")
# after offset
norm.adjc <- adjc - offsets/log(2)
norm.fc <- norm.adjc[,1]-norm.adjc[,4]
smoothScatter(win.ab, norm.fc, ylim=c(-6, 6), xlim=c(0, 5), xlab="Average abundance", ylab="Log-fold change")
dev.off()
# 2.3 edgeR filter
library(edgeR)
y <- asDGEList(filtered.data)
y$offset <- offsets
y <- estimateDisp(y, design) summary(y$trended.dispersion)

keep <- aveLogCPM(asDGEList(data)) >= -1
data <- data[keep,]




