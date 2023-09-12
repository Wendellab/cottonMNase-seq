########################################################
################ Regulation question 3 #################
########################################################
# Are ACRs conserved across genomes/subgenomes

# module load py-deeptools r/3.5.0-py2-ufvuwmm bedtools2
# cd /work/LAS/jfw-lab/hugj2006/cottonLeaf
# R
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

peakSumByChr=function(gr,chr_size){
    require(GenomicRanges)
    df=chr_size
    names(df)=c("chr","chr.length")
    rownames(df)=df$chr
    ## number of peaks
    nPeaks = tapply(as.character(seqnames(gr)), list(seqnames(gr)),length)
    df[names(nPeaks),"number"]= nPeaks
    ## Region bps called as peaks
    sRegion = tapply(width(gr), list(seqnames(gr)),sum)
    df[names(sRegion),"size"]= sRegion
    ## Proportion of regions called
    df$sizePerc=df$size/df$chr.length*100
    return(df)
}
###########################
## Deeptools correlation ##
###########################
# Qregulation1 already generate correlations between ATAC, SPO, SPfrenter, and DNS


###################
## Assemble ACRs ##
###################
library(rtracklayer)
library(GenomicAlignments)
library(GenomicFeatures)
library(plot.matrix)
library(gplots)

# import chr size
chr.size.A2 = read.table("mappingA_WHU/chr.size.txt",sep="\t")
chr.size.D5 = read.table("mappingD/chr.size.txt",sep="\t")
chr.size.AD1 = read.table("mappingM_UTX/chr.size.txt",sep="\t")
chr.size.F1 = read.table("mappingF2020/chr.size.txt",sep="\t")
chr.size.AD1a = chr.size.AD1[1:13, ]; chr.size.AD1a$V1 = gsub("A","Chr",chr.size.AD1a$V1)
chr.size.AD1d = chr.size.AD1[14:26,]; chr.size.AD1d$V1 = gsub("D","Chr",chr.size.AD1d$V1)
chr.size.F1a  = chr.size.A2
chr.size.F1d  = chr.size.D5

### DNase-seq DHS macs2, combine
DHS.A2.c = import("DNase/A2/macs2.combine.bed", format="bed")
DHS.D5.c = import("DNase/D5/macs2.combine.bed", format="bed")
DHS.AD1.c = import("DNase/AD1/macs2.combine.bed", format="bed")
DHSsum.c = list(D5=peakSumByChr(DHS.D5.c,chr.size.D5), A2=peakSumByChr(DHS.A2.c,chr.size.A2), AD1=peakSumByChr(DHS.AD1.c,chr.size.AD1))
# D5 1-2%, A2 0.9-1.6%, AD1 0.4-1%
write.table(as.data.frameDHSsum.c,file="Qregulation3/DHS_combined_sum.txt", sep="\t")

### DNase-seq DHS macs2, intersect
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")
# D5
peaks0 = import("DNase/D5/macs2.intersect.bed", format="bed",extraCols = extraCols_narrowPeak)
peaks0=peaks0[width(peaks0)>1,]
peaks0 #12960
DHS.D5=peaks0
# A2
peaks0 = import("DNase/A2/macs2.intersect.bed", format="bed",extraCols = extraCols_narrowPeak)
peaks0=peaks0[width(peaks0)>1,]
peaks0 #2200, quite wrong
DHS.A2=peaks0
# AD1
peaks0 = import("DNase/AD1/macs2.intersect.bed", format="bed",extraCols = extraCols_narrowPeak)
peaks0=peaks0[width(peaks0)>1,]
peaks0 #15471, also lower than D5
DHS.AD1=peaks0
# sum
DHSsum = list(D5=peakSumByChr(DHS.D5,chr.size.D5), A2=peakSumByChr(DHS.A2,chr.size.A2),
AD1=peakSumByChr(DHS.AD1,chr.size.AD1))
write.table(DHSsum,file="Qregulation3/DHS_intersect_sum.txt", sep="\t")

### MSFs
#
peaks0 = import("isegRes/iseg_v1.3.4_083120/DcD_bc6.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]; peaks # 63282
bl =GRanges("Chr01", IRanges(23162000, 23820000))  # blacklist region
h=findOverlaps(peaks,bl); length(h) # 1241 in black region
MSF.D5=peaks[-queryHits(h)]# 62041
peaks0 = import("isegRes/iseg_v1.3.4_083120/A6Dn_bc6.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]; peaks # 128707
MSF.A2=peaks
#
peaks0 = import("isegRes/iseg_v1.3.4_083120/FcDd_bc6.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]; peaks # 60754
bl =GRanges("Chr01", IRanges(23162000, 23820000))  # blacklist region
h=findOverlaps(peaks,bl); length(h) # 333 in black region
MSF.F1d=peaks[-queryHits(h)]# 60421
peaks0 = import("isegRes/iseg_v1.3.4_083120/FcDa_bc6.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]; peaks # 128417
MSF.F1a=peaks
#
peaks0 = import("isegRes/iseg_v1.3.4_083120/FcD_bc6.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]; peaks # 190187
MSF.F1=peaks
#
peaks0 = import("isegRes/iseg_v1.3.4_083120/McD_bc6.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]; peaks # 156668
MSF.AD1=peaks
#
peaks0 = import("isegRes/iseg_v1.3.4_083120/McDa_bc6.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]; peaks # 98722
MSF.AD1a=peaks
#
peaks0 = import("isegRes/iseg_v1.3.4_083120/McDd_bc6.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]; peaks # 58684
MSF.AD1d=peaks

#sum
MSFsum = list(D5=peakSumByChr(MSF.D5,chr.size.D5), A2=peakSumByChr(MSF.A2,chr.size.A2), F1d=peakSumByChr(MSF.F1d,chr.size.D5), F1a=peakSumByChr(MSF.F1a,chr.size.A2), F1=peakSumByChr(MSF.F1,chr.size.F1), AD1=peakSumByChr(MSF.AD1,chr.size.AD1), AD1d=peakSumByChr(MSF.AD1d,chr.size.AD1d), AD1a=peakSumByChr(MSF.AD1a,chr.size.AD1a)
)
write.table(MSFsum,file="Qregulation3/MSF_sum.txt", sep="\t")
#save(list=grep("MSF|chr[.]size",ls(),value=T), file="Qregulation3/ACRs_MSF.rdata")

# SPO Fr
peaks0 = import("isegRes/SPO/iseg_v1.3.4_060520/D_frenters_bc5.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]
SPO.D5=peaks # 140268
peaks0 = import("isegRes/SPO/iseg_v1.3.4_060520/A_frenters_bc4.5.Fus.bed",format="bed") # 2/24/21 change from 4 to 4.5
peaks = peaks0[peaks0$itemRgb=="#1414FF"]
SPO.A2=peaks # 324705 to 186730
peaks0 = import("isegRes/SPO/iseg_v1.3.4_060520/F_frenters_bc4.5.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]
SPO.F1=peaks # 421763
peaks0 = import("isegRes/SPO/iseg_v1.3.4_060520/M_frenters_bc4.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]
SPO.AD1=peaks # 347908
SPOsum = list(D5=peakSumByChr(SPO.D5,chr.size.D5), A2=peakSumByChr(SPO.A2,chr.size.A2), F1=peakSumByChr(SPO.F1,chr.size.F1),AD1=peakSumByChr(SPO.AD1,chr.size.AD1))
write.table(as.data.frame(SPOsum),file="Qregulation3/SPO_sum.txt", sep="\t")
save(list=grep("SPO|MSF|chr[.]size",ls(),value=T), file="Qregulation3/ACRs.rdata")

# combine MSF and SPO
ACR.A2=reduce(c(MSF.A2,SPO.A2))
ACR.D5=reduce(c(MSF.D5,SPO.D5))
seqlevels(MSF.F1d) = gsub("Chr","D",seqlevels(MSF.F1d))
seqlevels(MSF.F1a) = gsub("Chr","A",seqlevels(MSF.F1a))
seqlevels(MSF.AD1a) = gsub("Chr","A",seqlevels(MSF.AD1a))
seqlevels(MSF.AD1d) = gsub("Chr","D",seqlevels(MSF.AD1d))
ACR.F1a=reduce(c(MSF.F1a,SPO.F1[grep("A",seqnames(SPO.F1)),]))
ACR.F1d=reduce(c(MSF.F1d,SPO.F1[grep("D",seqnames(SPO.F1)),]))
ACR.AD1a=reduce(c(MSF.AD1a,SPO.AD1[grep("A",seqnames(SPO.AD1)),]))
ACR.AD1d=reduce(c(MSF.AD1d,SPO.AD1[grep("D",seqnames(SPO.AD1)),]))
ACRsum = list(D5=peakSumByChr(ACR.D5,chr.size.D5), A2=peakSumByChr(ACR.A2,chr.size.A2), F1=peakSumByChr(c(ACR.F1a,ACR.F1d),chr.size.F1),AD1=peakSumByChr(c(ACR.AD1a,ACR.AD1d),chr.size.AD1))
write.table(as.data.frame(ACRsum),file="Qregulation3/ACR_sum.txt", sep="\t")
save(list=grep("ACR|chr[.]size",ls(),value=T), file="Qregulation3/ACRsC.rdata")


####################
## ACR annotation ##
####################
# module load r-udunits2
library(ChIPseeker) # module load r-udunits2
library(GenomicFeatures)
library(ggthemes)
library(gridExtra)
options(stringsAsFactors=FALSE)
load("Qregulation3/ACRs.rdata")->l;l

### import genome annotation
#
txdb <- p("refGenomes/txdb.D5.sqlite")
gns= genes(txdb)
gns.D5= gns[grep("Chr",seqnames(gns))]# 37223
#
txdb <- loadDb("refGenomes/txdb.A2WHU.sqlite")
gns= genes(txdb)
gns.A2= gns[grep("Chr",seqnames(gns))]# 41739
#
txdb <- loadDb("refGenomes/txdb.AD1utx.sqlite")
gns= genes(txdb)
gns.AD1= gns[grep("scaffold",seqnames(gns),invert=T)]# 74902
gns.AD1a = gns.AD1[grep("A",seqnames(gns.AD1)),]; seqlevels(gns.AD1a)=gsub("A","Chr",seqlevels(gns.AD1a))
gns.AD1d = gns.AD1[grep("D",seqnames(gns.AD1)),]; seqlevels(gns.AD1d)=gsub("D","Chr",seqlevels(gns.AD1d));
#
txdb <- loadDb("refGenomes/txdb.F2020.sqlite")
gns= genes(txdb)
gns.F1= gns[grep("scaffold|Contig",seqnames(gns),invert=T)]# 78962
gns.F1a = gns.F1[grep("A",seqnames(gns.F1)),];
gns.F1d = gns.F1[grep("D",seqnames(gns.F1)),];

ss=data.frame(genome=c("A2","D5","F1a","F1d","F1-At","F1-Dt","AD1a","AD1d","AD1-At","AD1-Dt"), gsize=c(sum(chr.size.A2$V2),sum(chr.size.D5$V2),sum(chr.size.F1$V2[1:13]),sum(chr.size.F1$V2[14:26]),sum(chr.size.F1$V2[1:13]),sum(chr.size.F1$V2[14:26]),sum(chr.size.AD1$V2[1:13]),sum(chr.size.AD1$V2[14:26]),sum(chr.size.AD1$V2[1:13]),sum(chr.size.AD1$V2[14:26])), geneN=c(length(gns.A2),length(gns.D5),length(gns.F1[grep("A",seqnames(gns.F1))]),length(gns.F1[grep("D",seqnames(gns.F1))]),length(gns.F1[grep("A",seqnames(gns.F1))]),length(gns.F1[grep("D",seqnames(gns.F1))]),length(gns.AD1[grep("A",seqnames(gns.AD1))]),length(gns.AD1[grep("D",seqnames(gns.AD1))]),length(gns.AD1[grep("A",seqnames(gns.AD1))]),length(gns.AD1[grep("D",seqnames(gns.AD1))]))
)
fls=data.frame(genome=c("A2","D5","F1a","F1d","F1","AD1a","AD1d","AD1"),tag1=c("A6Dn","DcD","FcDa","FcDd","FcD","McDa","McDd","McD"),tag2=c("A","D","Fa","Fd","F","Ma","Md","M"),chr=c("chr.size.A2","chr.size.D5","chr.size.A2","chr.size.D5","chr.size.F1","chr.size.AD1a","chr.size.AD1d","chr.size.AD1"), chrPath=c("mappingA_WHU/chr.size.txt","mappingD/chr.size.txt","mappingA_WHU/chr.size.txt","mappingD/chr.size.txt","mappingF2020/chr.size.txt","mappingM_UTX/chr.size.A.txt","mappingM_UTX/chr.size.D.txt", "mappingM_UTX/chr.size.txt"), gns=c("gns.A2","gns.D5","gns.F1","gns.F1","gns.F1","gns.AD1a","gns.AD1d","gns.AD1"),ref=c("refGenomes/A2WHU_13.fasta","refGenomes/Dgenome2_13.fasta","refGenomes/A2WHU_13.fasta","refGenomes/Dgenome2_13.fasta","refGenomes/F2020_26.fasta", "refGenomes/TM1utx_26.fasta","refGenomes/TM1utx_26.fasta","refGenomes/TM1utx_26.fasta"))
rm(res)
### g-p-d categorization
source("FUN.acr.r") # ann


#### MSF
# 4 genomes, 6 BC, annotate ACR, summaize length, number, percentage genome; 126 rows
for(g in 1:8){
    genome=fls$genome[g]
    print(genome)
    for(bc in c("4.0","4.5","5.0","5.5","6.0","6.5","7.0"))
    {
        print(bc)
        f=paste0("isegRes/iseg_v1.3.4_083120/",fls$tag1[g],"_bc",bc,".Fus.bed")
        peaks0 = import(f,format="bed")
        peaks = peaks0[peaks0$itemRgb=="#1414FF"] # 68713
        if(genome=="D5"|genome=="F1d"){
            bl =GRanges("Chr01", IRanges(23162000, 23820000))  # blacklist region
            h=findOverlaps(peaks,bl); length(h) # 1174 in black region
            peaks=peaks[-queryHits(h)]# 67539
        }
        an=annotateACRs(peaks, get(fls$gns[g]), distance=2000)
        if(genome %in% c("F1","AD1")){
            rrA=getPeakSum(an[grep("A",an$seqnames),],get(paste0("chr.size.",genome))[1:13,])
            rrA$BC = bc
            rrA$genome=paste0(genome,"-At")
            rrD=getPeakSum(an[grep("D",an$seqnames),],get(paste0("chr.size.",genome))[14:26,])
            rrD$BC = bc
            rrD$genome=paste0(genome,"-Dt")
            res= rbind(res,rrA, rrD)
        }else
        {
            rr=getPeakSum(an,get(paste0("chr.size.",genome)))
            rr$BC = bc
            rr$genome=genome
            if(exists("res")){res= rbind(res,rr)}else{res=rr}
        }
        
    }
}
write.table(res,sep="\t",file="Qregulation3/compareMSFs.txt")

res=read.table("Qregulation3/compareMSFs.txt",sep="\t",header=T)
res$genome=factor(res$genome, levels=c("A2","D5","F1a","F1d","F1-At","F1-Dt","AD1a","AD1d","AD1-At","AD1-Dt"))
#plots
pdf("Qregulation3/compareMSFs_by_genome.pdf")
mycols=few_pal()(3)
### Does BC stringency affect categorization? NO
ggplot(res, aes(x = BC, y = sizeP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=5) + labs(title="Length Proportion")
ggplot(res, aes(x = BC, y = numberP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=5) + labs(title="Number Proportion")
ggplot(res, aes(x = BC, y = size/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=5) + labs(title="Length (Mbp)")
ggplot(res, aes(x = BC, y = number/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=5) + labs(title="Number (million)")
ggplot(res, aes(x = BC, y = genomeP*100, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=5) + labs(title="Genome Percentage")
dev.off()
### Does BC stringency affect categorization patterns across genome? NO
pdf("Qregulation3/compareMSFs_by_BC.pdf")
ggplot(res, aes(x = genome, y = sizeP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Length Proportion")+ coord_flip()
ggplot(res, aes(x = genome, y = numberP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Number Proportion")+ coord_flip()
ggplot(res, aes(x = genome, y = size/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Length (Mbp)")+ coord_flip()
ggplot(res, aes(x = genome, y = number/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Number (million)")+ coord_flip()
ggplot(res, aes(x = genome, y = genomeP*100, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Genome Percentage")+ coord_flip()
dev.off()
# focus on BC=6.5
res0=res
res=res0[res0$BC=="6",]
pdf("Qregulation3/compareMSFs_BC6.0.pdf")
ggplot(res, aes(x = genome, y = sizeP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Length Proportion")+ coord_flip()
ggplot(res, aes(x = genome, y = numberP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Number Proportion")+ coord_flip()
ggplot(res, aes(x = genome, y = size/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Length (Mbp)")+ coord_flip()
ggplot(res, aes(x = genome, y = number/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Number (million)")+ coord_flip()
ggplot(res, aes(x = genome, y = genomeP*100, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Genome Percentage")+ coord_flip()
dev.off()
# scale with genome size
res=merge(res,ss,by="genome")
res$ploidy="Diploid"
res$ploidy[grep("F1",res$genome)]="F1"
res$ploidy[grep("AD1",res$genome)]="AD1"
res$ploidy=factor(res$ploidy,levels=c("Diploid","F1","AD1"))
write.table(res,"Qregulation3/compareMSFs_BC6.0.txt",sep="\t")
# plot genome size at X-axis
pdf("Qregulation3/compareMSFs_XY_BC6.0_genome.pdf")
resG=res[res$genome%in%c("A2","D5","F1-At","F1-Dt","AD1-At","AD1-Dt"),]
p1=ggplot(resG, aes(x=gsize/10^9, y=number/10^6, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p2=ggplot(resG, aes(x=gsize/10^9, y=size/10^6, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p3=ggplot(resG, aes(x=gsize/10^9, y=numberP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p4=ggplot(resG, aes(x=gsize/10^9, y=sizeP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p5=ggplot(resG, aes(x=gsize/10^9, y=genomeP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
grid.arrange( p3 + theme(legend.position = "none",axis.text = element_text(size=8)), p1+theme(legend.position = "none",axis.text = element_text(size=8)), p2+theme(legend.position="none",axis.text = element_text(size=8)), p4+theme(legend.position="none",axis.text = element_text(size=8)), nrow = 2, ncol =2,layout_matrix = rbind(c(1,3), c(2,4)))
p1
p2
p3
p4
p5
dev.off()
pdf("Qregulation3/compareMSFs_XY_BC6.0_subOnly.pdf")
resS=res[res$genome%in%c("A2","D5","F1a","F1d","AD1a","AD1d"),]
p1=ggplot(resS, aes(x=gsize/10^9, y=number/10^6, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p2=ggplot(resS, aes(x=gsize/10^9, y=size/10^6, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p3=ggplot(resS, aes(x=gsize/10^9, y=numberP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p4=ggplot(resS, aes(x=gsize/10^9, y=sizeP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p5=ggplot(resS, aes(x=gsize/10^9, y=genomeP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
grid.arrange( p3 + theme(legend.position = "none",axis.text = element_text(size=8)), p1+theme(legend.position = "none",axis.text = element_text(size=8)), p2+theme(legend.position="none",axis.text = element_text(size=8)), p4+theme(legend.position="none",axis.text = element_text(size=8)), nrow = 2, ncol =2,layout_matrix = rbind(c(1,3), c(2,4)))
p1
p2
p3
p4
p5
dev.off()

# stats: more g and p ACRs in AD1 than A2/F1?
library("agricolae")
#
res=resG
summary(a<-aov(lm(size~ploidy,res[res$type=="Distal",])))
summary(a<-aov(lm(size~ploidy,res[res$type=="Genic",])))
summary(a<-aov(lm(size~ploidy,res[res$type=="Proximal",]))) #P=0.0339
TukeyHSD(a)
duncan.test(a,"ploidy")
print(duncan.test(a,"ploidy"))  # AD1>Diploids>F1
summary(a<-aov(lm(numberP~ploidy,res[res$type=="Distal",])))
summary(a<-aov(lm(numberP~ploidy,res[res$type=="Genic",])))
summary(a<-aov(lm(numberP~ploidy,res[res$type=="Proximal",]))) #0.0915
TukeyHSD(a)
duncan.test(a,"ploidy")
print(duncan.test(a,"ploidy"))  # AD1=Diploids=F1
#
res=resS
summary(a<-aov(lm(size~ploidy,res[res$type=="Distal",])))
summary(a<-aov(lm(size~ploidy,res[res$type=="Genic",])))
summary(a<-aov(lm(size~ploidy,res[res$type=="Proximal",]))) #P=0.123
TukeyHSD(a)
duncan.test(a,"ploidy")
print(duncan.test(a,"ploidy"))  # AD1>Diploids>F1
summary(a<-aov(lm(numberP~ploidy,res[res$type=="Distal",])))
summary(a<-aov(lm(numberP~ploidy,res[res$type=="Genic",])))
summary(a<-aov(lm(numberP~ploidy,res[res$type=="Proximal",]))) #0.0997
TukeyHSD(a)
duncan.test(a,"ploidy")
print(duncan.test(a,"ploidy"))  # AD1=Diploids=F1


#### SPO
# 4 genomes, 7 BC, annotate ACR, summaize length, number, percentage genome; 126 rows
for(g in c("A2","D5","F1","AD1")){
    print(g)
    for(bc in c("4.0","4.5","5.0","5.5","6.0","6.5"))
    {
        print(bc)
        f=paste0("isegRes/SPO/iseg_v1.3.4_060520/",fls$tag2[fls$genome==g],"_frenters_bc",bc,".Fus.bed")
        peaks0 = import(f,format="bed")
        peaks = peaks0[peaks0$itemRgb=="#1414FF"] # 290381
        if(g=="D5"){
            bl =GRanges("Chr01", IRanges(23162000, 23820000))  # blacklist region
            h=findOverlaps(peaks,bl);
            if(length(h)>0){ # 1174 in black region
                peaks=peaks[-queryHits(h)]}# 67539
        }
        an=annotateACRs(peaks, get(fls$gns[fls$genome==g]), distance=2000)
        if(g %in% c("A2","D5")){
            rr=getPeakSum(an,get(fls$chr[fls$genome==g]))
            rr$BC = bc
            rr$genome=g
            if(exists("res")){res= rbind(res,rr)}else{res=rr}
        }
        if(g %in% c("F1","AD1")){
            rrA=getPeakSum(an[grep("A",an$seqnames),],get(paste0("chr.size.",g))[1:13,])
            rrA$BC = bc
            rrA$genome=paste0(g,"-At")
            rrD=getPeakSum(an[grep("D",an$seqnames),],get(paste0("chr.size.",g))[14:26,])
            rrD$BC = bc
            rrD$genome=paste0(g,"-Dt")
            res= rbind(res,rrA, rrD)
        }
    }
}
write.table(res,sep="\t",file="Qregulation3/compareSPOs.txt")
pdf("Qregulation3/compareSPOs_by_genome.pdf")
mycols=few_pal()(3)
### Does BC stringency affect categorization? NO
ggplot(res, aes(x = BC, y = sizeP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=3) + labs(title="Length Proportion")
ggplot(res, aes(x = BC, y = numberP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=3) + labs(title="Number Proportion")
ggplot(res, aes(x = BC, y = size/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=3) + labs(title="Length (Mbp)")
ggplot(res, aes(x = BC, y = number/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=3) + labs(title="Number (million)")
ggplot(res, aes(x = BC, y = genomeP*100, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=3) + labs(title="Genome Percentage")
dev.off()
pdf("Qregulation3/compareSPOs_by_BC.pdf")
### Does BC stringency affect categorization patterns across genome? NO
res$genome=factor(res$genome,levels=rev(c("A2","D5","F1-At","F1-Dt","AD1-At","AD1-Dt")))
ggplot(res, aes(x = genome, y = sizeP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Length Proportion")+ coord_flip()
ggplot(res, aes(x = genome, y = numberP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Number Proportion")+ coord_flip()
ggplot(res, aes(x = genome, y = size/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Length (Mbp)")+ coord_flip()
ggplot(res, aes(x = genome, y = number/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Number (million)")+ coord_flip()
ggplot(res, aes(x = genome, y = genomeP*100, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Genome Percentage")+ coord_flip()
dev.off()
with(res, tapply(genomeP,list(genome,BC),sum))
#                4.0         4.5         5.0         5.5          6.0
# AD1-Dt 0.008055880 0.004947770 0.003104047 0.001997496 0.0012162494
# AD1-At 0.005842111 0.003494750 0.002180364 0.001396616 0.0008430932
# F1-Dt  0.016190405 0.009758967 0.006572444 0.003994044 0.0025905409
# F1-At  0.011255527 0.006608217 0.004385114 0.002573342 0.0015794952
# D5     0.016853256 0.011507408 0.006763354 0.003648476 0.0024885306
# A2     0.008330607 0.004799906 0.002591990 0.001617232 0.0009708739
#                 6.5
# AD1-Dt 0.0008370619
# AD1-At 0.0005803150
# F1-Dt  0.0017958310
# F1-At  0.0010366324
# D5     0.0015366295
# A2     0.0006512284
### SO, BC=4.0 for AD1, 4.5 for A2 and F1, 5 for D5
res0=res
res1=rbind(res[res$BC=="4.5"&res$genome=="A2",],res[res$BC=="4.0"&grepl("AD1",res$genome),],res[res$BC=="4.5"&grepl("F1",res$genome),],res[res$BC=="5.0"&res$genome=="D5",])
res=res1
#res$genome=factor(res$genome,levels=rev(c("A2","F1-At","AD1-At","D5","F1-Dt","AD1-Dt")))
res$genome=factor(res$genome,levels=rev(c("A2","D5","F1-At","F1-Dt","AD1-At","AD1-Dt")))
pdf("Qregulation3/compareSPOs_BC.pdf")
ggplot(res, aes(x = genome, y = sizeP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw()  + labs(title="Length Proportion")+ coord_flip()
ggplot(res, aes(x = genome, y = numberP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw()  + labs(title="Number Proportion")+ coord_flip()
ggplot(res, aes(x = genome, y = size/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + labs(title="Length (Mbp)")+ coord_flip()
ggplot(res, aes(x = genome, y = number/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + labs(title="Number (million)")+ coord_flip()
ggplot(res, aes(x = genome, y = genomeP*100, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + labs(title="Genome Percentage")+ coord_flip()
dev.off()
# scale with genome size
res=merge(res,ss,by="genome")
res$ploidy="Diploid"
res$ploidy[grep("F1",res$genome)]="F1"
res$ploidy[grep("AD1",res$genome)]="AD1"
res$ploidy=factor(res$ploidy,levels=c("Diploid","F1","AD1"))
write.table(res,"Qregulation3/compareSPOs_BC.txt",sep="\t")
pdf("Qregulation3/compareSPOs_XY_BC.pdf")
p1=ggplot(res, aes(x=gsize/10^9, y=number/10^6, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p2=ggplot(res, aes(x=gsize/10^9, y=size/10^7, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p3=ggplot(res, aes(x=gsize/10^9, y=numberP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p4=ggplot(res, aes(x=gsize/10^9, y=sizeP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p5=ggplot(res, aes(x=gsize/10^9, y=genomeP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
grid.arrange( p3 + theme(legend.position = "none",axis.text = element_text(size=8)), p2+theme(legend.position = "none",axis.text = element_text(size=8)), p1+theme(legend.position="none",axis.text = element_text(size=8)), p4+theme(legend.position="none",axis.text = element_text(size=8)), nrow = 2, ncol =2,layout_matrix = rbind(c(1,3), c(2,4)))
p1
p2
p3
p4
p5
dev.off()
## SPO pattern same as MSF, mostly, probably left it

### ACR, combined MSF and SPO, genomic and TE annotation
load("Qregulation3/ACRsC.rdata")
load("refGenomes/TEannotation.rdata")
# "gr.F1"  "gr.A2"  "gr.D5"  "gr.AD1"
gr.F1a=gr.F1[grep("^A",seqnames(gr.F1))]
gr.F1d=gr.F1[grep("^D",seqnames(gr.F1))]
gr.AD1a=gr.AD1[grep("^A",seqnames(gr.AD1))]
gr.AD1d=gr.AD1[grep("^D",seqnames(gr.AD1))]
seqlevels(gns.AD1a)<-gsub("Chr","A",seqlevels(gns.AD1a))
seqlevels(gns.AD1d)<-gsub("Chr","D",seqlevels(gns.AD1d))
#fix fls
fls$ref[fls$genome=="F1a"]="refGenomes/F2020_26.fasta"
fls$ref[fls$genome=="F1d"]="refGenomes/F2020_26.fasta"
fls$chrPath[fls$genome=="F1a"]="mappingF2020/chr.size.A.txt"
fls$chrPath[fls$genome=="F1d"]="mappingF2020/chr.size.D.txt"
fls$gns[fls$genome=="F1a"]="gns.F1a"
fls$gns[fls$genome=="F1d"]="gns.F1d"
# loop
for(g in c("A2","D5","F1a","F1d","AD1a","AD1d")){
    print(g)
    peaks=get(paste0("ACR.",g))
    if(g=="A2"){seqlevels(peaks)<-gsub("Chr","A",seqlevels(peaks))}
    if(g=="D5"){seqlevels(peaks)<-gsub("Chr","D",seqlevels(peaks))}
    tes=get(paste0("gr.",g))
    seqlevels(peaks) <- seqlevelsInUse(peaks)
    seqlevels(tes) <- seqlevelsInUse(tes)
    # TE categorization
    tes$class <- gsub(".*Classification=|;Sequence_ontology.*","",tes$family)

    # annotate ACRs first with TE
    peaksTE=annotateACRs_TE(peaks,tes=tes,upsetplot=TRUE)
    p.te.upset = peaksTE$plot
    
    # then with respect to nearby genes
    peaks=peaksTE$gr
    if(g %in% c("A2","D5")){
        seqlevels(peaks)<-gsub("A|D","Chr",seqlevels(peaks))
    }
    gns<-get(fls$gns[fls$genome==g])
    df = annotateACRs(peaks,gns, distance=2000)
    d_df = df[df$type=="Distal",c(1:3)]
    p_df = df[df$type=="Proximal",c(1:3)]
    g_df = df[df$type=="Genic",c(1:3)]
    write.table(d_df, file="temp/dACRs.bed",sep="\t",row.names=FALSE,col.names=FALSE, quote=FALSE)
    write.table(g_df, file="temp/gACRs.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(p_df, file="temp/pACRs.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    
    # to test whether ACR percentage in TEs is higher or lower than expected, generate null distribution by shuffling peak resgions 100 times
    res=peaksTE$enrich
    gc_permutation =data.frame(TE_family=res$TE_family)
    rownames(gc_permutation) =res$TE_family
    # export gene and 2kb flanking region
    flank = gns
    start(flank) = start(gns) - 2000
    end(flank) = end(gns) + 2000
    start(flank)[start(flank)<1] =1
    export(flank, file("temp/gene_with_2k_flanking.bed"), format ="bed")
    if(g %in% c("A2","D5")){
        seqlevels(tes)<-gsub("A|D","Chr",seqlevels(tes))
        }
    for(pi in 1:100){
        # shuffle peaks with BEDtools
        print(pi)
        cmd<-paste0("bedtools shuffle -i temp/dACRs.bed -excl temp/gene_with_2k_flanking.bed -g ",fls$chrPath[fls$genome==g]," >temp/shuffle.bed; bedtools shuffle -i temp/pACRs.bed -incl temp/gene_with_2k_flanking.bed -g ",fls$chrPath[fls$genome==g]," >>temp/shuffle.bed; bedtools shuffle -i temp/gACRs.bed -incl temp/gene_with_2k_flanking.bed -g ",fls$chrPath[fls$genome==g]," >>temp/shuffle.bed")
        system(cmd)
        #system("bedtools shuffle -i temp/dACRs.bed -excl temp/gene_with_2k_flanking.bed -g mappingD/chr.size.txt >temp/shuffle.bed; bedtools shuffle -i temp/pACRs.bed -incl temp/gene_with_2k_flanking.bed -g mappingD/chr.size.txt >>temp/shuffle.bed; bedtools shuffle -i temp/gACRs.bed -incl temp/gene_with_2k_flanking.bed -g mappingD/chr.size.txt >>temp/shuffle.bed")
        # import peaks and get TE overlap
        pp = import("temp/shuffle.bed", format="bed")
        gc_permutation[,pi] = annotateACRs_TE(peaks=pp,tes=tes,upsetplot=FALSE)$en[res$TE_family,"ACR_perc"]
    }
    pp$type= c(rep("Distal",nrow(d_df)),rep("Proximal",nrow(p_df)),rep("Genic",nrow(g_df)))
    res$null_perc_mean = rowMeans(gc_permutation)
    res$enrichment = log2(res$ACR_perc/res$null_perc_mean)
    empiricalP = function(x,y){sum(x>y)/length(y)}
    res$Probablity = 0
    for(i in 1:nrow(res)){
        res$Probablity[i] = empiricalP(res$ACR_perc[i],gc_permutation[i,])
    }
    res_TE =res
    # plot heatmap
    #labeledHeatmap(xx,xLabels=1:2,yLabels=res$TE_family,colors = brewer.pal(n=10,name="PRGn"))
    res0=res
    res0$TE_family=gsub("_.*","",res0$TE_family)
    res0$enrichment[res$enrichment==-Inf] = min(res0$enrichment[!is.infinite(res0$enrichment)])
    p.te = ggplot(data = res0, aes(x = 1, y = TE_family)) + geom_tile(aes(fill = enrichment),color="white") + scale_fill_gradient2(low="purple", high="darkgreen", guide="colorbar") +theme_bw()
    ## ACRs are depleted from TEs in general

    # categorization: ACR size
    res_Gene = aggregate(width(peaks),by=list(df$type),sum)
    names(res_Gene)=c("type","Bp")
    #  Group.1       x
    #   Distal 1108508
    #    Genic 1330851
    # Proximal 1059869
    
    #### pie plot of ACR numbers by category
    # pie(as.numeric(table(df$type) , labels = c("dACRs","gACRs","pACRs"), border="white", col = c("#999999", "#E69F00", "#56B4E9"))
    Prop = as.data.frame(table(df$type))
    names(Prop)=c("type","number")
    Prop$perc = round(as.numeric(Prop$number)/nrow(df)*100,1)
    mycols=few_pal()(3)
    # Add label position
    Prop <- Prop %>%
      arrange(desc(type)) %>%
      mutate(lab.ypos = cumsum(perc) - 0.5*perc)
    Prop$label = paste0(Prop$type,"\n", Prop$perc,"%")
    # pie & donut
    p.pie = ggplot(Prop, aes(x = "", y = perc, fill = type)) +
      geom_bar(stat = "identity", color = "white") +
      coord_polar("y", start = 0) +
      geom_text(aes(y = lab.ypos, label = label), color = "white",size=6) +
      scale_fill_manual(values = mycols) +
      theme_void()
    p.donut = ggplot(Prop, aes(x = 2, y = perc, fill = type)) +
      geom_bar(stat = "identity", color = "white") +
      coord_polar("y", start = 0) +
      geom_text(aes(y = lab.ypos, label = label), color = "white",size=6) +
      scale_fill_manual(values = mycols) +
      theme_void()+xlim(0.5, 2.5)
    p.donut.smallfont = ggplot(Prop, aes(x = 2, y = perc, fill = type)) +
      geom_bar(stat = "identity", color = "white") +
      coord_polar("y", start = 0) +
      geom_text(aes(y = lab.ypos, label = label), color = "white",size=3) +
      scale_fill_manual(values = mycols) +
      theme_void()+xlim(0.5, 2.5)
    res_Gene=merge(res_Gene,Prop, by="type")

    ### distribution of ACR lengths
    # hist(width(peaks), breaks=100, xlim =c(0,1000), freq=F)
    df$width2 = df$width
    df$width2[df$width>500] = 501
    p.width = ggplot(df, aes(x = width2))+ geom_histogram(aes(color = type, fill = type),alpha=0.4, position = "identity") + scale_fill_manual(values = mycols) +  scale_color_manual(values = mycols) + theme_classic() +theme(axis.text = element_text(size=15), legend.text = element_text(size=15))

    ### density of distance to nearest genes
    #meanweight
    density_peak=function(xx){
        x=log10(xx)
        d=density(x)
        i = which.max(d$y)
        peak = d$x[i]
        #h = hist(x, breaks=500)
        #i = which.max(h$counts)
        #peak = h$mids[i]
        return(10^peak)
    }
    md <- df[df$type!="Genic",] %>%
      group_by(type) %>%
      summarise(d.submit = density_peak(distance2nearest))
    res_Gene = merge(res_Gene, as.data.frame(md),by="type",all.x=T)
    require(scales)
    p.distance = ggplot(df, aes(x = distance2nearest))+ geom_density(aes( color=type, fill = type),alpha=0.4) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=trans_format("log10", math_format(10^.x))) + scale_fill_manual(values = mycols[c(1,3)]) +  scale_color_manual(values = mycols[c(1,3)]) + theme_classic() + geom_vline(aes(xintercept = d.submit, color = type),data = md, linetype = "dashed")+theme(axis.text = element_text(size=15), legend.text = element_text(size=15))
    #geom_label(aes(x=c(1311), y = 1.05, label = c("1.3 Kb"))) +
    #geom_label(aes(x=c(37661), y = 1.05, label = c("37.7 Mb")))

    ### get CG distribution in peaks and control
    write.table(df[,1:7], file="temp/peaks.bed",sep="\t",row.names=FALSE,col.names=FALSE, quote=FALSE)
    system(paste0("bedtools nuc -fi ",fls$ref[fls$genome==g]," -bed temp/peaks.bed >temp/peaksCG.bed"))
    df$gc_content=read.table("temp/peaksCG.bed",sep="\t")$V9
    system(paste0("bedtools nuc -fi ",fls$ref[fls$genome==g]," -bed temp/shuffle.bed >temp/shuffleCG.bed"))
    pp$gc_content=read.table("temp/shuffleCG.bed",sep="\t")$V5
    pp$group ="control"
    df$group ="observed"
    GC = rbind(df[,c("type","gc_content","group")],as.data.frame(pp)[,c("type","gc_content","group")])
    gc_anova = anova(f<-aov(gc_content~type+group,data=GC)) #  ***
    res_GC = rbind(TukeyHSD(f)$type,TukeyHSD(f)$group)
    p.gc=ggplot(GC, aes(x=type, y=gc_content, color=type, fill=type)) + geom_violin() + geom_boxplot(width=0.1, color="grey", alpha=0.2) + theme_bw() + theme(legend.position="none",axis.text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("GC content") + scale_fill_manual(values = mycols) +  scale_color_manual(values = mycols) + facet_grid(. ~ group)

    ### plot CG distribution: peaks vs control
    pdf(paste0("Qregulation3/ACRplots_",g,".pdf"))
    textplot(res_Gene)
    textplot(res_TE)
    textplot(res_GC)
    print(p.pie)
    print(p.donut)
    print(p.distance)
    print(p.width)
    print(p.gc)
    print(p.te.upset)
        
    p.donut.smallfont = ggplot(Prop, aes(x = 2, y = perc, fill = type)) +
      geom_bar(stat = "identity", color = "white") +
      coord_polar("y", start = 0) +
      geom_text(aes(y = lab.ypos, label = label), color = "white",size=2.5, fontface = "bold") +
      scale_fill_manual(values = mycols) + xlim(0.5, 2.5)+
      theme_void()+ theme(legend.position = "none")
    grid.arrange( p.donut.smallfont,
    p.distance+ theme(legend.position = "none",axis.text = element_text(size=8)),
    p.width+ theme(legend.position = "none",axis.text = element_text(size=8)),
    p.gc + theme(legend.position = "none",axis.text = element_text(size=8)),
    p.te + theme(axis.text = element_text(size=8)),
    nrow = 2, ncol =3, layout_matrix = rbind(c(1,2,5), c(3,4,5)))
    dev.off()
    save(df, pp, gc_permutation, res_TE, res_Gene, res_GC, gc_anova, file=paste0("Qregulation3/ACR_",g,".rdata"))
}

# compare ACRs
for(g in c("A2","D5","F1a","F1d","AD1a","AD1d")){
    print(g)
    load(paste0("Qregulation3/ACR_",g,".rdata"))
    # d, g, p
    rr=getPeakSum(df,get(fls$chr[fls$genome==g]))
    rr$genome=g
    if(exists("comp")){comp= rbind(comp,rr)}else{comp=rr}
    # what percentage of each dgp ACRs are located in TE
    tt=getPeakSumTE(df,get(fls$chr[fls$genome==g]))
    tt$genome=g
    if(exists("compTEp")){compTEp= rbind(compTEp,tt)}else{compTEp=tt}
    # collect TE enrichment
    en=res_TE
    en$genome=g
    if(exists("enTE")){enTE= rbind(enTE,en)}else{enTE=en}
}
    
pdf("Qregulation3/TEen.pdf")
res0=enTE
res0$enrichment[res0$enrichment==-Inf] = min(res0$enrichment[!is.infinite(res0$enrichment)])
res0$genome=factor(res0$genome,levels=c("A2","F1a","AD1a","D5","F1d","AD1d"))
ggplot(data = res0, aes(x = genome, y = TE_family)) + geom_tile(aes(fill = enrichment),color="white") + scale_fill_gradient2(low="purple", high="darkgreen", guide="colorbar") +theme_bw()
dev.off()
print(comp)
write.table(comp,sep="\t",file="Qregulation3/compareACRs.txt",row.names=F)
print(compTEp)
write.table(compTEp,sep="\t",file="Qregulation3/compareACRs_TE.txt",row.names=F)

res=read.table("Qregulation3/compareACRs.txt",sep="\t",header=T)
res=merge(res,ss,by="genome")
res$genome=factor(res$genome, levels=rev(c("A2","D5","F1a","F1d","AD1a","AD1d")))
res$ploidy="Diploid"
res$ploidy[grep("F1",res$genome)]="F1"
res$ploidy[grep("AD1",res$genome)]="AD1"
res$ploidy=factor(res$ploidy,levels=c("Diploid","F1","AD1"))
# plot
pdf("Qregulation3/compareACRs.pdf")
ggplot(res, aes(x = genome, y = sizeP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw()  + labs(title="Length Proportion")+ coord_flip()
ggplot(res, aes(x = genome, y = numberP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw()  + labs(title="Number Proportion")+ coord_flip()
ggplot(res, aes(x = genome, y = size/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw()  + labs(title="Length (Mbp)")+ coord_flip()
ggplot(res, aes(x = genome, y = number/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw()  + labs(title="Number (million)")+ coord_flip()
ggplot(res, aes(x = genome, y = genomeP*100, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw()  + labs(title="Genome Percentage")+ coord_flip()
dev.off()
# # scale with genome size, plot genome size at X-axis
pdf("Qregulation3/compareACRs_XY.pdf")
resS=res[res$genome%in%c("A2","D5","F1a","F1d","AD1a","AD1d"),]
p1=ggplot(resS, aes(x=gsize/10^9, y=number/10^6, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p2=ggplot(resS, aes(x=gsize/10^9, y=size/10^6, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p3=ggplot(resS, aes(x=gsize/10^9, y=numberP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p4=ggplot(resS, aes(x=gsize/10^9, y=sizeP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p5=ggplot(resS, aes(x=gsize/10^9, y=genomeP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
grid.arrange( p3 + theme(legend.position = "none",axis.text = element_text(size=8)), p1+theme(legend.position = "none",axis.text = element_text(size=8)), p2+theme(legend.position="none",axis.text = element_text(size=8)), p4+theme(legend.position="none",axis.text = element_text(size=8)), nrow = 2, ncol =2,layout_matrix = rbind(c(1,3), c(2,4)))
p1
p2
p3
p4
p5
dev.off()
# stats: more g and p ACRs in AD1 than A2/F1?, not significant for ACRs
library("agricolae")
#size
summary(a<-aov(lm(size~ploidy,res[res$type=="Distal",])))
summary(a<-aov(lm(size~ploidy,res[res$type=="Genic",])))
summary(a<-aov(lm(size~ploidy,res[res$type=="Proximal",]))) #0.0711 .
TukeyHSD(a)
duncan.test(a,"ploidy")
print(duncan.test(a,"ploidy"))  # AD1=Diploids=F1
#sizeP
summary(a<-aov(lm(sizeP~ploidy,res[res$type=="Distal",])))
summary(a<-aov(lm(sizeP~ploidy,res[res$type=="Genic",])))
summary(a<-aov(lm(sizeP~ploidy,res[res$type=="Proximal",])))
TukeyHSD(a)
duncan.test(a,"ploidy")
print(duncan.test(a,"ploidy"))  # AD1=Diploids=F1
# number
summary(a<-aov(lm(number~ploidy,res[res$type=="Distal",])))
summary(a<-aov(lm(number~ploidy,res[res$type=="Genic",])))
summary(a<-aov(lm(number~ploidy,res[res$type=="Proximal",])))
TukeyHSD(a)
duncan.test(a,"ploidy")
print(duncan.test(a,"ploidy"))  # AD1=Diploids=F1
# numberP
summary(a<-aov(lm(numberP~ploidy,res[res$type=="Distal",])))
summary(a<-aov(lm(numberP~ploidy,res[res$type=="Genic",])))
summary(a<-aov(lm(numberP~ploidy,res[res$type=="Proximal",])))
TukeyHSD(a)
duncan.test(a,"ploidy")
print(duncan.test(a,"ploidy"))  # AD1=Diploids=F1
# genomeP
summary(a<-aov(lm(genomeP~ploidy,res[res$type=="Distal",])))
summary(a<-aov(lm(genomeP~ploidy,res[res$type=="Genic",])))
summary(a<-aov(lm(genomeP~ploidy,res[res$type=="Proximal",])))
TukeyHSD(a)
duncan.test(a,"ploidy")
print(duncan.test(a,"ploidy"))  # AD1=Diploids=F1


###########################
## ACR - fine annotation ##
###########################
# ml r-udunits2
# ml r/3.5.0-py2-ufvuwmm
library(ChIPseeker) # module load r-udunits2
library(GenomicFeatures)
library(ggthemes)
library(gridExtra)
library(rtracklayer)
library(GenomicAlignments)
library(plot.matrix)
library(gplots)
load("Qregulation3/ACRsC.rdata")
# fls
fls=data.frame(genome=c("A2","D5","F1a","F1d","F1","AD1a","AD1d","AD1"),tag1=c("A6Dn","DcD","FcDa","FcDd","FcD","McDa","McDd","McD"),tag2=c("A","D","Fa","Fd","F","Ma","Md","M"),chr=c("chr.size.A2","chr.size.D5","chr.size.A2","chr.size.D5","chr.size.F1","chr.size.AD1a","chr.size.AD1d","chr.size.AD1"), chrPath=c("mappingA_WHU/chr.size.txt","mappingD/chr.size.txt","mappingA_WHU/chr.size.txt","mappingD/chr.size.txt","mappingF2020/chr.size.txt","mappingM_UTX/chr.size.A.txt","mappingM_UTX/chr.size.D.txt", "mappingM_UTX/chr.size.txt"), gns=c("gns.A2","gns.D5","gns.F1","gns.F1","gns.F1","gns.AD1a","gns.AD1d","gns.AD1"),ref=c("refGenomes/A2WHU_13.fasta","refGenomes/Dgenome2_13.fasta","refGenomes/A2WHU_13.fasta","refGenomes/Dgenome2_13.fasta","refGenomes/F2020_26.fasta", "refGenomes/TM1utx_26.fasta","refGenomes/TM1utx_26.fasta","refGenomes/TM1utx_26.fasta"))
# some adjustment
fls$txdb=c("refGenomes/txdb.A2WHU.sqlite","refGenomes/txdb.D5.sqlite","refGenomes/txdb.A2WHU.sqlite","refGenomes/txdb.D5.sqlite","refGenomes/txdb.F2020.sqlite","refGenomes/txdb.AD1utx.sqlite","refGenomes/txdb.AD1utx.sqlite","refGenomes/txdb.AD1utx.sqlite")
fls$gns=c("gns.A2","gns.D5","gns.A2","gns.D5","gns.F1","gns.AD1a","gns.AD1d","gns.AD1")
anL=list()
for(g in c("A2","D5","F1a","F1d","AD1a","AD1d")){
    print(g)
    peaks=get(paste0("ACR.",g))
    if(g %in% c("F1a","F1d")){
        seqlevels(peaks) <- seqlevelsInUse(peaks)
        seqlevels(peaks)=gsub("A|D","Chr",seqlevels(peaks))
    }
    gns=get(fls$gns[fls$genome==g])
    txdb=as.character(fls$txdb[fls$genome==g])
    an= annotateACRs(peaks,gns, distance=2000)
    # focus on g and p
    select=which(an$type!="Distal")
    # Chipseeker annotation, need to tune priority
    ann2 = annotatePeak(peaks[select], TxDb=loadDb(txdb), tssRegion=c(-3000, 3000), verbose=FALSE,genomicAnnotationPriority = c("Promoter","Exon", "Intron", "Downstream", "Intergenic","5UTR", "3UTR"))
    ann2@anno$type=an$type[select]
    ann2@anno$annotation2=gsub("on .*","on",ann2@anno$annotation)
    anL[[g]]=ann2
}
lapply(anL,function(x)length(x@anno))
    
# summary
getPeakInfo = function(peakAnno){
    anno =peakAnno@anno
    catLevel = c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Downstream (<1kb)", "Downstream (1-2kb)", "Downstream (2-3kb)", "Exon", "Intron", "Distal Intergenic" , "5' UTR","3' UTR")
    category = factor(anno$annotation2,levels = catLevel)
    df =data.frame(Feature=catLevel, Number = tapply(width(anno),category,sum), RegionSize = tapply(rep(1,length(category)),category,sum) )
    rownames(df)=NULL
    df[is.na(df)]=0
    df$NumberPerc = df$Number/sum(df$Number)
    df$RegionSizePerc = df$RegionSize/sum(df$RegionSize)
    return(df)
}
dflist2df = function(df.list)
{
    samples = names(df.list)
    df = df.list[[1]]
    df$Sample=samples[1]
    for(i in 2:length(samples)){
        x = df.list[[i]]
        x$Sample = samples[i]
        df = rbind(df, x)
    }
    # levels(df$Feature)=c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Exon", "Intron", "Downstream (<=300)","Distal Intergenic" )
    # df$Sample=factor(df$Sample, levels=c("A2","D5","F1.At","F1.Dt","AD1.At","AD1.Dt"))
    return(df)
}
plotPeakInfo = function(df, y=c("Number", "RegionSize", "NumberPerc", "RegionSizePerc") ,xlab="", ylab="Peak region (bp)", title="")
{
    p <- ggplot(df, aes_string(x = "Sample", fill = "Feature", y = y)) + geom_bar(stat="identity")
    p <- p + ylab(ylab) + xlab(xlab) + ggtitle(title)
    p <- p + coord_flip() + theme_bw()
    p <- p+scale_fill_manual(values=rev(getCols(nlevels(df$Feature))), guide=guide_legend(reverse=TRUE))
    return(p)
    print(p)
}
getCols <- function(n) {
    col <- c("#8dd3c7", "#ffffb3", "#bebada",
    "#fb8072", "#80b1d3", "#fdb462",
    "#b3de69", "#fccde5", "#d9d9d9",
    "#bc80bd", "#ccebc5", "#ffed6f")
    
    col2 <- c("#1f78b4", "#ffff33", "#c2a5cf",
    "#ff7f00", "#810f7c", "#a6cee3",
    "#006d2c", "#4d4d4d", "#8c510a",
    "#d73027", "#78c679", "#7f0000",
    "#41b6c4", "#e7298a", "#54278f")
    
    col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
    "#33a02c", "#fb9a99", "#e31a1c",
    "#fdbf6f", "#ff7f00", "#cab2d6",
    "#6a3d9a", "#ffff99", "#b15928")
    
    ## colorRampPalette(brewer.pal(12, "Set3"))(n)
    colorRampPalette(col3)(n)
}

# plot summary
df = dflist2df(lapply(anL, getPeakInfo))
df=df[df$Number>0,]
df$Feature= factor(df$Feature)
df$Sample=factor(df$Sample, levels=rev(c("A2","D5","F1a","F1d","AD1a","AD1d")))
df$ploidy=gsub("[.].*|a$|d$","",df$Sample)
df$ploidy[grep("A2|D5",df$Sample)]="Diploid"
save(anL, df, file="Qregulation3/ACRpg.rdata")

pdf("Qregulation3/annotACRgp.pdf")
p1<-plotPeakInfo(df, y="RegionSize",ylab ="Peak region (bp)")
p2<-plotPeakInfo(df, y="RegionSizePerc",ylab ="Peak region %")
p3<-plotPeakInfo(df, y="Number",ylab ="Peak number")
p4<-plotPeakInfo(df, y="NumberPerc",ylab ="Peak number %")
grid.arrange(p4+theme(legend.position="none"),p1+theme(legend.position="none"),p3+theme(legend.position="none"),p2+theme(legend.position="none"),nrow=2,ncol=2)
p1;p2;p3;p4
print(plotAnnoBar(anL, title="Footprints Distribution"))
dev.off()

# stats: which region have more ACRs in AD1 than A2/F1?
library("agricolae")
for(f in unique(df$Feature)){
    print(f)
    print(summary(a<-aov(lm(RegionSize~ploidy,df[df$Feature==f,]))))
}
# "Promoter (<=1kb), Downstream (<1kb), Distal Intergenic" with singal ploidy effect
summary(a<-aov(lm(RegionSize~ploidy,df[df$Feature=="Promoter (<=1kb)",])))
TukeyHSD(a)
duncan.test(a,"ploidy")
print(duncan.test(a,"ploidy"))  # AD1>Diploids=F1
#
summary(a<-aov(lm(RegionSize~ploidy,df[df$Feature=="Downstream (<1kb)",])))
TukeyHSD(a)
duncan.test(a,"ploidy")
print(duncan.test(a,"ploidy"))  # AD1>Diploids,AD1=F1,AD1=diploid
# save chipseeker annotation result


q("no")


## get ACRpg sequence
library(rtracklayer)
library(ChIPseeker)
library(AnnotationDbi)
load("Qregulation3/ACRpg.rdata")
ref=c("refGenomes/A2WHU_13.fasta","refGenomes/Dgenome2_13.fasta","refGenomes/A2WHU_13.fasta","refGenomes/Dgenome2_13.fasta","refGenomes/TM1utx_26.fasta","refGenomes/TM1utx_26.fasta")
names(ref)<-names(anL)
for(i in names(anL)){
    anL[[i]]
    gr<-anL[[i]]@anno
    print(table(gr$annotation2,gr$type))
    gr<-gr[gr$annotation2=="Promoter (<=1kb)"]
    #print(gr)
    bedN<-paste0("Qregulation3/",i,"_ACRs_promoter1k.bed")
    print(bedN)
    export(gr, bedN, format ="bed")
    cmd<-paste0("bedtools getfasta -fi ",ref[i]," -bed ",bedN," -fo ",gsub(".bed",".fasta",bedN))
    message(cmd)
    system(cmd)
}

getTSS=function(txdb, include=NULL, exclude.pattern=NULL, be=1000, af=0){
    tss =  genes(txdb)
    if(!is.null(exclude.pattern)) {tss = tss[grep(exclude.pattern,tss$gene_id,invert=T)]}
    if(!is.null(include)) {tss = tss[tss$gene_id %in% include]}
    # +
    end(tss[strand(tss)=="+",])  =start(tss[strand(tss)=="+",])+af
    start(tss[strand(tss)=="+",])  =start(tss[strand(tss)=="+",])-be
    # -
    start(tss[strand(tss)=="-",])=end(tss[strand(tss)=="-",])-af
    end(tss[strand(tss)=="-",])=end(tss[strand(tss)=="-",])+be
    # remove duplicated TSSes ie alternative transcripts
    # this keeps the first instance and removes duplicates
    tss=tss[!duplicated(tss),]
    seqlevels(tss)=unique(as.character(seqnames(tss)))
    return(tss)
}
ref=c("refGenomes/A2WHU_13.fasta","refGenomes/Dgenome2_13.fasta","refGenomes/TM1utx_26.fasta")
names(ref)<-c("A2","D5","AD1")
txdbs=c("refGenomes/txdb.A2WHU.sqlite","refGenomes/txdb.D5.sqlite","refGenomes/txdb.AD1utx.sqlite")
names(txdbs)<-c("A2","D5","AD1")
# system("grep ">" RNAseq/refTranscriptomes/Ghirsutum_527_v2.1.transcript_primaryTranscriptOnly.fa| sed s/' pacid=.*'//g  | sed s/'>'//g >primaryID.txt")
for(i in c("A2","D5")){
    txdb <- loadDb(txdbs[i])
    if(i=="A2"){pms = getTSS(txdb, exclude.pattern="GarUnG")}
    if(i=="D5"){pms = getTSS(txdb, exclude.pattern="Gorai.N")}
    start(pms)[start(pms)<0]=1
    bedN<-paste0("Qregulation3/",i,"_promoter1k.bed")
    export(pms, bedN, format ="bed")
    cmd<-paste0("bedtools getfasta -fi ",ref[i]," -bed ",bedN," -fo ",gsub(".bed",".fasta",bedN))
    message(cmd)
    system(cmd)
}

# bedtools getfasta -fi refGenomes/A2WHU_13.fasta -bed Qregulation3/A2_promoter1k.bed -fo Qregulation3/A2_promoter1k.fasta
# bedtools getfasta -fi refGenomes/Dgenome2_13.fasta -bed Qregulation3/D5_promoter1k.bed -fo Qregulation3/D5_promoter1k.fasta

i="AD1"
txdb <- loadDb(txdbs[i])
pms = getTSS(txdb, exclude.pattern="Gohir.1Z")
pms.A<-pms[grep("A",seqnames(pms))] #36116
pms.D<-pms[grep("D",seqnames(pms))] #38784
table(start(pms.A)<0)
table(start(pms.D)<0)
start(pms.A)[start(pms.A)<0]=1
start(pms.D)[start(pms.D)<0]=1
bedN.A<-paste0("Qregulation3/",i,"a_promoter1k.bed")
bedN.D<-paste0("Qregulation3/",i,"d_promoter1k.bed")
export(pms.A, bedN.A, format ="bed")
export(pms.D, bedN.D, format ="bed")
cmd<-paste0("bedtools getfasta -fi ",ref[i]," -bed ",bedN.A," -fo ",gsub(".bed",".fasta",bedN.A))
message(cmd)
system(cmd)
cmd<-paste0("bedtools getfasta -fi ",ref[i]," -bed ",bedN.D," -fo ",gsub(".bed",".fasta",bedN.D))
message(cmd)
system(cmd)

# bedtools getfasta -fi refGenomes/TM1utx_26.fasta -bed Qregulation3/AD1a_promoter1k.bed -fo Qregulation3/AD1a_promoter1k.fasta
# bedtools getfasta -fi refGenomes/TM1utx_26.fasta -bed Qregulation3/AD1d_promoter1k.bed -fo Qregulation3/AD1d_promoter1k.fasta
 
