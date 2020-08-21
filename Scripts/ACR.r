########################################################
################ Regulation question 3 #################
########################################################
# Are ACRs conserved across genomes/subgenomes

# module load py-deeptools
# module load r/3.5.0-py2-ufvuwmm
# module load bedtools2
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

### DNase-seq DHS macs2, combine
DHS.A2.c = import("DNase/A2/macs2.combine.bed", format="bed")
DHS.D5.c = import("DNase/D5/macs2.combine.bed", format="bed")
DHS.AD1.c = import("DNase/AD1/macs2.combine.bed", format="bed")
DHSsum.c = list(D5=peakSumByChr(DHS.D5.c,chr.size.D5), A2=peakSumByChr(DHS.A2.c,chr.size.A2), AD1=peakSumByChr(DHS.AD1.c,chr.size.AD1))
# D5 1-2%, A2 0.9-1.6%, AD1 0.4-1%
write.table(DHSsum.c,file="Qregulation3/DHS_combined_sum.txt", sep="\t")

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
peaks0 = import("isegRes/iseg_v1.3.4_062020/DcD_bc6.5.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"] # 68713
bl =GRanges("Chr01", IRanges(23162000, 23820000))  # blacklist region
h=findOverlaps(peaks,bl); length(h) # 1174 in black region
MSF.D5=peaks[-queryHits(h)]# 67539
peaks0 = import("isegRes/iseg_v1.3.4_062020/A6Dn_bc6.5.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"] # 138688
MSF.A2=peaks
peaks0 = import("isegRes/iseg_v1.3.4_062020/FcD_bc6.5.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"] # 204028
MSF.F1=peaks
peaks0 = import("isegRes/iseg_v1.3.4_062020/McD_bc6.5.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"] # 168789
MSF.AD1=peaks
#sum
MSFsum = list(D5=peakSumByChr(MSF.D5,chr.size.D5), A2=peakSumByChr(MSF.A2,chr.size.A2), F1=peakSumByChr(MSF.F1,chr.size.F1),AD1=peakSumByChr(MSF.AD1,chr.size.AD1))
write.table(MSFsum,file="Qregulation3/MSF_sum.txt", sep="\t")

# SPO Fr
peaks0 = import("isegRes/SPO/iseg_v1.3.4_060520/D_frenters_bc5.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]
SPO.D5=peaks # 140268
peaks0 = import("isegRes/SPO/iseg_v1.3.4_060520/A_frenters_bc4.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]
SPO.A2=peaks # 324705
peaks0 = import("isegRes/SPO/iseg_v1.3.4_060520/F_frenters_bc4.5.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]
SPO.F1=peaks # 421763
peaks0 = import("isegRes/SPO/iseg_v1.3.4_060520/M_frenters_bc4.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]
SPO.AD1=peaks # 347908
SPOsum = list(D5=peakSumByChr(SPO.D5,chr.size.D5), A2=peakSumByChr(SPO.A2,chr.size.A2), F1=peakSumByChr(SPO.F1,chr.size.F1),AD1=peakSumByChr(SPO.AD1,chr.size.AD1))
write.table(SPOsum,file="Qregulation3/SPO_sum.txt", sep="\t")

save(list=grep("SPO|MSF|chr[.]size",ls(),value=T), file="Qregulation3/ACRs.rdata")

####################
## ACR annotation ##
####################
# module load r-udunits2
library(ChIPseeker) # module load r-udunits2
library(GenomicFeatures)
library(ggthemes)
library(gridExtra)
load("Qregulation3/ACRs.rdata")->l;l

### import genome annotation
#
txdb <- loadDb("refGenomes/txdb.D5.sqlite")
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
#
txdb <- loadDb("refGenomes/txdb.F2020.sqlite")
gns= genes(txdb)
gns.F1= gns[grep("scaffold|Contig",seqnames(gns),invert=T)]# 78962

ss=data.frame(genome=c("A2","D5","F1-At","F1-Dt","AD1-At","AD1-Dt"),
gsize=c(sum(chr.size.A2$V2),sum(chr.size.D5$V2),sum(chr.size.F1$V2[1:13]),sum(chr.size.F1$V2[14:26]),sum(chr.size.AD1$V2[1:13]),sum(chr.size.AD1$V2[14:26])),
geneN=c(length(gns.A2),length(gns.D5),length(gns.F1[grep("A",seqnames(gns.F1))]),length(gns.F1[grep("D",seqnames(gns.F1))]),length(gns.AD1[grep("A",seqnames(gns.AD1))]),length(gns.AD1[grep("D",seqnames(gns.AD1))]))
)
fl=data.frame(genome=c("A2","D5","F1","AD1"),tag1=c("A6Dn","DcD","FcD","McD"),tag2=c("A","D","F","M"), gns=c("gns.A2","gns.D5","gns.F1","gns.AD1"))
rm(res)

### g-p-d categorization
source("FUN.acr.r") # ann
getPeakSum=function(acr_df,chr_size){
    x=data.frame(table(acr_df$type),size=tapply(acr_df$width,list(acr_df$type),sum))
    names(x)=c("type","number","size")
    x$numberP = x$number/sum(x$number)
    x$sizeP = x$size/sum(x$size)
    x$genomeP = x$size/sum(chr_size$V2)
    return(x)
}



#### MSF

# 4 genomes, 6 BC, annotate ACR, summaize length, number, percentage genome; 126 rows
for(g in 1:4){
    genome=fl$genome[g]
    print(genome)
    for(bc in c("4.0","4.5","5.0","5.5","6.0","6.5","7.0"))
    {
        print(bc)
        f=paste0("isegRes/iseg_v1.3.4_062020/",fl$tag1[g],"_bc",bc,".Fus.bed")
        peaks0 = import(f,format="bed")
        peaks = peaks0[peaks0$itemRgb=="#1414FF"] # 68713
        if(genome=="D5"){
            bl =GRanges("Chr01", IRanges(23162000, 23820000))  # blacklist region
            h=findOverlaps(peaks,bl); length(h) # 1174 in black region
            peaks=peaks[-queryHits(h)]# 67539
        }
        an=annotateACRs(peaks, get(paste0("gns.",genome)), distance=2000)
        if(genome %in% c("A2","D5")){
            rr=getPeakSum(an,get(paste0("chr.size.",genome)))
            rr$BC = bc
            rr$genome=genome
            if(exists("res")){res= rbind(res,rr)}else{res=rr}
        }
        if(genome %in% c("F1","AD1")){
            rrA=getPeakSum(an[grep("A",an$seqnames),],get(paste0("chr.size.",genome))[1:13,])
            rrA$BC = bc
            rrA$genome=paste0(genome,"-At")
            rrD=getPeakSum(an[grep("D",an$seqnames),],get(paste0("chr.size.",genome))[14:26,])
            rrD$BC = bc
            rrD$genome=paste0(genome,"-Dt")
            res= rbind(res,rrA, rrD)
        }
    }
}
write.table(res,sep="\t",file="Qregulation3/compareMSFs.txt")
pdf("Qregulation3/compareMSFs_by_genome.pdf")
mycols=few_pal()(3)
### Does BC stringency affect categorization? NO
ggplot(res, aes(x = BC, y = sizeP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=3) + labs(title="Length Proportion")
ggplot(res, aes(x = BC, y = numberP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=3) + labs(title="Number Proportion")
ggplot(res, aes(x = BC, y = size/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=3) + labs(title="Length (Mbp)")
ggplot(res, aes(x = BC, y = number/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=3) + labs(title="Number (million)")
ggplot(res, aes(x = BC, y = genomeP*100, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ genome,nrow=3) + labs(title="Genome Percentage")
dev.off()
pdf("Qregulation3/compareMSFs_by_BC.pdf")
### Does BC stringency affect categorization patterns across genome? NO
res$genome=factor(res$genome,levels=rev(c("A2","D5","F1-At","F1-Dt","AD1-At","AD1-Dt")))
ggplot(res, aes(x = genome, y = sizeP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Length Proportion")+ coord_flip()
ggplot(res, aes(x = genome, y = numberP, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Number Proportion")+ coord_flip()
ggplot(res, aes(x = genome, y = size/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Length (Mbp)")+ coord_flip()
ggplot(res, aes(x = genome, y = number/10^6, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Number (million)")+ coord_flip()
ggplot(res, aes(x = genome, y = genomeP*100, fill = type))+ geom_bar(stat = "identity", color = "white") + scale_fill_manual(values = mycols) + theme_bw() + facet_wrap(~ BC,nrow=3) + labs(title="Genome Percentage")+ coord_flip()
dev.off()
# focus on BC=6.5
res0=res
res=res[res$BC=="6.5",]
res$genome=factor(res$genome,levels=rev(c("A2","F1-At","AD1-At","D5","F1-Dt","AD1-Dt")))
res$genome=factor(res$genome,levels=rev(c("A2","D5","F1-At","F1-Dt","AD1-At","AD1-Dt")))
pdf("Qregulation3/compareMSFs_BC6.5.pdf")
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
write.table(res,"Qregulation3/compareMSFs_BC6.5.txt",sep="\t")
pdf("Qregulation3/compareMSFs_XY_BC6.5.pdf")
p1=ggplot(res, aes(x=gsize/10^9, y=number/10^6, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p2=ggplot(res, aes(x=gsize/10^9, y=size/10^6, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p3=ggplot(res, aes(x=gsize/10^9, y=numberP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p4=ggplot(res, aes(x=gsize/10^9, y=sizeP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
p5=ggplot(res, aes(x=gsize/10^9, y=genomeP, color=type, shape=ploidy)) + geom_point(aes(size=ploidy)) + scale_color_manual(values = mycols)+ theme_classic() + geom_smooth(aes(group=type),method=lm, se=FALSE, fullrange=TRUE) + scale_size_manual(values=c(3,3,3))+  xlim(0.65, 1.6)
grid.arrange( p3 + theme(legend.position = "none",axis.text = element_text(size=8)), p1+theme(legend.position = "none",axis.text = element_text(size=8)), p2+theme(legend.position="none",axis.text = element_text(size=8)), p4+theme(legend.position="none",axis.text = element_text(size=8)), nrow = 2, ncol =2,layout_matrix = rbind(c(1,3), c(2,4)))
p1
p2
p3
p4
p5
dev.off()

# stats: more g and p ACRs in AD1 than A2/F1?
library("agricolae")
summary(a<-aov(lm(size~ploidy,res[res$type=="Distal",])))
summary(a<-aov(lm(size~ploidy,res[res$type=="Genic",])))
summary(a<-aov(lm(size~ploidy,res[res$type=="Proximal",]))) #P=0.0339
TukeyHSD(a)
duncan.test(a,"ploidy")
print(duncan.test(a,"ploidy"))  # AD1>Diploids>F1

#### SPO
# 4 genomes, 7 BC, annotate ACR, summaize length, number, percentage genome; 126 rows
for(g in 1:4){
    genome=fl$genome[g]
    print(genome)
    for(bc in c("4.0","4.5","5.0","5.5","6.0","6.5"))
    {
        print(bc)
        f=paste0("isegRes/SPO/iseg_v1.3.4_060520/",fl$tag2[g],"_frenters_bc",bc,".Fus.bed")
        peaks0 = import(f,format="bed")
        peaks = peaks0[peaks0$itemRgb=="#1414FF"] # 290381
        if(genome=="D5"){
            bl =GRanges("Chr01", IRanges(23162000, 23820000))  # blacklist region
            h=findOverlaps(peaks,bl);
            if(length(h)>0){ # 1174 in black region
                peaks=peaks[-queryHits(h)]}# 67539
        }
        an=annotateACRs(peaks, get(paste0("gns.",genome)), distance=2000)
        if(genome %in% c("A2","D5")){
            rr=getPeakSum(an,get(paste0("chr.size.",genome)))
            rr$BC = bc
            rr$genome=genome
            if(exists("res")){res= rbind(res,rr)}else{res=rr}
        }
        if(genome %in% c("F1","AD1")){
            rrA=getPeakSum(an[grep("A",an$seqnames),],get(paste0("chr.size.",genome))[1:13,])
            rrA$BC = bc
            rrA$genome=paste0(genome,"-At")
            rrD=getPeakSum(an[grep("D",an$seqnames),],get(paste0("chr.size.",genome))[14:26,])
            rrD$BC = bc
            rrD$genome=paste0(genome,"-Dt")
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
### SO, BC=4.0 for AD1&A2, 4.5 for F1, 5 for D5
res0=res
res1=rbind(res[res$BC=="4.0"&res$genome=="A2",],res[res$BC=="4.0"&grepl("AD1",res$genome),],res[res$BC=="4.5"&grepl("F1",res$genome),],res[res$BC=="5.0"&res$genome=="D5",])
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

###########################
## MSF - fine annotation ##
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
load("Qregulation3/ACRs.rdata")->l;l

### g-p-d categorization
source("FUN.acr.r") # ann
getPeakSum=function(acr_df,chr_size){
    x=data.frame(table(acr_df$type),size=tapply(acr_df$width,list(acr_df$type),sum))
    names(x)=c("type","number","size")
    x$numberP = x$number/sum(x$number)
    x$sizeP = x$size/sum(x$size)
    x$genomeP = x$size/sum(chr_size$V2)
    return(x)
}

### import genome annotation
#
txdb <- loadDb("refGenomes/txdb.D5.sqlite")
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
#
txdb <- loadDb("refGenomes/txdb.F2020.sqlite")
gns= genes(txdb)
gns.F1= gns[grep("scaffold|Contig",seqnames(gns),invert=T)]# 78962

# loop
options(stringsAsFactors=FALSE)
fl=data.frame(genome=c("A2","D5","F1","AD1"),tag1=c("A6Dn","DcD","FcD","McD"),tag2=c("A","D","F","M"), gns=c("gns.A2","gns.D5","gns.F1","gns.AD1"), txdb=c("refGenomes/txdb.A2WHU.sqlite","refGenomes/txdb.D5.sqlite","refGenomes/txdb.F2020.sqlite","refGenomes/txdb.AD1utx.sqlite"))
anL=list()
for(g in 1:2){
    genome=fl$genome[g]
    print(genome)
    gns=get(as.character(fl$gns[g]))
    peaks=get(paste0("MSF.",genome))
    # annotate genic, proximal, distal
    an= annotateACRs(peaks,gns, distance=2000)
    txdb=as.character(fl$txdb[g])
    # focus on g and p
    select=which(an$type!="Distal")
    # Chipseeker annotation, need to tune priority
    # default genomicAnnotationPriority = c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Intergenic")
    #ann0 = annotatePeak(peaks[select], TxDb=loadDb(txdb), tssRegion=c(-2000, 2000), verbose=FALSE)
    #ann1 = annotatePeak(peaks[select], TxDb=loadDb(txdb), tssRegion=c(-2000, 2000), verbose=FALSE,genomicAnnotationPriority = c("Exon", "Intron","Promoter", "Downstream", "Intergenic","5UTR", "3UTR"))
    ann2 = annotatePeak(peaks[select], TxDb=loadDb(txdb), tssRegion=c(-3000, 3000), verbose=FALSE,genomicAnnotationPriority = c("Promoter","Exon", "Intron", "Downstream", "Intergenic","5UTR", "3UTR"))
    #table(an$type[select],gsub("on .*","on",ann2@anno$annotation))
    ann2@anno$type=an$type[select]
    ann2@anno$annotation2=gsub("on .*","on",ann2@anno$annotation)
    #print(plotAnnoBar(ann2))
    #print(plotDistToTSS(ann2, title="Distribution of Pos/Neg footprints relative to TSS"))
    #print(upsetplot(ann2))
    #print(plotAnnoPie(ann2))
    anL[[genome]]=ann2
}
for(g in 3:4){
    genome=fl$genome[g]
    print(genome)
    gns=get(as.character(fl$gns[g]))
    peaks0=get(paste0("MSF.",genome))
    txdb=as.character(fl$txdb[g])
    for(sub in c("A","D")){
        peaks=peaks0[grep(sub,seqnames(peaks0))]
        # annotate genic, proximal, distal
        an= annotateACRs(peaks,gns, distance=2000)
        # focus on g and p
        select=which(an$type!="Distal")
        # Chipseeker annotation, need to tune priority
        ann2 = annotatePeak(peaks[select], TxDb=loadDb(txdb), tssRegion=c(-3000, 3000), verbose=FALSE,genomicAnnotationPriority = c("Promoter","Exon", "Intron", "Downstream", "Intergenic","5UTR", "3UTR"))
        ann2@anno$type=an$type[select]
        ann2@anno$annotation2=gsub("on .*","on",ann2@anno$annotation)
        subgenome=paste0(genome,".",sub,"t")
        anL[[subgenome]]=ann2
    }
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
df$Sample=factor(df$Sample, levels=rev(c("A2","D5","F1.At","F1.Dt","AD1.At","AD1.Dt")))
p1<-plotPeakInfo(df, y="RegionSize",ylab ="Peak region (bp)")
p2<-plotPeakInfo(df, y="RegionSizePerc",ylab ="Peak region %")
p3<-plotPeakInfo(df, y="Number",ylab ="Peak number")
p4<-plotPeakInfo(df, y="NumberPerc",ylab ="Peak number %")
pdf("Qregulation3/annotMSFgp.pdf")
grid.arrange(p4+theme(legend.position="none"),p1+theme(legend.position="none"),p3+theme(legend.position="none"),p2+theme(legend.position="none"),nrow=2,ncol=2)
p1;p2;p3;p4
print(plotAnnoBar(anL, title="Footprints Distribution"))
dev.off()

# stats: which region have more ACRs in AD1 than A2/F1?
library("agricolae")
df$ploidy=gsub("[.].*","",df$Sample)
df$ploidy[grep("A2|D5",df$Sample)]="Diploid"
for(f in unique(df$Feature)){
    print(f)
    print(summary(a<-aov(lm(RegionSize~ploidy,df[df$Feature==f,]))))
}
# "Promoter (<=1kb)" with singal ploidy effect
summary(a<-aov(lm(RegionSize~ploidy,df[df$Feature=="Promoter (<=1kb)",])))
TukeyHSD(a)
duncan.test(a,"ploidy")
print(duncan.test(a,"ploidy"))  # AD1>Diploids=F1

# save chipseeker annotation result
save(anL, df, file="Qregulation3/MSFpg.rdata")


-----book
## Hybridization Effect
load("Qregulation3/MSFpg.rdata")
A2 = anL$A2@anno; seqlevels(A2)=gsub("Chr","A",seqlevels(A2))
At = anL$F1.At@anno
D5 = anL$D5@anno; seqlevels(D5)=gsub("Chr","D",seqlevels(D5))
Dt = anL$F1.Dt@anno
# few intersections, why???
length(h<-findOverlaps(At,A2)) #1149
length(h<-findOverlaps(Dt,D5)) #1170

############################
## Aggregate ACR on genes ##
############################

library(ChIPseeker)
library(gridExtra)
library(dplyr)
source("FUN.aggreNvisual.r")
source("chipseeker.utilities.R")
source("plotTagMatrix.r")

resC

# peak lists
resC$type
grl=list(MSF.bc6.5,SPO.bc5,ATAC.genrich,DHS.macs2.i)
names(grl) = c("MSF.bc6.5","SPO.bc5","ATAC.genrich","DHS.macs2.i")

# gene windows
tss = getPoint(gns,"start",be=3000, af=1000)
tes = getPoint(gns,"end",be=1000, af=3000)

# aggregate
pdf("Qregulation/ACRs_D/ACR_over_genes.pdf")
tagMatrixListS = lapply(grl, getTagMatrix, windows=tss)
tagMatrixListT = lapply(grl, getTagMatrix, windows=tes)
p1=plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS")
p2=plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TES")
ylim= range(c(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range))
grid.arrange(p1+ theme(legend.position = 'bottom') +  scale_y_continuous(limits=ylim), p2+ theme(legend.position = 'bottom')+  scale_y_continuous(limits=ylim), nrow = 1)
plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS",facet="row")
plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TES",facet="row")
dev.off()

### compare ATAC-seq ACRs
grl=list(ATAC.homer.c,ATAC.macs2.c,ATAC.genrich)
names(grl) = c("ATAC.homer.c","ATAC.macs2.c","ATAC.genrich")
# aggregate
pdf("Qregulation/ACRs_D/ACR_over_genes.atac.pdf")
tagMatrixListS = lapply(grl, getTagMatrix, windows=tss)
tagMatrixListT = lapply(grl, getTagMatrix, windows=tes)
p1=plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS")
p2=plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TES")
ylim= range(c(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range))
grid.arrange(p1+ theme(legend.position = 'bottom') +  scale_y_continuous(limits=ylim), p2+ theme(legend.position = 'bottom')+  scale_y_continuous(limits=ylim), nrow = 1)
plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS",facet="row")
plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TES",facet="row")
dev.off()

    


#######################
## Cross Aggregation ##
#######################

# peak lists
resC$type
# ACR homer
peakl=list(MSF.bc6,MSF.bc6.5,SPO.bc4.5,SPO.bc5,ATAC.homer.c,ATAC.homer.i,ATAC.genrich,ATAC.macs2.c,ATAC.macs2.i, DHS.macs2.c,DHS.macs2.i)
names(peakl)=resC$type

# ATAC summit by macs2.i
ATAC.summit <- ATAC.macs2.i
s =start(ATAC.summit)+ATAC.summit$peak-1
end(ATAC.summit)=s
start(ATAC.summit)=s
# DHS summit by mascs2.i
DHS.summit <- DHS.macs2.i
s =start(DHS.summit)+DHS.summit$peak-1
end(DHS.summit)=s
start(DHS.summit)=s
# window list
grl=list(MSF.bc6.5,SPO.bc5,ATAC.summit, DHS.summit)
names(grl) = c("MSF.bc6.5", "SPO.bc5","ATAC.summit", "DHS.summit")

# ACRs over summit/midpoint windows
for(i in 1:length(grl))  # loop
{
    print(flag<-names(grl)[i])
    peaks = get(flag)
    # get center window
    mid=round((start(peaks)+end(peaks))/2)
    cc = peaks
    start(cc)=mid-2000
    end(cc)=mid+2000
    tagMatrixListC = lapply(peakl, getTagMatrix, windows=cc)
    
    pdf(paste0("Qregulation/ACRs_D/ACR_over_",flag,".pdf"))
    print(plotAvgProf.internal(tagMatrixListC, xlim=c(-2000, 2000), origin_label = "Center"))
    print(plotAvgProf.internal(tagMatrixListC, xlim=c(-2000, 2000), origin_label = "Center",facet="row"))
    for(j in 1:length(peakl))
    {
        tagHeatmap(tagMatrixListC[[j]], xlim=c(-2000, 2000), color="red", title=names(tagMatrixListC)[j])
    }
    dev.off()
}

# Plot D BWs over ATAC summit
export(ATAC.summit,"Qregulation/ACRs_D/ATAC.summit.bed")
export(DHS.summit,"Qregulation/ACRs_D/DHS.summit.bed")
load("Qregulation/bws.rdata")
genome="D"
select = (df$genome==genome & df$feature!="NP" & (df$combo==TRUE| grepl("ATAC",df$feature)|grepl("DNase",df$feature)))
use=df[select,]
cfile= "Qregulation/command.aggregateSummit.sh"
cat("module load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
for(bed in c("Qregulation/ACRs_D/ATAC.summit.bed","Qregulation/ACRs_D/DHS.summit.bed")){
    flag=gsub("Qregulation/ACRs_D/|.summit.bed","",bed)
    for (each in use$feature){
        bws=use$bigwig[use$feature==each]
        labels=each
        cmd=paste0("computeMatrix reference-point --referencePoint TSS -S ",bws," -R ",bed," --samplesLabel ",labels," -o summit.gz -b 1500 -a 1500 --skipZeros\nplotHeatmap -m summit.gz --refPointLabel summit -out Qregulation/plot",flag,"Summit.",each,".",genome,".png")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
    }
}
system(paste0("cat ",cfile))
#system(paste0("bash ",cfile))
sbatch(cfile, name="JOB", time = "26:00:00", N=1, n=8)
# Profiles show peak but heatmaps reveal little signal enrichment

# examine nucleosome positioning over ACRs. Null hypothesis is that ATAC is more depleted with nucleosomes
nucleosome = import("Nucleosome/DcH_q20.nucleosome.gff", format="gff")
grl=list(MSF.bc6.5,SPO.bc5,ATAC.summit, DHS.summit)
names(grl) = c("MSF.bc6.5", "SPO.bc5","ATAC.summit", "DHS.summit")
tagMatrixList=list()
pdf("Qregulation/ACRs_D/NP_over_ACRs.pdf")
for(i in 1:length(grl))  # loop
{
    print(flag<-names(grl)[i])
    peaks = get(flag)
    # get center window
    mid=round((start(peaks)+end(peaks))/2)
    cc = peaks
    start(cc)=mid
    end(cc)=mid
    export(cc,paste0("Qregulation/",flag,".mid.bed"),format="bed")
    start(cc)=mid-2000
    end(cc)=mid+2000
    tagMatrix = getTagMatrix(nucleosome, windows=cc)
    tagHeatmap(tagMatrix, xlim=c(-2000, 2000), color="red", title=flag)
    tagMatrixList[[i]] = tagMatrix
}
names(tagMatrixList)=names(grl)
print(plotAvgProf.internal(tagMatrixList, xlim=c(-2000, 2000), origin_label = "Center"))
print(plotAvgProf.internal(tagMatrixList, xlim=c(-2000, 2000), origin_label = "Center",facet="row"))
dev.off()

# plot NP over summit/midpoint windows
load("Qregulation/bws.rdata")
genome="D"
select = (df$genome==genome&df$combo==TRUE)
use = df[select & df$feature=="NP", ]
bws=use$bigwig
bedL =list.files("Qregulation",pattern="mid.bed",full=T)

cfile= "Qregulation/command.aggregateNP.sh"
cat("module load py-deeptools", file=cfile,sep="\n")
cmd=paste0("computeMatrix reference-point --referencePoint center -S ",bws," -R ",paste(bedL,collapse=" ")," -o summit.gz -b 1500 -a 1500 --skipZeros\nplotHeatmap -m summit.gz --refPointLabel summit -out Qregulation/ACRs_D/plotCenter.NP.png\nplotProfile -m summit.gz --refPointLabel summit --perGroup -out Qregulation/ACRs_D/plotCenter.NP.profile.png")
cat(cmd, file =cfile,sep="\n", append=TRUE)
system(paste0("cat ",cfile))
#system(paste0("bash ",cfile))
sbatch(cfile, name="JOB", time = "06:00:00", N=1, n=8)
# Profiles show peak but heatmaps reveal little signal enrichment


## Aggregation on TEs ##
########################


