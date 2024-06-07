########################################################
################ Regulation question 3 #################
########################################################
# What is the origin and regulatory roles of TE-ACRs


##############################
## TE-ACR - fine annotation ##
##############################
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core
module load r/3.5.0-py2-ufvuwmm
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf
R


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

#gN<-c(length(gns.A2),length(gns.D5),length(gns.A2),length(gns.D5),length(gns.AD1a),length(gns.AD1d))
gN<-c(41739,37223,41739,37223,36118,38784)
## make better plots
res<-read.table("Qregulation3/ACRsfromTEs.distribution.txt",sep="\t",header=T)
res$TE<-factor(gsub(".[A|D]$","",res$Sample), levels=c("ambiguous", "LTR/Copia","LTR/Gypsy",         "LTR/unknown", "nonTIR/Helitron","TIR/CACTA" ,"TIR/hAT","TIR/Mutator", "TIR/PIF_Harbinger","TIR/Tc1_Mariner","nonTE" ))
res$genome<-gsub(".*[.]","",res$Sample)
res$ploidy=factor(res$ploidy,levels=c("Diploid","F1","AD1"))
res$Sample<-paste(res$ploidy,res$genome)
res$Sample = factor(res$Sample, levels=c("AD1 D","AD1 A", "Diploid D", "Diploid A", "F1 D", "F1 A"))
##
with(res,anova(lm(RegionSizePerc~Feature+TE+genome+ploidy))) # Feature ***
T=data.frame(with(res[res$Feature=="Promoter (<=1kb)",],aggregate(RegionSize,list(Sample),sum)))
rownames(T)=T$Group.1;T
T$perc<-T$x/gN
# Total ACR size, categorization by TE composition and by distribution
#     Group.1        x
# 1     AD1 A 16243202
# 2     AD1 D 11142951
# 3 Diploid A 16386709
# 4 Diploid D  9183655
# 5      F1 A 19359555
# 6      F1 D 11498496

