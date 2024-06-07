########################################################
################ Regulation question 3 #################
########################################################
# Are ACRs conserved across genomes/subgenomes
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core
module load r/3.5.0-py2-ufvuwmm
R

###########################################################
## ACR annotation and check overlaps between A2D5 and F1 ##
###########################################################
# module load r-udunits2
library(ChIPseeker) # module load r-udunits2
library(GenomicFeatures)
library(ggthemes)
library(gridExtra)
options(stringsAsFactors=FALSE)
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

library(VennDiagram)
ibrary(ChIPpeakAnno)

#### MSF
txdbF1 <- loadDb("refGenomes/txdb.F2020.sqlite")
rm(res)
for(bc in c("4.0","4.5","5.0","5.5","6.0","6.5","7.0"))
{
    print(bc)
    a = import(paste0("isegRes/iseg_v1.3.4_083120/A6Dn_bc",bc,".Fus.bed"),format="bed")
    d = import(paste0("isegRes/iseg_v1.3.4_083120/DcD_bc",bc,".Fus.bed"),format="bed")
    f1 = import(paste0("isegRes/iseg_v1.3.4_083120/FcD_bc",bc,".Fus.bed"),format="bed")
    seqlevels(a) = gsub("^Chr","A",seqlevels(a))
    seqlevels(d) = gsub("^Chr","D",seqlevels(d))
    di=c(a,d)
    
    o=countOverlaps(di,f1) #  the tabulated query overlap hits
    print(length(which(o>0))/length(di))
    o=countOverlaps(f1,di) #  the tabulated query overlap hits
    print(length(which(o>0))/length(di))
    
    h=findOverlaps(di,f1)
    opeaks=di[unique(queryHits(h))]
    an=annotateACRs(opeaks, gns.F1, distance=2000)
    rrA=getPeakSum(an[grep("A",an$seqnames),],chr.size.A2)
    rrA$BC = bc
    rrA$genome="A2o"
    rrD=getPeakSum(an[grep("D",an$seqnames),],chr.size.D5)
    rrD$BC = bc
    rrD$genome="D5o"
    if(exists("res")){res= rbind(res, rrA, rrD)}else{res= rbind(rrA, rrD)}
    
    pdf(paste("Qregulation3/venn.",bc,".pdf"))
    res <- makeVennDiagram(Peaks=list(di,f1), NameOfPeaks=c("diploids", "F1"))
    dev.off()
    print(res)
}


###########################################################
##  check DHS overlaps between A2D5 and F1 ##
###########################################################
# module load r-udunits2
library(ChIPseeker) # module load r-udunits2
library(GenomicFeatures)
library(ggthemes)
library(gridExtra)
options(stringsAsFactors=FALSE)

txdb <- loadDb("refGenomes/txdb.F2020.sqlite")
gns= genes(txdb)
gns.F1= gns[grep("scaffold|Contig",seqnames(gns),invert=T)]# 78962
gns.F1a = gns.F1[grep("A",seqnames(gns.F1)),];
gns.F1d = gns.F1[grep("D",seqnames(gns.F1)),];

### g-p-d categorization
source("FUN.acr.r") # ann

library(VennDiagram)
ibrary(ChIPpeakAnno)

#### DNase-seq by Han et al
a = import("DNase/A2/genrich.narrowPeak",format="narrowPeak")
d = import("DNase/D5/genrich.narrowPeak",format="narrowPeak")
f1 = import("DNase/F1/genrich.narrowPeak",format="narrowPeak")
at=f1[grep("A",seqnames(f1)),]
dt=f1[grep("D",seqnames(f1)),]
seqlevels(a) = gsub("^Chr","A",seqlevels(a))
seqlevels(d) = gsub("^Chr","D",seqlevels(d))
di=c(a,d)

length(a);sum(width(a))
# 28029
# 8,732,967
length(d);sum(width(d))
# 59763
# 34,144,955
length(at);sum(width(at))
# 18519
# 5,458,481
length(dt);sum(width(dt))
# 17948
# 5,697,692

length(f1);sum(width(f1))
# 36467
# 11,156,173
# 8732967
length(d);sum(width(d))
# 59763
# 34144955

    
o=countOverlaps(di,f1) #  the tabulated query overlap hits
length(which(o>0)) # 32168
print(length(which(o>0))/length(di)) #0.3664115
o=countOverlaps(f1,di) #  the tabulated query overlap hits
length(which(o>0))
print(length(which(o>0))/length(f1)) #0.944042????

all<-reduce(c(di,f1))
length(all) # 89628
sum(width(all)) # 44,160,833
median(width(all)) # 333
mean(width(all)) # 492.7125

#### Ous

load("Qregulation3/ACRsC.rdata")->l;l
a = ACR.A2
d = ACR.D5
at=ACR.F1a
dt=ACR.F1d
seqlevels(a) = gsub("^Chr","A",seqlevels(a))
seqlevels(d) = gsub("^Chr","D",seqlevels(d))
seqlevels(at) = gsub("^Chr","A",seqlevels(a))
seqlevels(dt) = gsub("^Chr","D",seqlevels(d))
di=c(a,d)
f1=c(at,dt)

o=countOverlaps(di,f1) #  the tabulated query overlap hits
print(length(which(o>0))/length(di)) #0.1047922
o=countOverlaps(f1,di) #  the tabulated query overlap hits
print(length(which(o>0))/length(f1)) #0.08798702

all<-reduce(c(di,f1))
length(all) # 1015558
sum(width(all)) # 54,400,419
median(width(all)) # 41
mean(width(all)) # 41
