# This script aims to identify cases to demonstrate mechanistics of duplicated gene regulation by assemble all evidences
# 1. genome level: TE composition
# 2. nucleosome organization
# 3. accessibility
# 4. motif
# 5. expression
# By Guanjing Hu


# module load r/3.5.0-py2-ufvuwmm
options(scipen=999)
library(GenomicRanges)
#library(RColorBrewer)
#library(ggbio)
library(GenomicFeatures)
library(genomation)

getTSS=function(txdb, include=NULL, exclude.pattern=NULL, be=1000, af=1000){
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

setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/")

## focus on OGs
ogQ<-read.table("orthohomoeolog052421.txt", sep="\t", header=FALSE,stringsAsFactor = FALSE)
ogQ$V3<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",ogQ$V3)))
ogQ$V4<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",ogQ$V4)))

## load annotation


## conserved regulatory patterns
con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")

## genic and proximal ACRs
load("Qregulation3/ACRsC.rdata")->k;k
# "ACR.A2"          "ACR.AD1a"        "ACR.AD1d"        "ACR.D5" "ACR.F1a"         "ACR.F1d"         "ACRsum"          "annotateACRs"  "annotateACRs_TE" "chr.size.A2"     "chr.size.AD1"    "chr.size.AD1a"  "chr.size.AD1d"   "chr.size.D5"     "chr.size.F1"     "chr.size.F1a" "chr.size.F1d"

## TE annotation, need to clean up categorization
load("refGenomes/TEannotation.rdata")->k;k #"gr.F1"  "gr.A2"  "gr.D5"  "gr.AD1"
##-- AD1 prep
txdb <- loadDb("refGenomes/txdb.AD1utx.sqlite")
gr=gr.AD1
gr=gr[grep("A|D",seqnames(gr)),]  # exclude scaffold
gr$class=  gsub(".*Classification=|;Sequence_ontology.*","",gr$family)
gr$class=gsub("DNA/Helitron","nonTIR/Helitron",gr$class)
gr$class=gsub("DNA|MITE","TIR",gr$class)
gr$class=gsub("DTC","CACTA",gr$class)
gr$class=gsub("DTM","Mutator",gr$class)
gr$class=gsub("DTH","PIF_Harbinger",gr$class)
gr$class=gsub("DTT","Tc1_Mariner",gr$class)
gr$class=gsub("DTA","hAT",gr$class)
gns<-genes(txdb)
gns=gns[grep("A|D",seqnames(gns)),]  # exclude scaffold
gr.AD1<-gr
gns.AD1<-gns
acr.AD1<-c(ACR.AD1a,ACR.AD1d)
##-- F1 prep
txdb <- loadDb("refGenomes/txdb.F2020.sqlite")
gr=gr.F1
gr=gr[grep("A|D",seqnames(gr)),]  # exclude scaffold
gr$class=  gsub(".*Classification=|;Sequence_ontology.*","",gr$family)
gr$class=gsub("DNA/Helitron","nonTIR/Helitron",gr$class)
gr$class=gsub("DNA|MITE","TIR",gr$class)
gr$class=gsub("DTC","CACTA",gr$class)
gr$class=gsub("DTM","Mutator",gr$class)
gr$class=gsub("DTH","PIF_Harbinger",gr$class)
gr$class=gsub("DTT","Tc1_Mariner",gr$class)
gr$class=gsub("DTA","hAT",gr$class)
gns<-genes(txdb)
gns=gns[grep("A|D",seqnames(gns)),]  # exclude scaffold
gr.F1<-gr
gns.F1<-gns
acr.F1<-c(ACR.F1a,ACR.F1d)
##-- Diploid prep
seqlevels(ACR.A2)<-gsub("Chr","A",seqlevels(ACR.A2))
seqlevels(ACR.D5)<-gsub("Chr","D",seqlevels(ACR.D5))
acr.Diploid <-c(ACR.A2,ACR.D5)

# function to annotate gene relative to feature and get distrance to nearest feature
annotateGENE=function(gns=genes,gr=features, promoterDistance=c(1000,2000,3000)){
    df=data.frame(gns)
    df$width = width(gns)
    # get distance to nearest feature
    d2n = distanceToNearest(gns, gr,ignore.strand=TRUE)
    df$distance2nearest = data.frame(d2n)$distance
    nearest_features = data.frame(gr[subjectHits(d2n)])
    names(nearest_features) = paste0(names(nearest_features),".nearest")
    df=cbind(df,nearest_features)
    # length of features in promoter
    for(d in promoterDistance){
        # make promoter range
        tss=gns
        # +
        end(tss[strand(tss)=="+",])  =start(tss[strand(tss)=="+",])
        start(tss[strand(tss)=="+",])  =start(tss[strand(tss)=="+",])-d
        # -
        start(tss[strand(tss)=="-",])=end(tss[strand(tss)=="-",])
        end(tss[strand(tss)=="-",])=end(tss[strand(tss)=="-",])+d
        # get overlap
        hits <- findOverlaps(query=tss, subject=reduce(gr))
        p <- findOverlapPairs(query=tss, subject=reduce(gr))
        pp<-pintersect(p)
        bp <- width(pp)
        df[paste0("promoter",d)] = 0
        df[queryHits(hits),paste0("promoter",d)] = bp
        }
    return(df)
}

# ACR and TE bp within 1k/2k/3k promoters
resTE.AD1 <- annotateGENE(gns=gns.AD1,gr=gr.AD1)
resACR.AD1 <- annotateGENE(gns=gns.AD1,gr=acr.AD1)
resTE.F1 <- annotateGENE(gns=gns.F1,gr=gr.F1)
resACR.F1 <- annotateGENE(gns=gns.F1,gr=acr.F1)
resACR.Diploid <- annotateGENE(gns=gns.F1,gr=acr.Diploid)
save(resTE.AD1,resTE.F1,resACR.AD1,resACR.F1,resACR.Diploid,file="Qregulation5.rdata")

#### Questions

## Are promoter TE and ACR content negatively correlated? r ranges -0.05 to -0.08
cor.test(resTE.AD1$promoter3000,resACR.AD1$promoter3000)
cor.test(resTE.F1$promoter3000,resACR.F1$promoter3000)
cor.test(resTE.F1$promoter3000,resACR.Diploid$promoter3000)
# -0.03 to -0.07
cor.test(resTE.AD1$promoter1000,resACR.AD1$promoter1000)
cor.test(resTE.F1$promoter1000,resACR.F1$promoter1000)
cor.test(resTE.F1$promoter1000,resACR.Diploid$promoter1000)

##
table(con$Pr)
#     -  Pr<0  Pr=0  Pr>0
#  2993   628 18741   527

select=(con$Pr=="Pr>0")
for(d in c(1000,2000,3000)){
    og=con[select,"V1"]
    at=con[select,"V3"]
    dt=con[select,"V4"]
    a2=paste0(con[select,"V2"],".gene")
    d5=gsub("[.]1$","",con[select,"V5"])
    tedf=data.frame( At=resTE.AD1[match(at,resTE.AD1$gene_id),paste0("promoter",d)], Dt=resTE.AD1[match(dt,resTE.AD1$gene_id),paste0("promoter",d)], A2=resTE.F1[match(a2,resTE.F1$gene_id),paste0("promoter",d)], D5=resTE.F1[match(d5,resTE.F1$gene_id),paste0("promoter",d)] )
    acrdf=data.frame( At=resACR.AD1[match(at,resACR.AD1$gene_id),paste0("promoter",d)], Dt=resACR.AD1[match(dt,resACR.AD1$gene_id),paste0("promoter",d)], A2=resACR.Diploid[match(a2,resACR.Diploid$gene_id),paste0("promoter",d)], D5=resACR.Diploid[match(d5,resACR.Diploid$gene_id),paste0("promoter",d)] )
    

}


---book

# flanking region
    flank = features
    start(flank) = start(features) - distance
    end(flank) = end(features) + distance
    # find gACRs that overlaps with gene and 2kb flanking genes
    df=data.frame(peaks)
    df$width = width(peaks)
    df$type = "Distal"
    hit=findOverlaps(peaks,flank)
    df[queryHits(hit),"type"] = "Proximal"
    hit=findOverlaps(peaks,features)
    df[queryHits(hit),"type"] = "Genic"
    # get distance to nearest genes
    d2n = distanceToNearest(peaks, features)
    df$distance2nearest = data.frame(d2n)$distance
    nearest_genes = data.frame(features[subjectHits(d2n)])
    names(nearest_genes) = paste0(names(nearest_genes),".nearest")
    df=cbind(df,nearest_genes)
    return(df)
}

peakAnnoList <- lapply(grl, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
resAD1<-dflist2df(lapply(peakAnnoList,getPeakInfo))


###############################
## identified ACR comparison ##
###############################
# load p and g ACRs
#setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions")
#load('ACRs_compare/ACRpg.rdata')
# load duplicated gene patterns
con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")

# Func
getGrInfo = function(gr){
    anno =gr
    catLevel = c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Downstream (<1kb)", "Downstream (1-2kb)", "Downstream (2-3kb)", "Exon", "Intron", "Distal Intergenic" , "5' UTR","3' UTR")
    category = factor(anno$annotation2,levels = catLevel)
    df =data.frame(Feature=catLevel, RegionSize = tapply(width(anno),category,sum), Number = tapply(rep(1,length(category)),category,sum) , geneN=tapply(anno$geneId,category,function(x)length(unique(x))))
    rownames(df)=NULL
    df[is.na(df)]=0
    df$NumberPerGene = df$Number/df$geneN
    df$RegionSizePerGene = df$RegionSize/df$geneN
    return(df)
}
dflist2df = function(df.list,label="Sample")
{
    samples = names(df.list)
    df = df.list[[1]]
    df[,label]=samples[1]
    for(i in 2:length(samples)){
        x = df.list[[i]]
        x[,label] = samples[i]
        df = rbind(df, x)
    }
    # levels(df$Feature)=c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Exon", "Intron", "Downstream (<=300)","Distal Intergenic" )
    # df$Sample=factor(df$Sample, levels=c("A2","D5","F1.At","F1.Dt","AD1.At","AD1.Dt"))
    return(df)
}
library(ggplot2)
library(gridExtra)
library(IRanges)

# get promoter ACR size per AD1 gene
# At
gr<-anL$AD1a@anno
gr=gr[gr$annotation2=="Promoter (<=1kb)"]
da<-data.frame(id=gr$geneId,feature=gr$annotation2, size=width(gr))
sa<-aggregate(da$size, list(da$id),sum)
names(sa)=c("id","AD1a")
sa$V1<- con$V1[match(sa$id,con$V3)]
# Dt
gr<-anL$AD1d@anno
gr=gr[gr$annotation2=="Promoter (<=1kb)"]
dd<-data.frame(id=gr$geneId,feature=gr$annotation2, size=width(gr))
sd<-aggregate(dd$size, list(dd$id),sum)
names(sd)=c("id","AD1d")
sd$V1<- con$V1[match(sd$id,con$V4)]
# assemble by OG
cc<-con
ccc<-merge(merge(cc,sa[,c("V1","AD1a")],allx=T),sd[,c("V1","AD1d")],all.x=T)
ccc$AD1a[is.na(ccc$AD1a)]=0
ccc$AD1d[is.na(ccc$AD1d)]=0
# examine ACR sizes, At always higher than Dt regardless of bias
boxplot(ccc$AD1a[ccc$Bp=="Bp<0"], ccc$AD1d[ccc$Bp=="Bp<0"], main="Bp<0")
boxplot(ccc$AD1a[ccc$Bp=="Bp>0"], ccc$AD1d[ccc$Bp=="Bp>0"], main="Bp>0")
boxplot(ccc$AD1a[ccc$Bp=="Bp=0"], ccc$AD1d[ccc$Bp=="Bp=0"], main="Bp=0")
# No clear pattern from scatterplot either
plot(ccc$AD1a,ccc$AD1d,col=factor(ccc$Bp))

# get average size (totalSize/geneNumber)

# AD1
gr<-anL$AD1a@anno
gr$test <-  con$Bp[match(gr$geneId,con$V3)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$AD1d@anno
gr$test <-  con$Bp[match(gr$geneId,con$V4)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pAD1bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1") + theme_bw()
pAD1bar
pAD1box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1") + theme_bw()
pAD1box


# F1
gr<-anL$F1a@anno
gr$test <-  con$B[match(gr$transcriptId,con$V2)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$F1d@anno
gr$test <-  con$B[match(paste0(gr$geneId,".1"),con$V5)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pF1bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1") + theme_bw()
pF1box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1") + theme_bw()
pF1bar
pF1box

# A2 D5
gr<-anL$A2@anno
gr$test <-  con$A[match(gr$transcriptId,con$V2)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$D5@anno
gr$test <-  con$A[match(paste0(gr$geneId,".1"),con$V5)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pA2D5bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("A2 vs D5") + theme_bw()
pA2D5bar
pA2D5box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("A2 vs D5") + theme_bw()
pA2D5box

pdf("Qregulation4/acr2DGE.pdf")
pA2D5bar
pF1bar
pAD1bar
grid.arrange(pA2D5bar+theme(legend.position="none"),pF1bar+theme(legend.position="none"),pAD1bar+theme(legend.position="none"),nrow=3,ncol=1)
pA2D5box
pF1box
pAD1box
grid.arrange(pA2D5box+theme(legend.position="none"),pF1box+theme(legend.position="none"),pAD1box+theme(legend.position="none"),nrow=3,ncol=1)
dev.off()


##########################
## DA result comparison ##
##########################
library(ChIPseeker)
library(ggplot2)
library(gridExtra)
# load duplicated gene patterns
con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")
# Func
getGrInfo = function(gr){
    anno =gr
    catLevel = c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Downstream (<1kb)", "Downstream (1-2kb)", "Downstream (2-3kb)", "Exon", "Intron", "Distal Intergenic" , "5' UTR","3' UTR")
    category = factor(anno$annotation2,levels = catLevel)
    df =data.frame(Feature=catLevel, RegionSize = tapply(width(anno),category,sum), Number = tapply(rep(1,length(category)),category,sum) , geneN=tapply(anno$geneId,category,function(x)length(unique(x))))
    rownames(df)=NULL
    df[is.na(df)]=0
    df$NumberPerGene = df$Number/df$geneN
    df$RegionSizePerGene = df$RegionSize/df$geneN
    return(df)
}
dflist2df = function(df.list,label="Sample")
{
    samples = names(df.list)
    df = df.list[[1]]
    df[,label]=samples[1]
    if(length(samples)>1){
        for(i in 2:length(samples)){
            x = df.list[[i]]
            x[,label] = samples[i]
            df = rbind(df, x)
        }
    }
    # levels(df$Feature)=c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Exon", "Intron", "Downstream (<=300)","Distal Intergenic" )
    # df$Sample=factor(df$Sample, levels=c("A2","D5","F1.At","F1.Dt","AD1.At","AD1.Dt"))
    return(df)
}

# load combined DA
#setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions")
#load('~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/DAanalysis/AD1ref/merged.annotation.rdata')
load("Qregulation_diffACRs/AD1ref/merged.annotation.rdata")
names(anL)
# "resDiA" "resFA"  "resMA"  "resHrA" "resWrA" "resPrA" "resDiD" "resFD"  "resMD"  "resHrD" "resWrD" "resPrD"
#
# AD1
gr<-anL$resMA@anno
gr$test <-  con$Bp[match(gr$geneId,con$V3)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$resMD@anno
gr$test <-  con$Bp[match(gr$geneId,con$V4)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pAD1bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1 ref: AD1") + theme_bw()
pAD1bar
pAD1box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1 ref: AD1") + theme_bw()
pAD1box

# F1
gr<-anL$resFA@anno
gr$test <-  con$B[match(gr$geneId,con$V3)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$resFD@anno
gr$test <-  con$B[match(gr$geneId,con$V4)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pF1bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1 ref: F1") + theme_bw()
pF1box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1 ref: F1") + theme_bw()
pF1bar
pF1box

# A2 D5
gr<-anL$resDiA@anno
gr$test <-  con$A[match(gr$geneId,con$V3)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$resDiD@anno
gr$test <-  con$A[match(gr$geneId,con$V4)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pA2D5bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1 ref: A2 vs D5") + theme_bw()
pA2D5bar
pA2D5box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1 ref: A2 vs D5") + theme_bw()
pA2D5box

pdf("Qregulation4/DAad1ref2DGE.pdf")
pA2D5bar
pF1bar
pAD1bar
grid.arrange(pA2D5bar+theme(legend.position="none"),pF1bar+theme(legend.position="none"),pAD1bar+theme(legend.position="none"),nrow=3,ncol=1)
pA2D5box
pF1box
pAD1box
grid.arrange(pA2D5box+theme(legend.position="none"),pF1box+theme(legend.position="none"),pAD1box+theme(legend.position="none"),nrow=3,ncol=1)
dev.off()

# Hr: DA up or down has no association with expression Hr
gr<-anL$resHrA@anno
gr$test <-  con$Hr[match(gr$geneId,con$V3)]
table(gr$test,gr$direction)
#       down up
# -       1  1
# Hr=0   10 10
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="A"
dfA<-df
gr<-anL$resHrD@anno
gr$test <-  con$Hr[match(gr$geneId,con$V4)]
table(gr$test) # NA
table(gr$direction)
# down   up
#   4    4
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
df$res="resHr"
ddf<-df
# Pr
gr<-anL$resPrA@anno
gr$test <-  con$Pr[match(gr$geneId,con$V3)]
table(gr$test,gr$direction)
#       down mixed up
# -       1     0 14
# Pr<0    1     0  6
# Pr=0   11     0 95
# Pr>0    0     0  3
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="A"
dfA<-df
gr<-anL$resPrD@anno
gr$test <-  con$Pr[match(gr$geneId,con$V4)]
table(gr$test,gr$direction)
#      down up
# -       0  2
# Pr=0    5 16
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
df$res="resPr"
ddf<-rbind(ddf,df)
# Wr
gr<-anL$resWrA@anno
gr$test <-  con$Wr[match(gr$geneId,con$V3)]
table(gr$test,gr$direction)
#      down up
# -       0 11
# Wr<0    0  5
# Wr=0    6 57
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="A"
dfA<-df
gr<-anL$resWrD@anno
gr$test <-  con$Wr[match(gr$geneId,con$V4)]
table(gr$test,gr$direction)
#       down up
# -       1  1
# Wr=0    6 18
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
df$res="resWr"
ddf<-rbind(ddf,df)

df<-ddf[ddf$Feature=="Promoter (<=1kb)"&ddf$direction!="mixed",]
p1<-ggplot(df, aes_string(x = "direction", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~res) + ylab("DA region per gene (bp)") + xlab("") + ggtitle("AD1 ref") + theme_bw()
p2<-ggplot(df, aes_string(x = "direction", fill = "genome", y = "RegionSize")) + geom_bar(position="dodge",stat="identity") + facet_grid(~res) + ylab("DA region (bp)") + xlab("") + ggtitle("AD1 ref") + theme_bw()
pdf("Qregulation4/DAad1ref2p1k.pdf")
grid.arrange(p1+theme(legend.position="none"),p2+theme(legend.position="none"),nrow=3,ncol=1)
p1
p2
dev.off()





# load combined DA
#setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions")
#load('~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/DAanalysis/F1ref/merged.annotation.rdata')
library(ChIPseeker)
load("Qregulation_diffACRs/Fref/merged.annotation.rdata")
names(anL)
# "resDiA" "resFA"  "resMA"  "resHrA" "resWrA" "resPrA" "resDiD" "resFD"  "resMD"  "resHrD" "resWrD" "resPrD"
#
# AD1
gr<-anL$resMA@anno
gr$test <-  con$Bp[match(gr$transcriptId,con$V2)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$resMD@anno
gr$test <-  con$Bp[match(gr$geneId,gsub("[.]1$","",con$V5))]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pAD1bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1 ref: AD1") + theme_bw()
pAD1bar
pAD1box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1 ref: AD1") + theme_bw()
pAD1box

# F1
gr<-anL$resFA@anno
gr$test <-  con$B[match(gr$transcriptId,con$V2)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$resFD@anno
gr$test <-  con$B[match(gr$geneId,gsub("[.]1$","",con$V5))]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pF1bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1 ref: F1") + theme_bw()
pF1box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1 ref: F1") + theme_bw()
pF1bar
pF1box

# A2 D5
gr<-anL$resDiA@anno
gr$test <-  con$A[match(gr$transcriptId,con$V2)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$resDiD@anno
gr$test <-  con$A[match(gr$geneId,gsub("[.]1$","",con$V5))]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pA2D5bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1 ref: A2 vs D5") + theme_bw()
pA2D5bar
pA2D5box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1 ref: A2 vs D5") + theme_bw()
pA2D5box

pdf("Qregulation4/DAf1ref2DGE.pdf")
pA2D5bar
pF1bar
pAD1bar
grid.arrange(pA2D5bar+theme(legend.position="none"),pF1bar+theme(legend.position="none"),pAD1bar+theme(legend.position="none"),nrow=3,ncol=1)
pA2D5box
pF1box
pAD1box
grid.arrange(pA2D5box+theme(legend.position="none"),pF1box+theme(legend.position="none"),pAD1box+theme(legend.position="none"),nrow=3,ncol=1)
dev.off()


# Hr: DA up or down has no association with expression Hr
gr<-anL$resHrA@anno
gr$test <-  con$Hr[match(gr$transcriptId,con$V2)]
table(gr$test) # Hr=0 1
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="A"
dfA<-df
gr<-anL$resHrD@anno
gr$test <-  con$Hr[match(gr$geneId,gsub("[.]1$","",con$V5))]
table(gr$test) # Hr=0 7
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
df$res="resHr"
ddf<-df
# Pr
gr<-anL$resPrA@anno
gr$test <-  con$Pr[match(gr$transcriptId,con$V2)]
table(gr$test,gr$direction)
#      down mixed  up
# -       0     0  15
# Pr<0    0     0   5
# Pr=0    8     0 125
# Pr>0    0     0   3
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="A"
dfA<-df
gr<-anL$resPrD@anno
gr$test <-  con$Pr[match(gr$geneId,gsub("[.]1$","",con$V5))]
table(gr$test,gr$direction)
#      down up
# -       1  5
# Pr<0    1  1
# Pr=0   17 33
# Pr>0    1  0
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
df$res="resPr"
ddf<-rbind(ddf,df)
# Wr
gr<-anL$resWrA@anno
gr$test <-  con$Wr[match(gr$transcriptId,con$V2)]
table(gr$test,gr$direction)
#      down up
# -       1 14
# Wr<0    0  2
# Wr=0   11 96
# Wr>0    0  3
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="A"
dfA<-df
gr<-anL$resWrD@anno
gr$test <-  con$Wr[match(gr$geneId,gsub("[.]1$","",con$V5))]
table(gr$test,gr$direction)
#       down up
# -       0  9
# Wr<0    0  2
# Wr=0   10 28
# Wr>0    0  1
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
df$res="resWr"
ddf<-rbind(ddf,df)

df<-ddf[ddf$Feature=="Promoter (<=1kb)"&ddf$direction!="mixed",]
p1<-ggplot(df, aes_string(x = "direction", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~res) + ylab("DA region per gene (bp)") + xlab("") + ggtitle("F1 ref") + theme_bw()
p2<-ggplot(df, aes_string(x = "direction", fill = "genome", y = "RegionSize")) + geom_bar(position="dodge",stat="identity") + facet_grid(~res) + ylab("DA region (bp)") + xlab("") + ggtitle("F1 ref") + theme_bw()
pdf("Qregulation4/DAf1ref2p1k.pdf")
grid.arrange(p1+theme(legend.position="none"),p2+theme(legend.position="none"),nrow=3,ncol=1)
p1
p2
dev.off()

