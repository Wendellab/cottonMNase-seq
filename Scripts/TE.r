load("TE-ACR.rdata")->ll;ll###############################################################
################ Aggregation and Visulization #################
###############################################################
## This script exploitd aggregated profiles over genomic features (gene body, TSS, TTS, repeats, etc.)
# make comparison across genomes and ploidy

# use speedy
# module load bedtools2 py-deeptools r-udunits2 r/3.5.0-py2-ufvuwmm
# cd /work/LAS/jfw-lab/hugj2006/cottonLeaf
# R
library(GenomicRanges)
#library(RColorBrewer)
#library(ggbio)
library(GenomicFeatures)
library(genomation)
library(ChIPseeker)
options(scipen=999)

g=c("A2WHU","D5","F2020","AD1utx")
ref = data.frame(genome = c("A2", "D5","F1","AD1"), pattern=c("^A","^D","^F","^M"), mapping=c("mappingA_WHU","mappingD","mappingF2020","mappingM_UTX"), txDB =  paste0("refGenomes/txdb.",g,".sqlite"), TE =c("gr.A2","gr.D5","gr.F1","gr.AD1"), subgenome =c(1,1,2,2), nuc = list.files("Nucleosome",pattern="*coverage.bw", full=TRUE) , stringsAsFactors =FALSE)
print(ref)

########################
## Genome composition ##
########################

# gene annotations
list.files("refGenomes",pattern="sqlite")
# [1] "txdb.A2WHU.sqlite"   "txdb.AD1utx.sqlite"  "txdb.AD1utxA.sqlite"  "txdb.AD1utxD.sqlite" "txdb.D5.sqlite"      "txdb.F2020.sqlite"
list.files("refGenomes",pattern="gene.*bed")
# "A2WHU.gene.bed"    "AD1utx.gene.A.bed" "AD1utx.gene.bed" "AD1utx.gene.D.bed" "D5.gene.bed"       "F2020.gene.A.bed" "F2020.gene.bed"    "F2020.gene.D.bed"
list.files("refGenomes",pattern="*og.bed")
# "A2WHU.og.bed"  "AD1utx.og.bed" "D5.og.bed"     "F2020.og.bed"

# load TE annotation
load("refGenomes/TEannotation.rdata")->k;k #"gr.F1"  "gr.A2"  "gr.D5"  "gr.AD1"

#AD1
txdb <- loadDb("refGenomes/txdb.AD1utx.sqlite")

pAD<-c(peakL$AD1a,peakL$AD1d)
type3sub= paste0(pAD$class,".",gsub("[0-9]","",seqnames(pAD)))
grl <- split(pAD,type3sub)
peakAnnoList <- lapply(grl, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
resAD1<-dflist2df(lapply(peakAnnoList,getPeakInfo))
#F1
txdb <- loadDb("refGenomes/txdb.F2020.sqlite")
pF1<-c(peakL$F1a,peakL$F1d)
type3sub= paste0(pF1$class,".",gsub("[0-9]","",seqnames(pF1)))
grl <- split(pF1,type3sub)
peakAnnoList <- lapply(grl, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
resF1<-dflist2df(lapply(peakAnnoList,getPeakInfo))
#
seqlevels(peakL$A2)=gsub("Chr","A",seqlevels(peakL$A2))
seqlevels(peakL$D5)=gsub("Chr","D",seqlevels(peakL$D5))
pDi<-c(peakL$A2,peakL$D5)
type3sub= paste0(pDi$class,".",gsub("[0-9]","",seqnames(pDi)))
grl <- split(pDi,type3sub)
peakAnnoList <- lapply(grl, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
plotAnnoBar(peakAnnoList)
plotDistToTSS(peakAnnoList)
resDiploid<-dflist2df(lapply(peakAnnoList,getPeakInfo))
dev.off()
# combine df
resAD1$ploidy="AD1"
resF1$ploidy ="F1"
resDiploid$ploidy="Diploid"
res<-rbind(resAD1,resF1,resDiploid)
write.table(res,"Qregulation3/ACRsfromTEs.distribution.txt",sep="\t")


## get genomic composition
content = list()
contentByChr = list()
grL=list()
for(i in c(1:2)){
    print(ref[i,])
    # load TE
    gr = get(ref$TE[i])
    gr=gr[grep("scf|tig|ctg|scaffold",seqnames(gr),invert=TRUE),]  # exclude scaffold
    gr$class=  gsub(".*Classification=|;Sequence_ontology.*","",gr$family)# extract classification
    gr$class=gsub("DNA/Helitron","nonTIR/Helitron",gr$class)
    gr$class=gsub("DNA|MITE","TIR",gr$class)
    gr$class=gsub("DTC","CACTA",gr$class)
    gr$class=gsub("DTM","Mutator",gr$class)
    gr$class=gsub("DTH","PIF_Harbinger",gr$class)
    gr$class=gsub("DTT","Tc1_Mariner",gr$class)
    gr$class=gsub("DTA","hAT",gr$class)
    s=aggregate(width(gr),by=list(gr$class),sum)
    ss=aggregate(width(gr),by=list(gr$class, as.character(seqnames(gr))),sum)
    s = rbind(s, c("TE", sum(width(reduce(gr))))) # NOTE the redudnancy!!!!!
    # load gene annot
    txdb <- loadDb(ref$txDB[i])
    gn = genes(txdb);
    gn =gn[grep("Chr",seqnames(gn)),]  # exclude scaffold
    s = rbind(s, c("gene", sum(width(gn))))
    names(s)=c("feature","Bp")
    names(ss)=c("feature","chr","Bp")
    # genome size
    chr.size= read.table(paste0(ref$mapping[i],"/chr.size.txt"),sep="\t")
    if(i==1){chr.size$V1= gsub("Chr","A",chr.size$V1)}
    if(i==2){chr.size$V1= gsub("Chr","D",chr.size$V1)}
    gsize= sum(chr.size$V2)
    s$Perc = as.numeric(s$Bp)/gsize*100
    ss$Perc = as.numeric(ss$Bp)/chr.size$V2[match(ss$chr,chr.size$V1)]*100
    content[[ref$genome[i]]]=s
    contentByChr[[ref$genome[i]]]=ss
    # ChIP-seq annotation
    grL[[ref$genome[i]]] <- gr
}
i=4
print(ref[i,])
gr = get(ref$TE[i])
gr=gr[grep("A|D",seqnames(gr)),]  # exclude scaffold
gr$class=  gsub(".*Classification=|;Sequence_ontology.*","",gr$family)
gr$class=gsub("DNA/Helitron","nonTIR/Helitron",gr$class)
gr$class=gsub("DNA|MITE","TIR",gr$class)
gr$class=gsub("DTC","CACTA",gr$class)
gr$class=gsub("DTM","Mutator",gr$class)
gr$class=gsub("DTH","PIF_Harbinger",gr$class)
gr$class=gsub("DTT","Tc1_Mariner",gr$class)
gr$class=gsub("DTA","hAT",gr$class)
s=aggregate(width(gr),by=list(gr$class),sum)
s = rbind(s, c("TE", sum(width(reduce(gr)))))
# load gene annot
txdb <- loadDb(ref$txDB[i])
gn = genes(txdb);
gn =gn[grep("A|D",seqnames(gn)),]  # exclude scaffold
s = rbind(s, c("gene", sum(width(genes(txdb)))))
names(s)=c("feature","Bp")
ss=aggregate(width(gr),by=list(gr$class, as.character(seqnames(gr))),sum)
names(ss)=c("feature","chr","Bp")
# genome size
chr.size= read.table(paste0(ref$mapping[i],"/chr.size.txt"),sep="\t")
gsize= sum(chr.size$V2)
s$Perc = as.numeric(s$Bp)/gsize*100
ss$Perc = as.numeric(ss$Bp)/chr.size$V2[match(ss$chr,chr.size$V1)]*100
content[[ref$genome[i]]]=s
contentByChr[[ref$genome[i]]]=ss
grL[[ref$genome[i]]] <- gr
## At
gr.sub=gr[grep("^A",seqnames(gr)),]
s=aggregate(width(gr.sub),by=list(gr.sub$class),sum)
names(s)=c("feature","Bp")
s = rbind(s, c("TE", sum(width(reduce(gr.sub)))))
gn.sub=gn[grep("^A",seqnames(gn)),]
s = rbind(s, c("gene", sum(width(gn.sub))))
gsize= sum(chr.size$V2[grep("^A",chr.size$V1)])
s$Perc = as.numeric(s$Bp)/gsize*100
content[["AD1.At"]]=s
# Dt
gr.sub=gr[grep("^D",seqnames(gr)),]
s=aggregate(width(gr.sub),by=list(gr.sub$class),sum)
names(s)=c("feature","Bp")
s = rbind(s, c("TE", sum(width(reduce(gr.sub)))))
gn.sub=gn[grep("^D",seqnames(gn)),]
s = rbind(s, c("gene", sum(width(gn.sub))))
gsize= sum(chr.size$V2[grep("^D",chr.size$V1)])
s$Perc = as.numeric(s$Bp)/gsize*100
content[["AD1.Dt"]]=s
###
content
#make dataframe
write.table(content,"Qregulation3/TEsummary.txt",sep="\t")
save(content,contentByChr, file="Qregulation3/TE.rdata")

#ChIP-seq annotation and plots
getPeakInfo = function(peakAnno){
    anno =peakAnno@anno
    anno$annotation2=gsub("on .*","on",anno$annotation)
    catLevel = c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Downstream (<1kb)", "Downstream (1-2kb)", "Downstream (2-3kb)", "Exon", "Intron", "Distal Intergenic" , "5' UTR","3' UTR")
    category = factor(anno$annotation2,levels = catLevel)
    df =data.frame(Feature=catLevel, RegionSize =tapply(width(anno),category,sum), Number =  tapply(rep(1,length(category)),category,sum) )
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
getCols3 <- function(n) {
    col3 <- c("#a6cee3", "#1f78b4", "#b2df8a",
    "#33a02c", "#fb9a99", "#e31a1c",
    "#fdbf6f", "#ff7f00", "#cab2d6",
    "#6a3d9a", "#ffff99", "#b15928")
    ## colorRampPalette(brewer.pal(12, "Set3"))(n)
    colorRampPalette(col3)(n)
}
getCols2 <- function(n) {
    col <- c("#8dd3c7", "#ffffb3", "#bebada",
    "#fb8072", "#80b1d3", "#fdb462",
    "#b3de69", "#fccde5", "#d9d9d9",
    "#bc80bd", "#ccebc5", "#ffed6f")
    col2 <- c("#1f78b4", "#ffff33", "#c2a5cf",
    "#ff7f00", "#810f7c", "#a6cee3",
    "#006d2c", "#4d4d4d", "#8c510a",
    "#d73027", "#78c679", "#7f0000",
    "#41b6c4", "#e7298a", "#54278f")
    colorRampPalette(col2)(n)
}
#AD1
txdb <- loadDb("refGenomes/txdb.AD1utx.sqlite")
pAD=grL$AD1
type3sub= paste0(pAD$class,".",gsub("[0-9]","",seqnames(pAD)))
grl <- split(pAD,type3sub)
peakAnnoList <- lapply(grl, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
resAD1<-dflist2df(lapply(peakAnnoList,getPeakInfo))
#F1
txdb <- loadDb("refGenomes/txdb.F2020.sqlite")
pF1<-c(grL$A2,grL$D5)
type3sub= paste0(pF1$class,".",gsub("[0-9]","",seqnames(pF1)))
grl <- split(pF1,type3sub)
peakAnnoList <- lapply(grl, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
resF1<-dflist2df(lapply(peakAnnoList,getPeakInfo))
#
resAD1$ploidy="AD1"
resF1$ploidy="Diploid"
res<-rbind(resAD1,resF1)
write.table(res,"Qregulation3/TEs.distribution.txt",sep="\t")
#
res<-read.table("Qregulation3/TEs.distribution.txt",sep="\t",header=T)
res$TE<-factor(gsub(".[A|D]$","",res$Sample), levels=c("ambiguous", "LTR/Copia","LTR/Gypsy",         "LTR/unknown", "nonTIR/Helitron","TIR/CACTA" ,"TIR/hAT","TIR/Mutator", "TIR/PIF_Harbinger","TIR/Tc1_Mariner","nonTE" ))
res$genome<-gsub(".*[.]","",res$Sample)
res$ploidy=factor(res$ploidy,levels=c("Diploid","F1","AD1"))
res$Sample<-paste(res$ploidy,res$genome)
res$Sample = factor(res$Sample, levels=c("AD1 D","AD1 A", "Diploid D", "Diploid A", "F1 D", "F1 A"))
##
with(res,anova(lm(RegionSizePerc~Feature+TE+genome+ploidy))) # Feature ***
T=data.frame(with(res,aggregate(RegionSize,list(Sample),sum)))
rownames(T)=T$Group.1;T
# Group.1          x
# AD1 D         AD1 D  639259338
# AD1 A         AD1 A 1384833057
# Diploid D Diploid D  489721782
# Diploid A Diploid A 1550757187
pdf("Qregulation3/TEs0.pdf")
ggplot(data=res, aes(x=ploidy, y=RegionSize, fill=Feature)) + geom_bar(stat="identity") +scale_fill_manual(values=rev(getCols3(nlevels(factor(res$Feature)))), guide=guide_legend(reverse=TRUE))+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle("TEs")
ggplot(data=res, aes(x=ploidy, y=RegionSize, fill=TE)) + geom_bar(stat="identity") + scale_fill_brewer(palette="Paired") + facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle("TEs")
ggplot(data=res, aes(x=ploidy, y=RegionSize, fill=TE)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired") + facet_grid(cols=vars(Feature),rows= vars(genome)) + theme_classic()+ ggtitle("TEs")
for(fe in unique(res$Feature)){
    use=res[res$Feature==fe,]
    p3<-ggplot(data=use, aes(x=ploidy, y=RegionSize, fill=TE)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle(fe)
    print(p3)
}
dev.off()


# plots
library(reshape2)
library(ggplot2)
library(RColorBrewer)
options(stringAsFactor=F)
df<-rbind(contentByChr$A2, contentByChr$D5, contentByChr$AD1)
df$genome=substring(df$chr,1,1)
df$ploidy=factor(rep(c("Diploid","AD1"),each=nrow(contentByChr$A2)*2),levels=c("Diploid","AD1"))
df$Bp=as.numeric(df$Bp)

pdf("Qregulation3/TEs.pdf")
ggplot(df,aes(x=ploidy, y=Bp,fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Dark2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6)) + ylab ("TEs (bp)") + facet_wrap(~ feature, ncol = 3)
ggplot(df,aes(x=ploidy, y=Perc,fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Dark2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6)) + ylab ("TEs (%)") + facet_wrap(~ feature, ncol = 3)
dev.off()
a= aov(Perc~feature+ploidy+genome, df)
TukeyHSD(a) # sig A>D

source("FUN.aggreNvisual.r") # getPoint
source("plotTagMatrix.r") # plotAvgProf.internal
source("chipseeker.utilities.R")
library(gridExtra)
library(dplyr)

## inspect TE in gene vicinity
for(i in 1:4){
    print(ref[i,])
    # load gene annot
    txdb <- loadDb(ref$txDB[i])
    gn = genes(txdb);
    # load TE
    gr = get(ref$TE[i])
    # exclude scaffold
    gn =gn[grep("scaffold|Contig",seqnames(gn),invert=TRUE),]
    gr=gr[grep("scf|tig|ctg|scaffold",seqnames(gr),invert=TRUE),]
    if(i==1){seqlevels(gn)= gsub("Chr","A",seqlevels(gn))}
    if(i==2){seqlevels(gn)= gsub("Chr","D",seqlevels(gn))}
    # gene windows and peak list
    tss = getPoint(gn,"start",be=3000, af=1000)
    tts = getPoint(gn,"end",be=1000, af=3000)
    gr$class=  gsub(".*Classification=|;Sequence_ontology.*","",gr$family)
    gr$class=gsub("DNA/Helitron","nonTIR/Helitron",gr$class)
    gr$class=gsub("DNA|MITE","TIR",gr$class)
    gr$class=gsub("DTC","CACTA",gr$class)
    gr$class=gsub("DTM","Mutator",gr$class)
    gr$class=gsub("DTH","PIF_Harbinger",gr$class)
    gr$class=gsub("DTT","Tc1_Mariner",gr$class)
    gr$class=gsub("DTA","hAT",gr$class)
    grl <- split(gr, gr$class)
    # aggregate
    pdf(paste0("Qregulation3/TEinGeneVicinity.",ref$genome[i],".pdf"))
    tagMatrixListS = lapply(grl, getTagMatrix, windows=tss)
    tagMatrixListT = lapply(grl, getTagMatrix, windows=tts)
    p1=plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS")
    p2=plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TTS")
    ylim= range(c(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range))
    grid.arrange(p1+ theme(legend.position = 'bottom') +  scale_y_continuous(limits=ylim), p2+ theme(legend.position = 'bottom')+  scale_y_continuous(limits=ylim), nrow = 1)
    plot(p1)
    plot(p2)
    plot(plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS",facet="row"))
    plot(plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TTS",facet="row"))
    if(i>2){
        gr$type3sub= paste0(gr$class,".",gsub("[0-9]","",seqnames(gr)))
        grl.sub <- split(gr, gr$type3sub)
        tagMatrixListS = lapply(grl.sub, getTagMatrix, windows=tss)
        tagMatrixListT = lapply(grl.sub, getTagMatrix, windows=tts)
        # for some reason error for Mutator
        for(type in names(grl)[c(1:4,6:9)]){
            print(type)
            p1= plotAvgProf.internal(tagMatrixListS[paste0(type,c(".A",".D"))], xlim=c(-3000, 1000), origin_label = "TSS", conf=0.95)
            p2= plotAvgProf.internal(tagMatrixListT[paste0(type,c(".A",".D"))], xlim=c(-1000, 3000), origin_label = "TTS", conf=0.95)
            ylim= range(c(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range))
            grid.arrange(p1+ theme(legend.position = 'bottom') +  scale_y_continuous(limits=ylim), p2+ theme(legend.position = 'bottom')+  scale_y_continuous(limits=ylim), nrow = 1)
        }
        # annotate TE relatve to genes, and compare distance to TSS profiles among TE families.
        peakAnnoList <- lapply(grl.sub, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
        plotAnnoBar(peakAnnoList)
        plotDistToTSS(peakAnnoList)
    }
    dev.off()
}

## Perform gene-associated genomic annotate for each type of TE using ChIPseeker, which shows the composition of TE annation in regions off gene body, 3', 5', intergenic, etc.
# --- TODO: espcially for gypsy


