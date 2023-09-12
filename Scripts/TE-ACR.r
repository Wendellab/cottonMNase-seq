########################################################
################ Regulation question 3 #################
########################################################
# What is the origin and regulatory roles of TE-ACRs


##############################
## TE-ACR - fine annotation ##
##############################
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

load("Qregulation3/ACRsC.rdata")->k;k
# "ACR.A2"          "ACR.AD1a"        "ACR.AD1d"        "ACR.D5" "ACR.F1a"         "ACR.F1d"         "ACRsum"          "annotateACRs"  "annotateACRs_TE" "chr.size.A2"     "chr.size.AD1"    "chr.size.AD1a"  "chr.size.AD1d"   "chr.size.D5"     "chr.size.F1"     "chr.size.F1a" "chr.size.F1d"

# load TE annotation
load("refGenomes/TEannotation.rdata")->k;k
# "gr.F1"  "gr.A2"  "gr.D5"  "gr.AD1"


# fls
fls=data.frame(genome=c("A2","D5","F1a","F1d","F1","AD1a","AD1d","AD1"),tag1=c("A6Dn","DcD","FcDa","FcDd","FcD","McDa","McDd","McD"),tag2=c("A","D","Fa","Fd","F","Ma","Md","M"),chr=c("chr.size.A2","chr.size.D5","chr.size.A2","chr.size.D5","chr.size.F1","chr.size.AD1a","chr.size.AD1d","chr.size.AD1"), chrPath=c("mappingA_WHU/chr.size.txt","mappingD/chr.size.txt","mappingA_WHU/chr.size.txt","mappingD/chr.size.txt","mappingF2020/chr.size.txt","mappingM_UTX/chr.size.A.txt","mappingM_UTX/chr.size.D.txt", "mappingM_UTX/chr.size.txt"), gns=c("gns.A2","gns.D5","gns.F1","gns.F1","gns.F1","gns.AD1a","gns.AD1d","gns.AD1"),ref=c("refGenomes/A2WHU_13.fasta","refGenomes/Dgenome2_13.fasta","refGenomes/A2WHU_13.fasta","refGenomes/Dgenome2_13.fasta","refGenomes/F2020_26.fasta", "refGenomes/TM1utx_26.fasta","refGenomes/TM1utx_26.fasta","refGenomes/TM1utx_26.fasta"))
# some adjustment
fls$txdb=c("refGenomes/txdb.A2WHU.sqlite","refGenomes/txdb.D5.sqlite","refGenomes/txdb.A2WHU.sqlite","refGenomes/txdb.D5.sqlite","refGenomes/txdb.F2020.sqlite","refGenomes/txdb.AD1utx.sqlite","refGenomes/txdb.AD1utx.sqlite","refGenomes/txdb.AD1utx.sqlite")
fls$gns=c("gns.A2","gns.D5","gns.A2","gns.D5","gns.F1","gns.AD1a","gns.AD1d","gns.AD1")


getPeakSumTE=function(acr_df,chr_size){
    acr_df$TE=gsub("DNA/Helitron","nonTIR/Helitron",acr_df$TE)
    acr_df$TE=gsub("DNA|MITE","TIR",acr_df$TE)
    acr_df$TE=gsub("/.*","",acr_df$TE)
    acr_df$TE[is.na(acr_df$TE)]="nonTE"
    x=data.frame(table(acr_df$TE),size=tapply(acr_df$width,list(acr_df$TE),sum))
    names(x)=c("type","number","size")
    x$numberP = x$number/sum(x$number)
    x$sizeP = x$size/sum(x$size)
    return(x)
}

# review previous TE-ACR annotation
content = list()
contentByChr = list()
for(g in c("A2","D5","F1a","F1d","AD1a","AD1d"))
{
    print(g)
    load(paste0("Qregulation3/ACR_",g,".rdata"))
    #  "df"             "pp"             "gc_permutation" "res_TE" "res_Gene"       "res_GC"         "gc_anova"
   
   # what percentage of each dgp ACRs are located in TE
    tt=getPeakSumTE(df,get(fls$chr[fls$genome==g]))
    tt$genome=g
    if(exists("compTEp")){compTEp= rbind(compTEp,tt)}else{compTEp=tt}
   
   # collect TE enrichment
    en=res_TE
    en$genome=g
    if(exists("enTE")){enTE= rbind(enTE,en)}else{enTE=en}
    
   # extract TE-ACR
   dfte<-df[!is.na(df$TE),]
   class=gsub("DNA/Helitron","nonTIR/Helitron",dfte$TE)
   class=gsub("DNA|MITE","TIR",class)
   class=gsub("DTC","CACTA",class)
   class=gsub("DTM","Mutator",class)
   class=gsub("DTH","PIF_Harbinger",class)
   class=gsub("DTT","Tc1_Mariner",class)
   class=gsub("DTA","hAT",class)
   class[grepl(";",class)]<-"ambiguous"
   s=aggregate(dfte$width,by=list(class),sum)
   names(s)=c("feature","Bp")
   ss=aggregate(dfte$width,by=list(class, as.character(dfte$seqnames)),sum)
   names(ss)=c("feature","chr","Bp")
   content[[g]]<-s
   contentByChr[[g]]<-ss
   
   grr<-makeGRangesFromDataFrame(dfte, keep.extra.columns=TRUE)
   grr$class=class
   grr$type3sub= paste0(grr$class,".",gsub("[0-9]","",seqnames(grr)))
   
}

save(content,contentByChr, file="Qregulation3/TE-ACR.rdata")

# plots
library(reshape2)
library(ggplot2)
library(RColorBrewer)
options(stringAsFactor=F)
load("Qregulation3/TE-ACR.rdata")
pdf("Qregulation3/TE-ACRs.pdf")
# by  chr
df<-rbind(contentByChr$A2, contentByChr$D5, contentByChr$F1a, contentByChr$F1d, contentByChr$AD1a, contentByChr$AD1d)
df$genome=rep(rep(c("A","D"),each=nrow(contentByChr$A2)),3)
df$ploidy=factor(rep(c("Diploid","F1","AD1"),each=nrow(contentByChr$A2)*2),levels=c("Diploid","F1","AD1"))
df$Bp=as.numeric(df$Bp)
ggplot(df[df$feature!="ambiguous",],aes(x=ploidy, y=Bp,fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Dark2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6)) + ylab ("TEs (bp)") + facet_wrap(~ feature, ncol = 3)
ggplot(df,aes(x=ploidy, y=Bp,fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Dark2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6)) + ylab ("TEs (bp)") + facet_wrap(~ feature, ncol = 3)
# proportion
content=lapply(content,function(x){x$Prop=x$Bp/sum(x$Bp);return(x)})
df<-rbind(content$A2, content$D5, content$F1a, content$F1d, content$AD1a, content$AD1d)
df$genome=rep(rep(c("A","D"),each=nrow(content$A2)),3)
df$ploidy=factor(rep(c("Diploid","F1","AD1"),each=nrow(content$A2)*2),levels=c("Diploid","F1","AD1"))
df$Bp=as.numeric(df$Bp)
ggplot(data=df, aes(x=ploidy, y=Bp, fill=feature)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_wrap(~ genome, nrow = 1) + theme_classic()
ggplot(data=df, aes(x=ploidy, y=Prop, fill=feature)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_wrap(~ genome, nrow = 1) + theme_classic()
dev.off()


peakL<-list()
for(g in c("A2","D5","F1a","F1d","AD1a","AD1d"))
{
    print(g)
    load(paste0("Qregulation3/ACR_",g,".rdata"))
    #  "df"             "pp"             "gc_permutation" "res_TE" "res_Gene"       "res_GC"         "gc_anova"
       
   # classify TE-ACR
   class=df$TE
   class[is.na(df$TE)] = "nonTE"
   class=gsub("DNA/Helitron","nonTIR/Helitron",class)
   class=gsub("DNA|MITE","TIR",class)
   class=gsub("DTC","CACTA",class)
   class=gsub("DTM","Mutator",class)
   class=gsub("DTH","PIF_Harbinger",class)
   class=gsub("DTT","Tc1_Mariner",class)
   class=gsub("DTA","hAT",class)
   class[grepl(";",class)]<-"ambiguous"
   
   grr<-makeGRangesFromDataFrame(df, keep.extra.columns=TRUE)
   grr$class=class
   peakL[[g]]<-grr
}
# custom function
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
pdf("Qregulation3/ACRsfromTEs.distribution.pdf")
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

## make better plots
res<-read.table("Qregulation3/ACRsfromTEs.distribution.txt",sep="\t",header=T)
res$TE<-factor(gsub(".[A|D]$","",res$Sample), levels=c("ambiguous", "LTR/Copia","LTR/Gypsy",         "LTR/unknown", "nonTIR/Helitron","TIR/CACTA" ,"TIR/hAT","TIR/Mutator", "TIR/PIF_Harbinger","TIR/Tc1_Mariner","nonTE" ))
res$genome<-gsub(".*[.]","",res$Sample)
res$ploidy=factor(res$ploidy,levels=c("Diploid","F1","AD1"))
res$Sample<-paste(res$ploidy,res$genome)
res$Sample = factor(res$Sample, levels=c("AD1 D","AD1 A", "Diploid D", "Diploid A", "F1 D", "F1 A"))
##
with(res,anova(lm(RegionSizePerc~Feature+TE+genome+ploidy))) # Feature ***
T=data.frame(with(res,aggregate(RegionSize,list(Sample),sum)))
rownames(T)=T$Group.1;T
# Total ACR size, categorization by TE composition and by distribution
#     Group.1        x
# 1     AD1 A 16243202
# 2     AD1 D 11142951
# 3 Diploid A 16386709
# 4 Diploid D  9183655
# 5      F1 A 19359555
# 6      F1 D 11498496
Perc = res$RegionSize/T[as.character(res$Sample),'x']
# plot Total ACR region and region percentage, colored by genomic feature
p1<-ggplot(data=res, aes(x=ploidy, y=RegionSize, fill=Feature)) + geom_bar(stat="identity") +scale_fill_manual(values=rev(getCols3(nlevels(factor(use$Feature)))), guide=guide_legend(reverse=TRUE))+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle("ACRs")
p2<-ggplot(data=res, aes(x=ploidy, y=Perc, fill=Feature)) + geom_bar(stat="identity") +scale_fill_manual(values=rev(getCols3(nlevels(factor(use$Feature)))), guide=guide_legend(reverse=TRUE))+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle("ACRs")
p3<-ggplot(data=res, aes(x=ploidy, y=RegionSize, fill=Feature)) + geom_bar(stat="identity") +scale_fill_manual(values=rev(getCols3(nlevels(factor(use$Feature)))), guide=guide_legend(reverse=TRUE))+ facet_grid(cols=vars(TE),rows= vars(genome)) + theme_classic()+ ggtitle("ACRs")
p3p<-ggplot(data=res[res$TE!="nonTE",], aes(x=ploidy, y=RegionSize, fill=Feature)) + geom_bar(stat="identity") +scale_fill_manual(values=rev(getCols3(nlevels(factor(use$Feature)))), guide=guide_legend(reverse=TRUE))+ facet_grid(cols=vars(TE),rows= vars(genome)) + theme_classic()+ ggtitle("ACRs")
# plot Total ACR region and region percentage, colored by TE annotation
p4<-ggplot(data=res, aes(x=ploidy, y=RegionSize, fill=TE)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle("ACRs")
p5<-ggplot(data=res, aes(x=ploidy, y=Perc, fill=TE)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle("ACRs")
p6<-ggplot(data=res, aes(x=ploidy, y=RegionSize, fill=TE)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_grid(cols=vars(Feature),rows= vars(genome)) + theme_classic()+ ggtitle("ACRs")
p6p<-ggplot(data=res[res$TE!="nonTE",], aes(x=ploidy, y=RegionSize, fill=TE)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_grid(cols=vars(Feature),rows= vars(genome)) + theme_classic()+ ggtitle("ACRs")
pdf("Qregulation3/ACRs.distribution.pdf")
print(p1)
print(p2)
print(p3)
print(p3p)
print(p4)
print(p5)
print(p6)
print(p6p)
dev.off()


pdf("Qregulation3/ACRs.distribution.perTE.pdf")
for(te in unique(res$TE)){
    use=res[res$TE==te,]
    use$Sample = factor(use$Sample, levels=c("AD1 D","AD1 A", "Diploid D", "Diploid A", "F1 D", "F1 A"))
    use$ploidy=factor(use$ploidy,levels=c("Diploid","F1","AD1"))
    #p1<-ggplot(use, aes_string(x = "Sample", fill = "Feature", y = "RegionSizePerc")) + geom_bar(stat="identity")+ coord_flip() + ylab("RegionSizePerc") + ggtitle(te) + theme_bw()+scale_fill_manual(values=rev(getCols3(nlevels(factor(use$Feature)))), guide=guide_legend(reverse=TRUE))
    p1<-ggplot(data=use, aes(x=ploidy, y=RegionSize, fill=Feature)) + geom_bar(stat="identity") +scale_fill_manual(values=rev(getCols3(nlevels(factor(use$Feature)))), guide=guide_legend(reverse=TRUE))+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle(paste("ACRs overlapped with TEs: ", te))
    p2<-ggplot(data=use, aes(x=ploidy, y=RegionSizePerc, fill=Feature)) + geom_bar(stat="identity") +scale_fill_manual(values=rev(getCols3(nlevels(factor(use$Feature)))), guide=guide_legend(reverse=TRUE))+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle(paste("ACRs overlapped with TEs: ", te))
    print(p1)
    print(p2)
}
dev.off()

pdf("Qregulation3/ACRs.distribution.perFeature.pdf")
for(fe in unique(res$Feature)){
    use=res[res$Feature==fe,]
    use$Sample = factor(use$Sample, levels=c("AD1 D","AD1 A", "Diploid D", "Diploid A", "F1 D", "F1 A"))
    use$ploidy=factor(use$ploidy,levels=c("Diploid","F1","AD1"))
    #p1<-ggplot(use, aes_string(x = "Sample", fill = "TE", y = "RegionSizePerc")) + geom_bar(stat="identity")+ coord_flip() + ylab("RegionSizePerc") + ggtitle(fe) + theme_bw()+ scale_fill_manual(values=rev(getCols2(nlevels(factor(use$TE)))), guide=guide_legend(reverse=TRUE))
    #p2<-ggplot(use, aes_string(x = "Sample", fill = "TE", y = "RegionSize")) + geom_bar(stat="identity")+ coord_flip() + ylab("RegionSize") + ggtitle(fe) + theme_bw()+ scale_fill_manual(values=rev(getCols2(nlevels(factor(use$TE)))), guide=guide_legend(reverse=TRUE))
    use0=use[use$TE!="nonTE",]
    sumF<-data.frame(with(use0,aggregate(RegionSize,list(Sample),sum)))
    names(sumF)<-c("Sample","sumF")
    use0=merge(use0,sumF)
    use0$Prop<-use0$RegionSize/use0$sumF
    p3<-ggplot(data=use, aes(x=ploidy, y=RegionSize, fill=TE)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle(fe)
    print(p3)
    p4<-ggplot(data=use0, aes(x=ploidy, y=RegionSize, fill=TE)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle(fe)
    print(p4)
    p5<-ggplot(data=use0, aes(x=ploidy, y=Prop, fill=TE)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle(fe)
    print(p5)
}
dev.off()

pdf("figureUse2.pdf")
res$TE<-factor(res$TE, levels=c( "LTR/Copia","LTR/Gypsy", "LTR/unknown", "nonTIR/Helitron","TIR/CACTA" ,"TIR/hAT","TIR/Mutator", "TIR/PIF_Harbinger","TIR/Tc1_Mariner","ambiguous","nonTE" ))
print(p4+theme(text = element_text(size = 15)))
p66<-ggplot(data=res, aes(x=ploidy, y=RegionSize, fill=TE)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_wrap(Feature~1, nrow = 3) + theme_classic()+ ggtitle("ACRs") + theme(legend.position="none")
print(p66)
dev.off()

pdf("figureUse3.pdf")
for(fe in unique(res$Feature)){
    use=res[res$Feature==fe,]
    use$Sample = factor(use$Sample, levels=c("AD1 D","AD1 A", "Diploid D", "Diploid A", "F1 D", "F1 A"))
    use$ploidy=factor(use$ploidy,levels=c("Diploid","F1","AD1"))
    #p1<-ggplot(use, aes_string(x = "Sample", fill = "TE", y = "RegionSizePerc")) + geom_bar(stat="identity")+ coord_flip() + ylab("RegionSizePerc") + ggtitle(fe) + theme_bw()+ scale_fill_manual(values=rev(getCols2(nlevels(factor(use$TE)))), guide=guide_legend(reverse=TRUE))
    #p2<-ggplot(use, aes_string(x = "Sample", fill = "TE", y = "RegionSize")) + geom_bar(stat="identity")+ coord_flip() + ylab("RegionSize") + ggtitle(fe) + theme_bw()+ scale_fill_manual(values=rev(getCols2(nlevels(factor(use$TE)))), guide=guide_legend(reverse=TRUE))
    use0=use[use$TE!="nonTE",]
    sumF<-data.frame(with(use0,aggregate(RegionSize,list(Sample),sum)))
    names(sumF)<-c("Sample","sumF")
    use0=merge(use0,sumF)
    use0$Prop<-use0$RegionSize/use0$sumF
    p3<-ggplot(data=use, aes(x=ploidy, y=RegionSize, fill=TE)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle(fe)
    print(p3+theme(text = element_text(size = 15)))
    p4<-ggplot(data=use0, aes(x=ploidy, y=RegionSize, fill=TE)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle(fe)
    print(p4+theme(text = element_text(size = 15)))
    p5<-ggplot(data=use0, aes(x=ploidy, y=Prop, fill=TE)) + geom_bar(stat="identity") +scale_fill_brewer(palette="Paired")+ facet_wrap(~ genome, nrow = 1) + theme_classic()+ ggtitle(fe)
    print(p5+theme(text = element_text(size = 15)))
}
dev.off()

# correlate ACR-TE with overall TEs
load('TE.rdata')->l;l
all<-content
load("TE-ACR.rdata")->ll;ll
cor.test(as.numeric(all$D5$Bp[1:9]),content$F1d$Bp[2:10])
cor.test(as.numeric(all$A2$Bp[1:9]),content$F1a$Bp[2:10])
cor.test(as.numeric(all$D5$Bp[1:9]),content$D5$Bp[2:10])
cor.test(as.numeric(all$A2$Bp[1:9]),content$A2$Bp[2:10])
cor.test(as.numeric(all$AD1.At$Bp[1:9]),content$AD1a$Bp[2:10])
cor.test(as.numeric(all$AD1.Dt$Bp[1:9]),content$AD1d$Bp[2:10])


library(reshape2);
aggregate(res$RegionSize,list(res$Sample,res$TE),sum)

####################################
## promoter TE - ACR - Expression ##
####################################
# for each OG, ask whether expression pattern is associated with TE presentation (amount with 1-2-3kb, distance to nearest TEs)



#################################
## DA analysis relative to TEs ##
#################################
# ml r-udunits2
# ml r/3.5.0-py2-ufvuwmm
options(scipen=999)
library(ChIPseeker) # module load r-udunits2
library(GenomicFeatures)
library(ggthemes)
library(gridExtra)
library(rtracklayer)
library(GenomicAlignments)
library(plot.matrix)
library(gplots)
library(ggplot2)


source("FUN.acr.r")
# load TE annotation
load("refGenomes/TEannotation.rdata")->k;k
# "gr.F1"  "gr.A2"  "gr.D5"  "gr.AD1"
addTEclass<-function(gr){
    gr$class<- gsub(".*Classification=|;Sequence_ontology.*","",gr$family)
    gr$class=gsub("DNA/Helitron","nonTIR/Helitron",gr$class)
    gr$class=gsub("DNA|MITE","TIR",gr$class)
    gr$class=gsub("DTC","CACTA",gr$class)
    gr$class=gsub("DTM","Mutator",gr$class)
    gr$class=gsub("DTH","PIF_Harbinger",gr$class)
    gr$class=gsub("DTT","Tc1_Mariner",gr$class)
    gr$class=gsub("DTA","hAT",gr$class)
    gr$class[grepl(";",gr$class)]<-"ambiguous"
    return(gr)
}

# load combined DA based on F1 reference, which resulted into more changes than using AD1 ref
#setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions")
#load('~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/DAanalysis/F1ref/merged.annotation.rdata')
library(ChIPseeker)
load("Qregulation_diffACRs/Fref/merged.annotation.rdata")
names(anL)
# "resDiA" "resFA"  "resMA"  "resHrA" "resWrA" "resPrA" "resDiD" "resFD"  "resMD"  "resHrD" "resWrD" "resPrD"
# plotAnnoBar(anL)

pdf("Qregulation4/DAf1ref.width.pdf")
# line plot of DA region width
cl<-brewer.pal(length(anL),"Paired")
medianWidth<-c()
for(i in 1:length(anL)){
    print(names(anL)[i])
    gr<-anL[[i]]@anno
    print(quantile(width(gr)))
    medianWidth[i] = median(width(gr))
    if(i==1){
        plot(density(width(gr)),ylim=c(0,0.00013),col=cl[i])
    }else{
        lines(density(width(gr)),col=cl[i])
    }
}
legend("topright",legend=paste(names(anL),medianWidth),col=cl, lty=1)

# boxplot of DA width by category
for(i in 1:length(anL)){
    print(names(anL)[i])
    gr<-anL[[i]]@anno
    p<-ggplot(as.data.frame(gr), aes(x=annotation2, y=width(gr), fill=direction)) +  geom_violin() + ylab("DA region width (bp)") + xlab("") + ggtitle(names(anL)[i]) + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6))
    print(p)
}
dev.off()


# genic and TE annotation
tes<-addTEclass(gr.F1)
txdb <- loadDb("refGenomes/txdb.F2020.sqlite")
gns= genes(txdb)
gns= gns[grep("scaffold|Contig",seqnames(gns),invert=T)]# 78962
anL2<-list()
for(i in 1:length(anL)){
    print(names(anL)[i])
    gr<-anL[[i]]@anno
    # annotate ACRs first with TE
    peaksTE=annotateACRs_TE(peaks=gr,tes=tes,upsetplot=TRUE)
    p.te.upset = peaksTE$plot
    # then with respect to nearby genes
    peaks=peaksTE$gr
    df = annotateACRs(peaks,gns, distance=2000)
    anL2[[i]]=df
}
names(anL2)=names(anL)
save(anL2, file="Qregulation4/Fref.DA.rdata")

# modified above scripts to make "Qregulation4/AD1ref.DA.rdata"

## plot DA regions relative genes and TEs
