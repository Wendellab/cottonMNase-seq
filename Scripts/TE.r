###############################################################
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
load("refGenomes/TEannotation.rdata")

## get genomic composition
content = list()
contentByChr = list()
for(i in c(1:2)){
    print(ref[i,])
    # load TE
    gr = get(ref$TE[i])
    gr=gr[grep("scf|tig|ctg|scaffold",seqnames(gr),invert=TRUE),]  # exclude scaffold
    gr$class=  gsub(".*Classification=|;Sequence_ontology.*","",gr$family)# extract classification
    s=aggregate(width(gr),by=list(gr$class),sum)
    ss=aggregate(width(gr),by=list(gr$class, as.character(seqnames(gr))),sum)
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
}
i=4
print(ref[i,])
gr = get(ref$TE[i])
gr=gr[grep("A|D",seqnames(gr)),]  # exclude scaffold
gr$class=  gsub(".*Classification=|;Sequence_ontology.*","",gr$family)
s=aggregate(width(gr),by=list(gr$class),sum)
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
## At
gr.sub=gr[grep("^A",seqnames(gr)),]
s=aggregate(width(gr.sub),by=list(gr.sub$class),sum)
names(s)=c("feature","Bp")
gn.sub=gn[grep("^A",seqnames(gn)),]
s = rbind(s, c("gene", sum(width(gn.sub))))
gsize= sum(chr.size$V2[grep("^A",chr.size$V1)])
s$Perc = as.numeric(s$Bp)/gsize*100
content[["AD1.At"]]=s
# Dt
gr.sub=gr[grep("^D",seqnames(gr)),]
s=aggregate(width(gr.sub),by=list(gr.sub$class),sum)
names(s)=c("feature","Bp")
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
# plots
library(reshape2)
library(ggplot2)
library(RColorBrewer)
options(stringAsFactor=F)
df<-rbind(contentByChr$A2, contentByChr$D5, contentByChr$AD1)
df$genome=substring(df$chr,1,1)
df$ploidy=factor(rep(c("Diploid","AD1"),each=nrow(contentByChr$A2)*2),levels=c("Diploid","AD1"))
df$Bp=as.numeric(df$Bp)
df$type=gsub("_.*","",df$feature)
pdf("Qregulation3/TEs.pdf")
ggplot(df,aes(x=ploidy, y=Bp,fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Dark2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6)) + ylab ("TEs (bp)") + facet_wrap(~ type, ncol = 6)
ggplot(df,aes(x=ploidy, y=Perc,fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Dark2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6)) + ylab ("TEs (%)") + facet_wrap(~ type, ncol = 6)
dev.off()
a= aov(Perc~feature+ploidy+genome, df)
TukeyHSD(a) # sig D>A


## inspect TE in gene vicinity
ifor(i in 1:4){
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
    grl <- split(gr, gr$class)
    # aggregate
    pdf(paste0("Qregulation3/TEinGeneVicinity.",ref$genome[i],".pdf"))
    tagMatrixListS = lapply(grl, getTagMatrix, windows=tss)
    tagMatrixListT = lapply(grl, getTagMatrix, windows=tts)
    p1=plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS")
    p2=plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TTS")
    ylim= range(c(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range))
    grid.arrange(p1+ theme(legend.position = 'bottom') +  scale_y_continuous(limits=ylim), p2+ theme(legend.position = 'bottom')+  scale_y_continuous(limits=ylim), nrow = 1)
    plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS",facet="row")
    plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TTS",facet="row")
    if(i>2){
        gr$type3sub= paste0(gr$class,".",gsub("[0-9]","",seqnames(gr)))
        grl.sub <- split(gr, gr$type3sub)
        tagMatrixListS = lapply(grl.sub, getTagMatrix, windows=tss)
        tagMatrixListT = lapply(grl.sub, getTagMatrix, windows=tts)
        for(type in names(grl)){
            p1= plotAvgProf.internal(tagMatrixListS[paste0(type,c(".A",".D"))], xlim=c(-3000, 1000), origin_label = "TSS", conf=0.95)
            p2= plotAvgProf.internal(tagMatrixListT[paste0(type,c(".A",".D"))], xlim=c(-1000, 3000), origin_label = "TTS", conf=0.95)
            ylim= range(c(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range))
            grid.arrange(p1+ theme(legend.position = 'bottom') +  scale_y_continuous(limits=ylim), p2+ theme(legend.position = 'bottom')+  scale_y_continuous(limits=ylim), nrow = 1)
        }
    }
    dev.off()
}

## Perform gene-associated genomic annotate for each type of TE using ChIPseeker, which shows the composition of TE annation in regions off gene body, 3', 5', intergenic, etc.
# --- TODO: espcially for gypsy


