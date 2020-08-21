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

ref = data.frame(genome = c("A2", "D5","F1","AD1"), pattern=c("^A","^D","^F","^M"), mapping=c("mappingA_new","mappingD","mappingF","mappingMnew"), txDB =  list.files("refGenomes",pattern="txdb",full.names =T), TE =c("gr.A2","gr.D5","gr.F1","gr.AD1"), subgenome =c(1,1,2,2), nuc = list.files("nucleosome",pattern="*coverage.bw", full=TRUE) , stringsAsFactors =FALSE)
print(ref)

########################
## Genome composition ##
########################

# gene annotations
list.files("refGenomes",pattern="sqlite")
# [1] "txdb.A2du.sqlite"     "txdb.D5.sqlite"       "txdb.F1.sqlite"   "txdb.TM1saski.sqlite"
list.files("refGenomes",pattern="gene.*bed")
# "A2.gene.bed"    "AD1.gene.A.bed" "AD1.gene.D.bed" "D5.gene.bed" "F1.gene.A.bed"  "F1.gene.D.bed"
list.files("refGenomes",pattern="*og.bed")
# "A2.og.bed"  "AD1.og.bed" "D5.og.bed"  "F1.og.bed"

# load TE annotation
load("refGenomes/TEannotation.rdata")

## get genomic composition
content = list()
for(i in c(1:2)){
    print(ref[i,])
    # load TE
    gr = get(ref$TE[i])
    gr=gr[grep("Chr",seqnames(gr)),]  # exclude scaffold
    s=aggregate(width(gr),by=list(gr$type),sum)
    # load gene annot
    txdb <- loadDb(ref$txDB[i])
    gn = genes(txdb);
    gn =gn[grep("Chr",seqnames(gn)),]  # exclude scaffold
    s = rbind(s, c("gene", sum(width(gn))))
    names(s)=c("feature","Bp")
    # genome size
    chr.size= read.table(paste0(ref$mapping[i],"/chr.size.txt"),sep="\t")
    gsize= sum(chr.size$V2)
    s$Perc = as.numeric(s$Bp)/gsize*100
    content[[ref$genome[i]]]=s
}
i=4
print(ref[i,])
gr = get(ref$TE[i])
gr=gr[grep("scf",seqnames(gr),invert=TRUE),]  # exclude scaffold
s=aggregate(width(gr),by=list(gr$type),sum)
# load gene annot
txdb <- loadDb(ref$txDB[i])
gn = genes(txdb);
gn =gn[grep("scaffold",seqnames(gn),invert=TRUE),]  # exclude scaffold
s = rbind(s, c("gene", sum(width(genes(txdb)))))
names(s)=c("feature","Bp")
# genome size
chr.size= read.table(paste0(ref$mapping[i],"/chr.size.txt"),sep="\t")
gsize= sum(chr.size$V2)
s$Perc = as.numeric(s$Bp)/gsize*100
content[[ref$genome[i]]]=s
## At
gr.sub=gr[grep("^A",seqnames(gr)),]
s=aggregate(width(gr.sub),by=list(gr.sub$type),sum)
names(s)=c("feature","Bp")
gn.sub=gn[grep("^A",seqnames(gn)),]
s = rbind(s, c("gene", sum(width(gn.sub))))
gsize= sum(chr.size$V2[grep("^A",chr.size$V1)])
s$Perc = as.numeric(s$Bp)/gsize*100
content[["AD1.At"]]=s
# Dt
gr.sub=gr[grep("^D",seqnames(gr)),]
s=aggregate(width(gr.sub),by=list(gr.sub$type),sum)
names(s)=c("feature","Bp")
gn.sub=gn[grep("^D",seqnames(gn)),]
s = rbind(s, c("gene", sum(width(gn.sub))))
gsize= sum(chr.size$V2[grep("^D",chr.size$V1)])
s$Perc = as.numeric(s$Bp)/gsize*100
content[["AD1.Dt"]]=s
###
content
#make dataframe
write.table(content,"TEsummary.txt",sep="\t")
save(content, file="aggreNvisual.rdata")

## inspect TE in gene vicinity
source("aggreNvisual.FUN.r")
source("/work/LAS/jfw-lab/hugj2006/cottonLeaf/chipseeker.utilities.R")
source("/work/LAS/jfw-lab/hugj2006/cottonLeaf/plotTagMatrix.r")
library(gridExtra)
library(dplyr)
for(i in 1:4){
    print(ref[i,])
    # load gene annot
    txdb <- loadDb(ref$txDB[i])
    gn = genes(txdb);
    # load TE
    gr = get(ref$TE[i])
    # exclude scaffold
    if(i>2){
        gn =gn[grep("scaffold",seqnames(gn),invert=TRUE),]
        gr=gr[grep("scf|tig",seqnames(gr),invert=TRUE),]  # exclude scaffold
    }else{
        gn =gn[grep("Chr",seqnames(gn)),]
        gr=gr[grep("Chr",seqnames(gr)),]  # exclude scaffold
    }
    # gene windows and peak list
    tss = getPoint(gn,"start",be=3000, af=1000)
    tts = getPoint(gn,"end",be=1000, af=3000)
    grl <- split(gr, gr$type3)
    # aggregate
    pdf(paste0("aggregationPlots/TEinGeneVicinity.",ref$genome[i],".pdf"))
    tagMatrixListS = lapply(grl, getTagMatrix, windows=tss)
    tagMatrixListT = lapply(grl, getTagMatrix, windows=tts)
    p1=plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS")
    p2=plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TTS")
    ylim= range(c(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range))
    grid.arrange(p1+ theme(legend.position = 'bottom') +  scale_y_continuous(limits=ylim), p2+ theme(legend.position = 'bottom')+  scale_y_continuous(limits=ylim), nrow = 1)
    plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS",facet="row")
    plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TTS",facet="row")
    if(i>2){
        gr$type3sub= paste0(gr$type3,".",gsub("[0-9]","",seqnames(gr)))
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
# F1 reflects A2vsD5


## another thing i can do is to perform gene-associated genomic annotate for each type of TE using ChIPseeker, which shows the composition of TE annation in regions off gene body, 3', 5', intergenic, etc.
# but it doesn't feel to add more



#############################
## Aggregation per species ##
#############################
source("aggreNvisual.FUN.r")
source("/work/LAS/jfw-lab/hugj2006/cottonLeaf/chipseeker.utilities.R")

## Prepare genomic annotation
ref = data.frame(genome = c("A2", "D5","F1","AD1"), pattern=c("^A","^D","^F","^M"), txDB =  list.files("refGenomes",pattern="txdb",full.names =T), TE =c("gr.A2","gr.D5","gr.F1","gr.AD1"), subgenome =c(1,1,2,2), nuc = list.files("nucleosome",pattern="*coverage.bw", full=TRUE) )

# Plot aggregation
for(i in 1:4)
{
    # genomic features
    print(ref[i,])
    # load gene annot
    txdb <- loadDb(ref$txDB[i])
    gn = genes(txdb);
    # load TE
    #gr = get(ref$TE[i])
    # exclude scaffold
    if(i==4|i==3){
        gn =gn[grep("scaffold",seqnames(gn),invert=TRUE),]
        #   gr=gr[grep("scf",seqnames(gr),invert=TRUE),]  # exclude scaffold
    }else{
        gn =gn[grep("Chr",seqnames(gn)),]
        #gr=gr[grep("Chr",seqnames(gr)),]  # exclude scaffold
    }
    # gene windows and peak list
    tss = getPoint(gn,"start",be=3000, af=1000)
    tts = getPoint(gn,"end",be=1000, af=3000)
    grl <- split(gr, gr$type3)
    
    # list all BWs
    bws = list.files("aggregation",pattern=as.character(ref$pattern[i]),full=T)
    
    pdf(paste0("aggregationPlots/chromationProfile.",ref$genome[i],".pdf"))
    
    # full occupancy
    title= paste0(ref$genome[i],": nucleosome occupancy")
    bw.files = grep("full",bws,value=TRUE)
    print(bw.files)
    tag = gsub("aggregation/|_.*","",bw.files)
    plotTSSnTTS(bw.files,tss,tts,tag,title)
    
    # small size fractioned occupancy
    title= paste0(ref$genome[i],": 0-130bp fragment occupancy")
    bw.files = grep("0-130",bws,value=TRUE)
    print(bw.files)
    tag = gsub("aggregation/|_.*","",bw.files)
    plotTSSnTTS(bw.files,tss,tts,tag,title)
    
    # big size fractioned occupancy
    title= paste0(ref$genome[i],": 130-260bp fragment occupancy")
    bw.files = grep("130-260",bws,value=TRUE)
    print(bw.files)
    tag = gsub("aggregation/|_.*","",bw.files)
    plotTSSnTTS(bw.files,tss,tts,tag,title)
    
    # nucleosome
    title= paste0(ref$genome[i],": 3bp nucleosome centers")
    bw.files = grep("mnase",bws,value=TRUE)
    print(bw.files)
    tag = gsub("aggregation/|_.*","",bw.files)
    plotTSSnTTS(bw.files,tss,tts,tag,title)
    
    # MOA
    title= paste0(ref$genome[i],": MOA frenters")
    bw.files = grep("moa",bws,value=TRUE)
    print(bw.files)
    tag = gsub("aggregation/|_.*","",bw.files)
    plotTSSnTTS(bw.files,tss,tts,tag,title)
    
    dev.off()
}


# subsection DNS scores for different genomic regions

oepness 
