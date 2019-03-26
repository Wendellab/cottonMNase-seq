# ChIPseeker is quite easy to use for peak annotation
###########################################
## Genomic annotation of segments with ChIPseeker
# 1. visualizing segment locations on chromosomes
# 2. plot coverage profile and heatmap over TSS
# 3. annotate genomic location of segments with nearby genes.

# mkdir SegAnnotation
# cd SegAnnotation
# module load r-udunits2
# module load r/3.5.0-py2-ufvuwmm
# R

# Load libraries
library(GenomicFeatures)
library(ChIPseeker) # module load r-udunits2
library(data.table)

## Load data
sampleInfo = data.frame( file= list.files("../isegv1.3.2",pattern="bc6.0.Fus.bed",full.names =T),
genome = c("A2", "D5","F1","AD1"),
txDB =  list.files("../refGenomes",pattern="txdb",full.names =T),
subgenome =c(1,1,2,2),
bigwig= list.files("../iseg", pattern="bg",full.names =T)[1:4])
# print info
sampleInfo

## run peak annotation
for(i in 1:4)
{
    print(sampleInfo[i,])

    filePath = as.character(sampleInfo$file[i])
    pdf(gsub(".*/","",gsub("_.*","_plotAnno.pdf",filePath)))
    
    # read peaks
    gr = readPeakFile(filePath)
    gr$V4 = ifelse(gr$V9=="20,20,255", "MSF","MRF")
    grl = split(gr, gr$V4)
    
    # Plot covplot - sample top 1% for visualization, otherwise completely satuated
    covplot(grl$MSF, weightCol="height.abs", title = "MSFs over chromosomes",lower=quantile(gr$height.abs,0.99))
    covplot(grl$MRF, weightCol="height.abs", title = "MRFs over chromosomes",lower=quantile(gr$height.abs,0.99))
    
    # load txdb
    dbPath= as.character(sampleInfo$txDB[i])
    txdb = loadDb(dbPath)
    # txdb = loadDb("refGenomes/txdb.TM1saski.sqlite")
    txdb
    
    # split At and Dt
    if(sampleInfo$subgenome[i]==2)
    {
        gr$genome = gsub("..$","",seqnames(gr))
        grl = list(gr[gr$V4=="MSF" & gr$genome=="A"], gr[gr$V4=="MSF" & gr$genome=="D"], gr[gr$V4=="MRF" & gr$genome=="A"], gr[gr$V4=="MRF" & gr$genome=="D"] )
        names(grl) <- c("MSF.A","MSF.D","MRF.A","MRF.D")
    }
    
    # plot coverage profile and heatmap over TSS
    # Prepare the promotor regions
    promoter = getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
    # Calculate the tag matrix
    tagMatrixList = lapply(grl, getTagMatrix, windows=promoter)
    ## Profile plots
    plotAvgProf(tagMatrixList, xlim=c(-1000, 1000), conf=0.95,resample=1000, facet="row")
    # Or simply from file
    # plotAvgProf2(segList, TxDb=txdb, upstream=1000, downstream=1000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
    # Plot heatmap
    tagHeatmap(tagMatrixList, xlim=c(-1000, 1000), color=NULL)
    
    # Annotation of peak location with genes
    peakAnnoList <- lapply(grl, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
    peakAnnoList
    # comparative plots
    plotAnnoBar(peakAnnoList)
    plotDistToTSS(peakAnnoList, title="Distribution of MSFs/MRFs relative to TSS")
    # single plot, not very useful though
    lapply(peakAnnoList, upsetplot)
    lapply(peakAnnoList, plotAnnoPie)
    # lapply(peakAnnoList, vennpie)
    dev.off()
    
    
    # save annotation results
    assign( paste0("peakAnnoList.",sampleInfo$genome[i]), peakAnnoList)

}
ls(pattern="peakAnnoList.")
save(list=c("sampleInfo", ls(pattern="peakAnnoList.")), file="genomicAnnotation.rdata")

## annotation list content
names(peakAnnoList.AD1)  # MSF.A, MSF.D, MRF.A, MRF.D
peakAnnoList.AD1$"MSF.D"@anno
peakAnnoList.AD1$"MSF.D"@detailGenomicAnnotation
peakAnnoList.AD1$"MSF.D"@annoStat
peakAnnoList.AD1$"MSF.D"@peakNum

## ref genome info
sl=list()
sl.nonzero = list()
for(i in 1:4)
{
    bg =as.character(sampleInfo$bigwig[i])
    dt=fread(bg,sep="\t")
    chrLen = tapply(dt$V3,dt$V1,max)
    dt.nonzero= dt[dt$V4!=0,]
    chrLen.nonzero = tapply(dt.nonzero$V3-dt.nonzero$V2,dt.nonzero$V1,sum)
    # split At and Dt
    if(sampleInfo$subgenome[i]==2)
    {
        tag =as.character(sampleInfo$genome[i])
        g = factor(gsub("..$","",names(chrLen)))
        sl[[paste0(tag,".At")]] =  sum(chrLen[g=="A"])
        sl[[paste0(tag,".Dt")]] =  sum(chrLen[g=="D"])
        sl.nonzero[[paste0(tag,".At")]] =  sum(chrLen.nonzero[g=="A"])
        sl.nonzero[[paste0(tag,".Dt")]] =  sum(chrLen.nonzero[g=="D"])
    }else{
        tag =as.character(sampleInfo$genome[i])
        sl[[tag]] = sum(chrLen)
        sl.nonzero[[tag]] = sum(chrLen.nonzero)
    }
}
# chr length
sl
# total nonzero region
sl.nonzero


## collect peak annotation results
MSFpeakAnnoList = list( A2 = peakAnnoList.A2$"MSF", D5 = peakAnnoList.D5$"MSF", F1.At = peakAnnoList.F1$"MSF.A", F1.Dt = peakAnnoList.F1$"MSF.D", AD1.At = peakAnnoList.AD1$"MSF.A", AD1.Dt = peakAnnoList.AD1$"MSF.D")
MRFpeakAnnoList = list( A2 = peakAnnoList.A2$"MRF", D5 = peakAnnoList.D5$"MRF", F1.At = peakAnnoList.F1$"MRF.A", F1.Dt = peakAnnoList.F1$"MRF.D", AD1.At = peakAnnoList.AD1$"MRF.A", AD1.Dt = peakAnnoList.AD1$"MRF.D")

save(sampleInfo, MSFpeakAnnoList, MRFpeakAnnoList, sl, sl.nonzero, file="iseg/genomicAnnotation.byType.rdata")

##
l = load("iseg/genomicAnnotation.byType.rdata")

## summarize peak number, region, ref genome,
# MSF
MSFsummary= data.frame( peakNum = unlist(lapply(MSFpeakAnnoList,function(x)x@peakNum)),
peakRegion = unlist(lapply(MSFpeakAnnoList,function(x){sum(width(x@anno))})),
peakRegionProportion = NA,
refGenome = unlist(sl),
refNonzero = unlist(sl.nonzero))
MSFsummary$peakRegionProportion = MSFsummary$peakRegion/ MSFsummary$refGenome
MSFsummary
# MRF
MRFsummary= data.frame( peakNum = unlist(lapply(MRFpeakAnnoList,function(x)x@peakNum)),
peakRegion = unlist(lapply(MRFpeakAnnoList,function(x){sum(width(x@anno))})),
peakRegionProportion = NA,
refGenome = unlist(sl),
refNonzero = unlist(sl.nonzero))
MRFsummary$peakRegionProportion = MRFsummary$peakRegion/ MRFsummary$refGenome
MRFsummary


## plot actual area
# functions
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
getPeakInfo = function(peakAnno){
    anno =peakAnno@anno
    catLevel = c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Exon", "Intron", "Downstream (<=300)","Distal Intergenic" , "5' UTR","3' UTR")
    category = factor(gsub("Downstream.*","Downstream (<=300)",gsub("on [(].*","on",anno$annotation)),levels = catLevel)
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
    print(p)
}
#
df = dflist2df(lapply(MSFpeakAnnoList, getPeakInfo))
df$Feature= factor(df$Feature, rev(c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Exon", "Intron", "Downstream (<=300)","Distal Intergenic" , "5' UTR","3' UTR")))
df$Sample=factor(df$Sample, levels=rev(c("A2","D5","F1.At","F1.Dt","AD1.At","AD1.Dt")))
pdf("iseg/plotAnnotMSF.pdf")
plotPeakInfo(df, y="RegionSize",ylab ="Peak region (bp)", title = "MSF")
plotPeakInfo(df, y="RegionSizePerc",ylab ="Peak region %", title = "MSF")
plotPeakInfo(df, y="Number",ylab ="Peak number", title = "MSF")
plotPeakInfo(df, y="NumberPerc",ylab ="Peak number %", title = "MSF")
plotAnnoBar(MSFpeakAnnoList, title="MSFs Distribution")
dev.off()
# because of the genome size difference, the same 1% of genomic region in D genome is more likely to be near genes than those in A genome. How to correct for the genome size differences, in order to tell whether MSFs are differentially located in A vs D genome?

df = dflist2df(lapply(MRFpeakAnnoList, getPeakInfo))
df$Feature= factor(df$Feature, rev(c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Exon", "Intron", "Downstream (<=300)","Distal Intergenic" , "5' UTR","3' UTR")))
df$Sample=factor(df$Sample, levels=rev(c("A2","D5","F1.At","F1.Dt","AD1.At","AD1.Dt")))
pdf("iseg/plotAnnotMRF.pdf")
plotPeakInfo(df, y="RegionSize",ylab ="Peak region (bp)", title = "MRF")
plotPeakInfo(df, y="RegionSizePerc",ylab ="Peak region %", title = "MRF")
plotPeakInfo(df, y="Number",ylab ="Peak number", title = "MRF")
plotPeakInfo(df, y="NumberPerc",ylab ="Peak number %", title = "MRF")
plotAnnoBar(MRFpeakAnnoList, title="MRFs Distribution")
dev.off()


#########################################################################################
## Inspect ortho-homoeolog quadruplets for the presence of MSFs in 1kb promoter region ##
#########################################################################################

l = load("iseg/genomicAnnotation.byType.rdata")

## collect gene containing MSFs in 1kb promoter
getGenes=function(csAnno, feature="Promoter (<=1kb)")
{
    anno = csAnno@anno
    genes = anno$geneId[anno$annotation==feature]
    return(genes)
}
geneL = lapply(MSFpeakAnnoList, getGenes)
str(geneL)
# List of 6
#  $ A2    : chr [1:12009] "Ga01G0007" "Ga01G0012" "Ga01G0017" "Ga01G0018" ...
#  $ D5    : chr [1:10128] "Gorai.001G001000" "Gorai.001G001700" "Gorai.001G002500" "Gorai.001G002600" ...
#  $ F1.At : chr [1:8452] "Ga01G0003" "Ga01G0005" "Ga01G0005" "Ga01G0005" ...
#  $ F1.Dt : chr [1:10035] "Gorai.001G000500" "Gorai.001G001300" "Gorai.001G002400" "Gorai.001G002900" ...
#  $ AD1.At: chr [1:17351] "Gohir.A01G000100" "Gohir.A01G000300" "Gohir.A01G000300" "Gohir.A01G000500" ...
#  $ AD1.Dt: chr [1:16626] "Gohir.D01G000400" "Gohir.D01G000700" "Gohir.D01G000800" "Gohir.D01G000900" ...

## load 20362 ortholog quadruplets derived from Justin's orthoMCL results version3_081618
ogQ<-read.table("/lss/research/jfw-lab/Projects/MNase-seq/orthohomoeologQuadruplets101218.txt", sep="\t", header=TRUE)
head(ogQ)
# correspondce
og6 = data.frame(groupID = ogQ$groupID, A2=ogQ$A2, D5=ogQ$D5, F1.At=ogQ$A2, F1.Dt=ogQ$D5, AD1.At=ogQ$At, AD1.Dt=ogQ$Dt)
# within in each quadruplets, which contains MSFs and which doesn't
paMSF = data.frame(groupID = ogQ$groupID, A2=ogQ$A2 %in% geneL$A2, D5=ogQ$D5 %in% geneL$D5, F1.At=ogQ$A2 %in% geneL$F1.At, F1.Dt=ogQ$D5 %in% geneL$F1.Dt, AD1.At=ogQ$At %in% geneL$AD1.At, AD1.Dt=ogQ$Dt %in% geneL$AD1.Dt)
# convert to 6 digit string code
paMSFn = apply(apply(paMSF[,-1],2,as.numeric), 1, function(x)paste0(x,collapse=""))
# examine the pattern of p/a
ddply(as.data.frame(paMSF[,-1]),.(A2,D5,F1.At,F1.Dt,AD1.At,AD1.Dt),nrow)
table(paMSFn)
