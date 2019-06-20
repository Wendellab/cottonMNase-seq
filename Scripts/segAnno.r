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
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/isegRes/SegAnnotation")
# Load libraries
library(GenomicFeatures)
library(ChIPseeker) # module load r-udunits2
library(data.table)
library(ggplot2)

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
peakAnnotation = function(input=sampleInfo, outDir=""){
    for(i in 1:nrow(sampleInfo))
    {
        print(sampleInfo[i,])
        
        filePath = as.character(sampleInfo$file[i])
        pdf(paste0(outDir,gsub("_.*","_plotAnno.pdf",gsub(".*/","",filePath))))
        
        # read peaks
        gr = readPeakFile(filePath)
        gr$V4 = ifelse(gr$V9=="20,20,255", "Pos","Neg")
        grl = split(gr, gr$V4)
        
        # Plot covplot - sample top 1% for visualization, otherwise completely satuated
        ## need height from iseg.py generated fus.txt
        #covplot(grl$MSF, weightCol="height.abs", title = "MSFs over chromosomes",lower=quantile(gr$height.abs,0.99))
        #covplot(grl$MRF, weightCol="height.abs", title = "MRFs over chromosomes",lower=quantile(gr$height.abs,0.99))
        
        # load txdb
        dbPath= as.character(sampleInfo$txDB[i])
        txdb = loadDb(dbPath)
        # txdb = loadDb("refGenomes/txdb.TM1saski.sqlite")
        txdb
        
        # split At and Dt
        if(sampleInfo$subgenome[i]==2)
        {
            gr$genome = gsub("..$","",seqnames(gr))
            gr$type = paste0(gr$V4,".",gr$genome)
            grl = split(gr, gr$type)
        }
        
        # plot coverage profile and heatmap over TSS
        # Prepare the promotor regions
        promoter = getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
        # Calculate the tag matrix
        tagMatrixList = lapply(grl, getTagMatrix, windows=promoter)
        ## Profile plots
        p1=plotAvgProf(tagMatrixList, xlim=c(-1000, 1000), conf=0.95,resample=1000, facet="row")
        # Or simply from file
        # plotAvgProf2(segList, TxDb=txdb, upstream=1000, downstream=1000, xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
        # Plot heatmap
        p2=tagHeatmap(tagMatrixList, xlim=c(-1000, 1000), color=NULL)
        print(p1)
        print(p2)
        
        # Annotation of peak location with genes
        peakAnnoList <- lapply(grl, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
        peakAnnoList
        # comparative plots
        p=plotAnnoBar(peakAnnoList);print(p)
        p=plotDistToTSS(peakAnnoList, title="Distribution of Pos/Neg footprints relative to TSS");print(p)
        # single plot, not very useful though
        lapply(peakAnnoList, upsetplot)
        lapply(peakAnnoList, plotAnnoPie)
        # lapply(peakAnnoList, vennpie)
        dev.off()
        
        # save annotation results
        assign( paste0("peakAnnoList.",sampleInfo$genome[i]), peakAnnoList)
        
    }
    ls(pattern="peakAnnoList.")
    save(list=c("sampleInfo", ls(pattern="peakAnnoList.")), file=paste0(outDir,"genomicAnnotation.rdata"))
    
    ## ref genome info
    sl=list() # chr length
    sl.nonzero = list()  # total nonzero region
    for(i in 1:nrow(sampleInfo))
    {
        bg =as.character(sampleInfo$inputBG[i])
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
    
    ## collect peak annotation results
    # t is either Pos or Neg
    Sums=list()
    for(t in names(peakAnnoList.A2)){
        resL=list( A2 = peakAnnoList.A2[[t]], D5 = peakAnnoList.D5[[t]], F1.At = peakAnnoList.F1[[paste0(t,".A")]], F1.Dt = peakAnnoList.F1[[paste0(t,".D")]], AD1.At = peakAnnoList.AD1[[paste0(t,".A")]], AD1.Dt = peakAnnoList.AD1[[paste0(t,".D")]])
        # save 
        assign(paste0(t, "_peakAnnoList"), resL)
       
        # generate summary
        sumRes = data.frame( peakNum = unlist(lapply(resL,function(x)x@peakNum)), peakRegion = unlist(lapply(resL,function(x){sum(width(x@anno))})), peakRegionProportion = NA, refGenome = unlist(sl), refNonzero = unlist(sl.nonzero))
        sumRes$peakRegionProportion = sumRes$peakRegion/ sumRes$refGenome
        Sums[[t]]=sumRes
       
        # plot summary 
        df = dflist2df(lapply(resL, getPeakInfo))
        df$Feature= factor(df$Feature, rev(c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Exon", "Intron", "Downstream (<=300)","Distal Intergenic" , "5' UTR","3' UTR")))
        df$Sample=factor(df$Sample, levels=rev(c("A2","D5","F1.At","F1.Dt","AD1.At","AD1.Dt")))
        pdf(paste0(outDir,"plotAnnot",t,".pdf"))
        plotPeakInfo(df, y="RegionSize",ylab ="Peak region (bp)", title = t)
        plotPeakInfo(df, y="RegionSizePerc",ylab ="Peak region %", title = t)
        plotPeakInfo(df, y="Number",ylab ="Peak number", title = t)
        plotPeakInfo(df, y="NumberPerc",ylab ="Peak number %", title = t)
        print(plotAnnoBar(resL, title="Footprints Distribution"))
        dev.off()
        # because of the genome size difference, the same 1% of genomic region in D genome is more likely to be near genes than those in A genome. How to correct for the genome size differences, in order to tell whether MSFs are differentially located in A vs D genome?
    }
    print(Sums)
    save(list=c("sampleInfo", "sl", "sl.nonzero","Sums", grep("_peakAnnoList",ls(),value=TRUE)), file=paste0(outDir,"genomicAnnotation.byType.rdata"))

} 



## run DNS segment annotation
system("mkdir dns")
segDir="../isegv1.3.2_032519"
## Load data
sampleInfo = data.frame( file= list.files(segDir,pattern="bc6.0.Fus.bed",full.names =T), genome = c("A2", "D5","F1","AD1"), txDB =  list.files("../../refGenomes",pattern="txdb",full.names =T), subgenome =c(1,1,2,2), inputBG= list.files("..", pattern="bg",full.names =T) )
# print info
sampleInfo
peakAnnotation(input=sampleInfo, outDir="dns/")

## run Small fragment annotation
system("mkdir small")
segDir="../MOA/iseg_v1.3.2_041219"
## Load data
sampleInfo = data.frame( file= paste0(segDir,"/",c("A6Ln_bc7.0.Fus.bed","DcL_bc4.0.Fus.bed","FcL_bc4.0.Fus.bed","McL_bc6.0.Fus.bed")), genome = c("A2", "D5","F1","AD1"), txDB =  list.files("../../refGenomes",pattern="txdb",full.names =T), subgenome =c(1,1,2,2), inputBG= list.files("../MOA", pattern="w20.bg",full.names =T)  )
# print info
sampleInfo
peakAnnotation(input=sampleInfo, outDir="smallL/")
