library(genomation)
library(GenomicRanges)
library(data.table)
library(ggbio)
options(scipen=10)

bgs = list.files("iseg/", pattern="bg")
sl=list()
nonzero = list()
for(bg in bgs)
{
    flag = gsub("_.*","",bg)
    filePath = paste0("iseg/",bg)
    dt=fread(filePath,sep="\t")
    sl[[flag]] = tapply(dt$V3,dt$V1,max)
    dt$width = dt$V3-dt$V2
    select= (dt$V4!=0 )
    x = aggregate(dt$width[select],by=list(dt$V1[select]),sum)
    nonzero[[flag]] =  x$x
    names(nonzero[[flag]]) =x$Group.1
}
# chr length
sl
# nonzero region length
nonzero

# after normalization, check "manifest.txt"
# MAD= 0.3038
# SD = 1.4826 * MAD
SD=0.045071

## plot function
log10hist <- function(x, main = filePath, xlab="",breaks = "Sturges")
{
    h<-hist(x,plot=F,breaks = breaks);
    h$counts[h$counts>0] <- log10(h$counts[h$counts>0])
    plot(h, ylab='log10(Freq)', xlab =xlab,main=main)
}

## prepare working list
# previous version
# fileL = grep("Fus.txt",list.files("iseg/"),value=TRUE)
# v1.3.2
filedir = "isegv1.3.2/"
fileL = grep("0.txt",list.files(filedir),value=TRUE)

## loop import and report stats for each flag
flags = unique(gsub("_.*","",fileL))
for(flag in flags){
    # collect files under same flag
    flagFiles = fileL[grep(flag,fileL)]
    # prep result table
    res = data.frame(ChrName=names(sl[[flag]]), Length = sl[[flag]], Chromosome = 1:length(sl[[flag]]))
    res.nz = data.frame(ChrName=names(nonzero[[flag]]), Length = nonzero[[flag]], Chromosome = 1:length(nonzero[[flag]]))
    res.peak = data.frame(ChrName=names(sl[[flag]]), Length = sl[[flag]], Chromosome = 1:length(sl[[flag]]))
    # loop each file
    for(file in flagFiles){
        bc = gsub(".*bc|[.].*","",file)
        # read into grange
        filePath = paste0(filedir,file)
        gr = readGeneric(file=filePath, chr = 1, start = 2, end = 3, strand = NULL, meta.col=list(height=4,tstat=5,pvalue=6), sep=" ")
        gr$pvalue=fread(filePath,sep=" ",select=6)$V6
        gr$dns = ifelse(gr$height>0,"MSF","MRF")
        ## Region bps called as peaks
        sRegion = data.frame(tapply(width(gr), list(seqnames(gr), gr$dns),sum))
        names(sRegion) =paste0(names(sRegion),".bc",bc)
        ## Proportion of regions called
        pRegion = sweep(sRegion,1, nonzero[[flag]],"/")
        ## reduce range to count peaks
        s =reduce(gr[gr$dns=="MSF"]); s$dns = "MSF"
        r =reduce(gr[gr$dns=="MRF"]); r$dns = "MRF"
        rgr = c(s,r)
        ## number of peaks
        nPeaks = data.frame(tapply(as.character(seqnames(rgr)), list(seqnames(rgr), rgr$dns),length))
        names(nPeaks) =paste0(names(nPeaks),".bc",bc)
        ## make report table
        res = cbind(res,sRegion)
        res.nz = cbind(res.nz,pRegion)
        res.peak = cbind(res.peak,nPeaks)
        
        ## plot segment width, using reduced
        pdf(paste0(filedir,"/plotDis.",gsub(".txt","",file),".pdf"))
        # hist of segment width
        log10hist(width(rgr), main = "Segment Width", xlab="bp")
        par(mfrow=c(4,4 ))
        rgrl <- split(rgr, seqnames(rgr))
        for(i in names(rgrl)) log10hist(width(rgrl[[i]]),main=i,xlab="bp")
        ## hist of height
        par(mfrow=c(1,1))
        log10hist(gr$height/SD, main = "Segment Biological Cutoff", xlab="BC", breaks=100)
        par(mfrow=c(4,4))
        for(i in names(grl)) log10hist(grl[[i]]$height/SD,main=i,xlab="BC", breaks=100)
        # Making karyograms and circos plots
        # par(mfrow=c(1,1))
        # autoplot(seqinfo(gr)) + layout_karyogram(gr, geom="point",aes(x=start,y=height,color=dns))
        dev.off()
    }
    
    cname = c(names(res)[1:3],sort(names(res[-(1:3)])))
    out=paste0(filedir,"sumStat.",flag,".txt")
    ## write results
    write.table("# segment bp",file=out,row.names=FALSE, sep="\t",quote=FALSE)
    write.table(res[,cname],file=out, row.names=FALSE, sep="\t",quote=FALSE, append=TRUE)
    write.table("# segment bp %",file=out, row.names=FALSE, sep="\t",quote=FALSE, append=TRUE)
    write.table(res.nz[,cname], file=out, row.names=FALSE, sep="\t",quote=FALSE, append=TRUE)
    write.table("# segment number",file=out, row.names=FALSE, sep="\t",quote=FALSE, append=TRUE)
    write.table(res.peak[,cname],file=out, row.names=FALSE, sep="\t",quote=FALSE, append=TRUE)
}




##################################### delete below
    --
    ## read chr info
    # genome(gr) <- "D5"
    if( unique(seqlevels(gr) == as.character(sl$V1))){
        seqlengths(gr)<-sl$V2}else{
            stop("Chromosome info doesn't match iSeg results.")
        }
    
    ## split by chr and by MSF/MRF
    # seqnames(gr)
    # gr[seqnames(gr)=="Chr01"]
    grl <- split(gr, seqnames(gr))
    grl_pos <- split(gr[gr$height>0], seqnames(gr[gr$height>0]))
    grl_neg <- split(gr[gr$height<0], seqnames(gr[gr$height<0]))


## pair replicates and check overlap
samples = data.frame(matrix(unlist(strsplit(fileL,"_")),byrow=TRUE,ncol=3)[,1:2])
rownames(samples)=fileL
names(samples)=c("sample","BC")
samples$raw = paste0("iseg/",rownames(samples))
samples$bed = paste0("iseg/",gsub("txt","bed",rownames(samples)))
samples

getOverlaps<-function(x){
    gr1<- readGeneric(file=x[1], chr = 1, start = 2, end = 3, strand = NULL, meta.col=list(height=4,tstat=5,pvalue=6),sep=" ")
    gr1_pos<-gr1[gr1$height>0]
    gr1_neg<-gr1[gr1$height<0]
    gr2<- readGeneric(file=x[2], chr = 1, start = 2, end = 3, strand = NULL, meta.col=list(height=4,tstat=5,pvalue=6),sep=" ")
    gr2_pos<-gr2[gr2$height>0]
    gr2_neg<-gr2[gr2$height<0]
    res = data.frame(raw=x,
    segN = c(length(gr1),length(gr2)),
    segN_pos = c(length(gr1_pos),length(gr2_pos)),
    segN_pos_overlap = c( table(countOverlaps(gr1_pos, gr2_pos)>0)["TRUE"],table(countOverlaps(gr2_pos, gr1_pos)>0)["TRUE"] ),
    segN_neg = c(length(gr1_neg),length(gr2_neg)),
    segN_neg_overlap = c( table(countOverlaps(gr1_neg, gr2_neg)>0)["TRUE"],table(countOverlaps(gr2_neg, gr1_neg)>0)["TRUE"] ) )
    return(res)
}
sumRep<-tapply(samples$raw,samples$BC,getOverlaps)
print(sumRep)
sumRep.df<-sumRep[[1]]
for(i in 2:length(sumRep))
{sumRep.df<-rbind(sumRep.df,sumRep[[i]])}


write.table(sumTbl,"iseg/sumTbl.txt",sep="\t",quote=FALSE)
write.table(sumRep.df,"iseg/sumReplicate.txt",sep="\t",quote=FALSE)
