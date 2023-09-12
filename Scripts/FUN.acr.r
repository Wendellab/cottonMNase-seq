# categorized ACRs on the basis of their distance to the nearest gene
# genic gACRs; overlapping a gene
# proximal pACRs; within 2 kb of a gene
# distal dACRs; >2 kb from a gene.

# Load libraries
library(rtracklayer)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(gridExtra)
library(gplots)
library(dplyr)
library(GenomicFeatures)
library(viridis)
library(UpSetR)
library(RColorBrewer)

# function to categorizate peaks and get distrance to nearest genes
annotateACRs=function(peaks,features, distance=2000){
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

# function to annotate peaks located within TEs (completely encompassed by TE intervals)
annotateACRs_TE=function(peaks, tes,upsetplot=TRUE){
    require(reshape2)
    require(UpSetR)
    # add column of TE family annotation, whether each peak overlap with TE
    hit=findOverlaps(peaks,tes)
    d=data.frame(hit)
    d$class =tes$type[subjectHits(hit)]
    type =  aggregate(d$class,list(d$queryHits),function(x)paste0(unique(sort(x)),collapse="; "))
    peaks$TE =NA
    peaks$TE[type$Group.1] =type$x
    # annotation by family, how each TE annotated peak correspond to multiple TE families
    number = dcast(d,queryHits~class,length)
    rownames(number) =number$queryHits
    number=number[,-1]
    number[number>1]=1
    # ACR proportion co in TE families, whether certain TE families enriched by peaks
    en = aggregate(width(tes),list(tes$type),sum)
    ACRinTEsize  = apply(number,2,function(x)sum(x*width(peaks)[as.numeric(rownames(number))]))
    names(en) =c("TE_family","Bp")
    en$TE_family=as.character(en$TE_family)
    en$ACR_in_TE.Bp = 0
    en$ACR_in_TE.Bp = ACRinTEsize[en$TE_family]
    en$ACR_in_TE.Bp[is.na(en$ACR_in_TE.Bp)]=0
    nn=nrow(en)+1
    en[nn,1] ="All"
    en[nn,2:3] = c(sum(width(reduce(tes))),sum(width(reduce(peaks[!is.na(peaks$TE)]))))
    en$ACR_perc = en$ACR_in_TE.Bp/en$Bp
    row.names(en)=en$TE_family
    #print(en)
    if(upsetplot){
        p=upset(number,nset=nrow(number), order.by="freq",number.angles = 30, point.size = 3.5, line.size = 2, mainbar.y.label = "TE Family Intersections", sets.x.label = "Per TE Family", sets.bar.color = "#56B4E9", text.scale = c(1.3, 1.3, 1, 1, 1.5, 0.65))
        return(list(gr=peaks,by.family=number, enrich = en, plot=p))
    }else{return(list(gr=peaks,by.family=number, enrich = en))}
}
getPeakSum=function(acr_df,chr_size){
    x=data.frame(table(acr_df$type),size=tapply(acr_df$width,list(acr_df$type),sum))
    names(x)=c("type","number","size")
    x$numberP = x$number/sum(x$number)
    x$sizeP = x$size/sum(x$size)
    x$genomeP = x$size/sum(chr_size$V2)
    return(x)
}
getPeakSumTE=function(acr_df,chr_size){
    acr_df$TE=gsub("/.*","",acr_df$TE)
    acr_df$TE[is.na(acr_df$TE)]="nonTE"
    x=data.frame(table(acr_df$TE),size=tapply(acr_df$width,list(acr_df$TE),sum))
    names(x)=c("type","number","size")
    x$numberP = x$number/sum(x$number)
    x$sizeP = x$size/sum(x$size)
    return(x)
}

