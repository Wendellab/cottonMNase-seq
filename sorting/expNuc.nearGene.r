library(genomation)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
options(scipen=999)

# FUN
getTSS=function(txdb, include=NULL, exclude.pattern=NULL, be=1000, af=1000){
    tss =  genes(txdb)
    if(!is.null(exclude.pattern)) {tss = tss[grep(exclude.pattern,tss$gene_id,invert=T)]}
    if(!is.null(include)) {tss = tss[tss$gene_id %in% include]}
    # end of the + strand genes must be equalized to start pos
    end(tss[strand(tss)=="+",])  =start(tss[strand(tss)=="+",])
    # start of the - strand genes must be equalized to end pos
    start(tss[strand(tss)=="-",])=end(tss[strand(tss)=="-",])
    # define flanking region
    start(tss)=start(tss)-be
    end(tss)  = end(tss)+ af
    # remove duplicated TSSes ie alternative transcripts
    # this keeps the first instance and removes duplicates
    tss=tss[!duplicated(tss),]
    seqlevels(tss)=unique(as.character(seqnames(tss)))
    return(tss)
}
prepWindows = function(goraiL, genome){
    exp=ogExp[[genome]]
    TSS = getTSS(loadDb(txdbs[which(genomes==genome)]))
    grl = lapply(goraiL, function(goraiID){TSS[TSS$gene_id %in% exp[goraiID,"id"]]})
}
plotQuantile = function(bw, main="Profile of expression quantiles"){
    sml = list(Q0 = ScoreMatrix(target = bw, windows = Q0, type="bigWig",strand.aware=TRUE),
    Q1 = ScoreMatrix(target = bw, windows = Q1, type="bigWig",strand.aware=TRUE),
    Q2 = ScoreMatrix(target = bw, windows = Q2, type="bigWig",strand.aware=TRUE),
    Q3 = ScoreMatrix(target = bw, windows = Q3, type="bigWig",strand.aware=TRUE),
    Q4 = ScoreMatrix(target = bw, windows = Q4, type="bigWig",strand.aware=TRUE) )
    names(sml) =c("Q0","Q1","Q2","Q3","Q4")
    sml= new("ScoreMatrixList",sml)
    plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="",dispersion = "se", xlab="bases around TSS", main=main)
    abline(v=0,lty=2)
}
plotBWsWindows = function(bws, windowList, compare=c("bw","window","paired"), xlab="bases around TSS", main=""){
    nbw=length(bws)
    if(compare=="window"){
        for(b in 1:nbw)
        {
            sml=list()
            for(w in names(windowList)){
                sml[[w]] =  ScoreMatrix(target = bws[b], windows = windowList[[w]], type="bigWig",strand.aware=TRUE)
            }
            sml= new("ScoreMatrixList",sml)
            plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="",dispersion = "se", xlab=xlab, main=paste(main, names(bws)[b]))
            abline(v=0,lty=2)
        }
    }
    if(compare=="bw"){
        for(w in names(windowList)){
            sml=ScoreMatrixList(target = bws, windowList[[w]], type="bigWig",strand.aware=TRUE)
            names(sml) = names(bws)
            plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="",dispersion = "se", xlab=xlab, main=w)
            abline(v=0,lty=2)
        }
    }
    if(compare=="paired"){
        if(nbw!=length(windowList)) {message("Error: input bws and windows are not paired!!")}else
        {
            sml=list()
            for(b in 1:nbw)
            {
                sml[[b]] =  ScoreMatrix(target = bws[b], windows = windowList[[b]], type="bigWig",strand.aware=TRUE)
            }
            sml= new("ScoreMatrixList",sml)
            plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(bws), ylab="",dispersion = "se", xlab=xlab, main=main)
            abline(v=0,lty=2)
        }
    }
}



grl.A2 = preWindows(goraiL=str(split(rownames(res), res$category)),"A2")
grl.D5 = preWindows(goraiL=str(split(rownames(res), res$category)),"D5")


# load gene expression results
load("../RNAseq/Ranalysis/expression.Indiv.rdata") # "A2"  "D5"  "F1"  "AD1"
load("../RNAseq/Ranalysis/expression.Dref.rdata") # "all"  "info" "A2"   "D5"   "F1"   "AD1"

# load ortholog info
ogQ<-read.table("/lss/research/jfw-lab/Projects/MNase-seq/orthohomoeologQuadruplets101218.txt", sep="\t", header=TRUE)
head(ogQ)
ogL =list()
ogL$AD1 <- data.frame(id = c(as.character(ogQ$At),as.character(ogQ$Dt)), index=c(paste0(ogQ$D5,".A"),paste0(ogQ$D5,".D")))
ogL$F1 <- data.frame(id = c(as.character(ogQ$A2),as.character(ogQ$D5)), index=c(paste0(ogQ$D5,".A"),paste0(ogQ$D5,".D")))
ogL$A2 <- data.frame(id = as.character(ogQ$A2), index=as.character(ogQ$D5))
ogL$D5 <- data.frame(id = as.character(ogQ$D5), index=as.character(ogQ$D5))

# prepare loop
genomes =c("AD1","F1","A2","D5")
tags=c("M","F","A","D")
dirs = c("../mappingMnew","../mappingF","../mappingA_new","../mappingD")
txdbs=c("../refGenomes/txdb.TM1saski.sqlite","../refGenomes/txdb.F1.sqlite","../refGenomes/txdb.A2du.sqlite","../refGenomes/txdb.D5.sqlite")
TSSs = sapply(txdbs, getTSS)
ogExp = list()

##############################################################################
## plot D/H/L expression quantiles for each species, save up 500 bp MS area ##
##############################################################################
for(i in 1:4)
{
    genome = genomes[i]
    dir=dirs[i]
    txdb=loadDb(txdbs[i])
    TSS = getTSS(txdb)
    
    # bws
    bws = list.files(dir,pattern="^(.c|...n).*genic.depth_qnorm.bw",full.names=T)
    bwsID= gsub("_.*|.*[/]","",bws)

    # rpkm
    rpkm = get(genome)[["rpkm"]]
    
    # ortholog expression
    exp=rpkm[rownames(rpkm) %in% ogL[[genome]]$index,]
    exp$id = ogL[[genome]]$id[match(rownames(exp), ogL[[genome]]$index)]
    
    # quantile
    log2mean = apply(exp[,1:ncol(rpkm)],1,function(x)log2(mean(x)+1))
    breaks=quantile(log2mean[log2mean>0], probs=0:4/4)
    exp$quantile=cut(log2mean, breaks, include.lowest=TRUE, labels=FALSE)
    exp$quantile[is.na(exp$quantile)]=0
    print(table(exp$quantile))
    exp$log2mean=log2mean
    
    # upstream 500bp area of D>0
    gr=import(bws[grep("..D",bwsID)])
    up=getTSS(txdb, be=500, af=0, include=exp$id)
    ov = findOverlaps(up, gr)
    area = aggregate(subjectHits(ov),by=list(queryHits(ov)), function(x){scores=gr[x]$score; sum(scores[scores>0])})
    exp$MSarea = area$x[match(exp$id,up$gene_id)]
    
    # save
    ogExp[[genome]]=exp
    
    # prepare windows
    Q0 = TSS[TSS$gene_id %in% exp$id[exp$quantile==0],]
    Q1 = TSS[TSS$gene_id %in% exp$id[exp$quantile==1],]
    Q2 = TSS[TSS$gene_id %in% exp$id[exp$quantile==2],]
    Q3 = TSS[TSS$gene_id %in% exp$id[exp$quantile==3],]
    Q4 = TSS[TSS$gene_id %in% exp$id[exp$quantile==4],]
    
    # plot
    pdf(paste0("plotQuantile.",genome,".pdf"))
    plotQuantile(bws[1],main=paste0("Profile of expression quantiles: ",bwsID[1]))
    plotQuantile(bws[2],main=paste0("Profile of expression quantiles: ",bwsID[2]))
    plotQuantile(bws[3],main=paste0("Profile of expression quantiles: ",bwsID[3]))
    dev.off()
}
str(ogExp)

load("../RNAseq/Ranalysis/diffExprPattern.rdata")
# "A"        "B"        "Bp"       "W"        "total"    "AD1vsA2"
# "AD1vsD5"  "AD1vsMid" "F1vsA2"   "F1vsD5"   "F1vsMid"  "MidvsA2"  "MidvsD5"
load("../RNAseq/Ranalysis/regPattern.rdata") # "res"           "dominance.AD1" "dominance.F1"
# restrict
res=res[as.character(ogQ$D5),]


####################################
## differential expression A vs D ##
####################################

# AD1
AoD = rownames(res)[which(res$Bp>0 & res$Bp.padj<0.05 & !is.na(res$Bp.padj))]
DoA = rownames(res)[which(res$Bp<0 & res$Bp.padj<0.05 & !is.na(res$Bp.padj))]
eqAD = setdiff(rownames(res),c(AoD,DoA))
grl.A=prepWindows(list(AoD=paste0(AoD,".A"), DoA=paste0(DoA,".A"), eqAD=paste0(eqAD,".A")),"AD1")
grl.D=prepWindows(list(AoD=paste0(AoD,".D"), DoA=paste0(DoA,".D"), eqAD=paste0(eqAD,".D")),"AD1")
bws=c("../mappingMnew/McD_q20.genic_chr.size_w20.genic.depth_qnorm.bw")
names(bws)="McD"
pdf("plotHomoeDE.AD1.pdf")
plotBWsWindows(bws,list(At=grl.A$AoD,Dt=grl.D$AoD), compare="window", paste0("At>Dt: ",length(AoD)))
plotBWsWindows(bws,list(At=grl.A$DoA,Dt=grl.D$DoA), compare="window",paste0("At>Dt: ",length(DoA)))
plotBWsWindows(bws,list(At=grl.A$eqAD,Dt=grl.D$eqAD), compare="window",paste0("At=Dt: ",length(eqAD)))
dev.off()

# F1
AoD = rownames(res)[which(res$B>0 & res$B.padj<0.05 & !is.na(res$B.padj))]
DoA = rownames(res)[which(res$B<0 & res$B.padj<0.05 & !is.na(res$B.padj))]
eqAD = setdiff(rownames(res),c(AoD,DoA))
grl.A=prepWindows(list(AoD=paste0(AoD,".A"), DoA=paste0(DoA,".A"), eqAD=paste0(eqAD,".A")),"F1")
grl.D=prepWindows(list(AoD=paste0(AoD,".D"), DoA=paste0(DoA,".D"), eqAD=paste0(eqAD,".D")),"F1")
bws=c("../mappingF/FcD_q20.genic_chr.size_w20.genic.depth_qnorm.bw")
names(bws)="FcD"
pdf("plotHomoeDE.F1.pdf")
plotBWsWindows(bws,list(At=grl.A$AoD,Dt=grl.D$AoD), compare="window", paste0("At>Dt: ",length(AoD)))
plotBWsWindows(bws,list(At=grl.A$DoA,Dt=grl.D$DoA), compare="window",paste0("At>Dt: ",length(DoA)))
plotBWsWindows(bws,list(At=grl.A$eqAD,Dt=grl.D$eqAD), compare="window",paste0("At=Dt: ",length(eqAD)))
dev.off()

# A2 vs D5
AoD = rownames(res)[which(res$A>0 & res$A.padj<0.05 & !is.na(res$A.padj))]
DoA = rownames(res)[which(res$A<0 & res$A.padj<0.05 & !is.na(res$A.padj))]
eqAD = setdiff(rownames(res),c(AoD,DoA))
bws=c("../mappingA_new/A6Dn_q20.genic_chr.size_w20.genic.depth_qnorm.bw","../mappingD/DcD_q20.genic_chr.size_w20.genic.depth_qnorm.bw" )
names(bws)=c("A2","D5")
grl.A=prepWindows(list(AoD=AoD, DoA=DoA, eqAD=eqAD),"A2")
grl.D=prepWindows(list(AoD=AoD, DoA=DoA, eqAD=eqAD),"D5")
pdf("plotHomoeDE.A2D5.pdf")
plotBWsWindows(bws,list(At=grl.A$AoD,Dt=grl.D$AoD), compare="paired", paste0("A2>D5: ",length(AoD)))
plotBWsWindows(bws,list(At=grl.A$DoA,Dt=grl.D$DoA), compare="paired",paste0("A2>D5: ",length(DoA)))
plotBWsWindows(bws,list(At=grl.A$eqAD,Dt=grl.D$eqAD), compare="paired",paste0("A2=D5: ",length(eqAD)))
dev.off()


## try again above using qnormed bws of A2, D5, F1.At, F1.Dt
library(travis)
options(threads=8)
options(verbose=T)
options(threads=8)
options(verbose=T)

setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf")
system("grep '^A' mappingF/FcD_q20.genic_chr.size_w20.genic.depth_qnorm.bg| sed 's/^A/Chr/g' >nearGenic/FcDa_q20.genic_chr.size_w20.genic.depth_qnorm.bg")
system("grep '^D' mappingF/FcD_q20.genic_chr.size_w20.genic.depth_qnorm.bg| sed 's/^D/Chr/g' >nearGenic/FcDd_q20.genic_chr.size_w20.genic.depth_qnorm.bg")
bgA = "mappingA_new/A6Dn_q20.genic_chr.size_w20.genic.depth_qnorm.bg"
bgD = "mappingD/DcD_q20.genic_chr.size_w20.genic.depth_qnorm.bg"
bgAt = "nearGenic/FcDa_q20.genic_chr.size_w20.genic.depth_qnorm.bg"
bgDt = "nearGenic/FcDd_q20.genic_chr.size_w20.genic.depth_qnorm.bg"
bgM = "mappingMnew/McD_q20.genic_chr.size_w20.genic.depth_qnorm.bg"
qbgs=bgQuantileNorm(c(bgD,bgA,bgAt,bgDt,bgM))
names(qbgs) = c("bgD","bgA","bgAt","bgDt","bgM")
# convert bedGraph files to bigWig
bwA=bedGraphToBigWig(qbgs['bgA'],"mappingA_new/chr.size.txt")
bwAt=bedGraphToBigWig(qbgs['bgAt'],"mappingA_new/chr.size.txt")
bwD=bedGraphToBigWig(qbgs['bgD'],"mappingD/chr.size.txt")
bwDt=bedGraphToBigWig(qbgs['bgDt'],"mappingD/chr.size.txt")
bwM=bedGraphToBigWig(qbgs['bgM'],"mappingMnew/chr.size.txt")
sapply(c(qbgs, bwA, bwAt, bwD, bwDt, bwM), function(x)system(paste0("mv ",x," nearGenic/")))
setwd("nearGenic/")

# F1
AoD = rownames(res)[which(res$B>0 & res$B.padj<0.05 & !is.na(res$B.padj))]
DoA = rownames(res)[which(res$B<0 & res$B.padj<0.05 & !is.na(res$B.padj))]
eqAD = setdiff(rownames(res),c(AoD,DoA))
bws=c(bwAt,bwDt)
names(bws)=c("At","Dt")
grl.A=prepWindows(list(AoD=AoD, DoA=DoA, eqAD=eqAD),"A2")
grl.D=prepWindows(list(AoD=AoD, DoA=DoA, eqAD=eqAD),"D5")
pdf("plotHomoeDE.F1.qnorm.pdf")
plotBWsWindows(bws,list(At=grl.A$AoD,Dt=grl.D$AoD), compare="paired", paste0("At>Dt: ",length(AoD)))
plotBWsWindows(bws,list(At=grl.A$DoA,Dt=grl.D$DoA), compare="paired",paste0("At>Dt: ",length(DoA)))
plotBWsWindows(bws,list(At=grl.A$eqAD,Dt=grl.D$eqAD), compare="paired",paste0("At=Dt: ",length(eqAD)))
dev.off()

# A2 vs D5
AoD = rownames(res)[which(res$A>0 & res$A.padj<0.05 & !is.na(res$A.padj))]
DoA = rownames(res)[which(res$A<0 & res$A.padj<0.05 & !is.na(res$A.padj))]
eqAD = setdiff(rownames(res),c(AoD,DoA))
bws=c(bwA,bwD)
names(bws)=c("A2","D5")
grl.A=prepWindows(list(AoD=AoD, DoA=DoA, eqAD=eqAD),"A2")
grl.D=prepWindows(list(AoD=AoD, DoA=DoA, eqAD=eqAD),"D5")
pdf("plotHomoeDE.A2D5.qnorm.pdf")
plotBWsWindows(bws,list(At=grl.A$AoD,Dt=grl.D$AoD), compare="paired", paste0("A2>D5: ",length(AoD)))
plotBWsWindows(bws,list(At=grl.A$DoA,Dt=grl.D$DoA), compare="paired",paste0("A2>D5: ",length(DoA)))
plotBWsWindows(bws,list(At=grl.A$eqAD,Dt=grl.D$eqAD), compare="paired",paste0("A2=D5: ",length(eqAD)))
dev.off()

###################
## Hybridization ##
###################

# A2 vs At
bws=c(bwA,bwAt)
names(bws)=c("A2","At")
pdf("plotHybridization.A2At.qnorm.pdf")
grl = prepWindows(list(all=rownames(res)),"A2")

plotBWsWindows(bws,grl, compare="bw")
plotBWsWindows(bws,grl, compare="bw")
plotBWsWindows(bws,grl, compare="bw")
dev.off()

# D5 vs Dt
Hr > 0
Hr < 0


# F1
AoD = rownames(res)[which(res$B>0 & res$B.padj<0.05 & !is.na(res$B.padj))]
DoA = rownames(res)[which(res$B<0 & res$B.padj<0.05 & !is.na(res$B.padj))]
eqAD = setdiff(rownames(res),c(AoD,DoA))
titles = c(paste0("At>Dt: ",length(AoD)), paste0("At<Dt: ",length(DoA)), paste0("At=Dt: ",length(eqAD)))
grl.F1A =prepWindows(list(AoD=AoD, DoA=DoA, eqAD=eqAD),"A2")
grl.F1D =prepWindows(list(AoD=AoD, DoA=DoA, eqAD=eqAD),"D5")
pdf("plotHomoeDE.F1.qnorm.pdf")
mapply(plotA2D5, bwAt, bwDt, grl.F1A, grl.F1D, titles)
dev.off()
# A2D5
AoD = rownames(res)[which(res$A>0 & res$A.padj<0.05 & !is.na(res$A.padj))]
DoA = rownames(res)[which(res$A<0 & res$A.padj<0.05 & !is.na(res$A.padj))]
eqAD = setdiff(rownames(res),c(AoD,DoA))
titles = c(paste0("A2>D5: ",length(AoD)), paste0("A2<D5: ",length(DoA)), paste0("A2=D5: ",length(eqAD)))
grl.F1A =prepWindows(list(AoD=AoD, DoA=DoA, eqAD=eqAD),"A2")
grl.F1D =prepWindows(list(AoD=AoD, DoA=DoA, eqAD=eqAD),"D5")
pdf("plotHomoeDE.A2D5.qnorm.pdf")
mapply(plotA2D5, bwA, bwD, grl.F1A, grl.F1D, titles)
dev.off()

## generate log2ratio At/A2, Dt/D5
cd ../nearGenic/hybridization/
bigwigCompare -b1 FDcD_q20.genic_chr.size_w20.genic.depth_qnorm.bw -b2 DcD_q20.genic_chr.size_w20.genic.depth_qnorm.bw --binSize 20 -o FDvsD5.log2ratio.bw
bigwigCompare -b1 FAcD_q20.genic_chr.size_w20.genic.depth_qnorm.bw -b2 A6Dn_q20.genic_chr.size_w20.genic.depth_qnorm.bw --binSize 20 -o FAvsA2.log2ratio.bw
#R script
bwDr = "FDvsD5.log2ratio.bw"
bwAr = "FAvsA2.log2ratio.bw"
plotA2D5(bwAr, bwDr, grl.A, grl.D)



