library(genomation)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
options(scipen=999)
library(travis)
options(threads=8)
options(verbose=T)
library(data.table)
# FUN
getEffectiveRefSize.bg=function(bedGraphFile){
    rpm <- fread(bedGraphFile, select = 4)$V4;
    start <- fread(bedGraphFile, select = 2)$V2;
    end <- fread(bedGraphFile, select = 3)$V3;
    ts <- end-start
    cs <- sum(ts[rpm!=0])
    return(cs)
}
getSD <- function(bg)
{
    val<-fread(bg,select=4)$V4
    return(sd(val))
}
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
            plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="",dispersion = "se", xlab=xlab, main=paste(main,w))
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
dirs = c("/work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingMnew","/work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingF","/work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingA_new","/work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingD")
txdbs=c("/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/txdb.TM1saski.sqlite","/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/txdb.F1.sqlite","/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/txdb.A2du.sqlite","/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/txdb.D5.sqlite")
ogExp = list()

##############################################################################
## plot D/H/L expression quantiles for each species, save up 500 bp MS area ##
##############################################################################
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/nearGenic")
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

####################################
## differential expression A vs D ##
####################################
load("../RNAseq/Ranalysis/diffExprPattern.rdata")
# "A"        "B"        "Bp"       "W"        "total"    "AD1vsA2"
# "AD1vsD5"  "AD1vsMid" "F1vsA2"   "F1vsD5"   "F1vsMid"  "MidvsA2"  "MidvsD5"
load("../RNAseq/Ranalysis/regPattern.rdata") # "res"           "dominance.AD1" "dominance.F1"
# restrict
res=res[as.character(ogQ$D5),]
save(genomes,tags,dirs,txdbs,ogExp,res,file="nearGenic.rdata")

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

####################################
## differential expression A vs D, using qnormed across A2/D5/ ##
####################################

## try again above using qnormed bws of A2, D5, F1.At, F1.Dt


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
plotBWsWindows(bws,list(At=grl.A$DoA,Dt=grl.D$DoA), compare="paired",paste0("At<Dt: ",length(DoA)))
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
plotBWsWindows(bws,list(At=grl.A$DoA,Dt=grl.D$DoA), compare="paired",paste0("A2<D5: ",length(DoA)))
plotBWsWindows(bws,list(At=grl.A$eqAD,Dt=grl.D$eqAD), compare="paired",paste0("A2=D5: ",length(eqAD)))
dev.off()

###################
## Hybridization ##
###################


--------------------start
# mappiable ref region, faster with bg


system("mkdir nearGenic/hybridization")
# beds
bedA = list.files("mappingA_new",pattern="A6.*n_q20.genic.bed",full.name=T)
bedD = list.files("mappingD",pattern="Dc.*_q20.genic.bed",full.name=T)
bedF = list.files("mappingF",pattern="Fc.*_q20.genic.bed",full.name=T)
bedFA = gsub("mappingF/Fc","nearGenic/hybridization/FAc",bedF)
bedFD = gsub("mappingF/Fc","nearGenic/hybridization/FDc",bedF)
cmdRun( paste0("grep '^A' ",bedF," | sed 's/^A/Chr/g' > ",bedFA ) )
cmdRun( paste0("grep '^D' ",bedF," | sed 's/^D/Chr/g' > ",bedFD ) )
beds = c(bedA,bedD,bedFA,bedFD)
# bgs
bgsA = list.files("mappingA_new",pattern="A6.*n_q20.genic_chr.size_w20.genic.bg",full.name=T)
bgsD = list.files("mappingD",pattern="Dc.*_q20.genic_chr.size_w20.genic.bg",full.name=T)
bgsF = list.files("mappingF",pattern="Fc.*_q20.genic_chr.size_w20.genic.bg",full.name=T)
bgsFA = gsub("mappingF/Fc","nearGenic/hybridization/FAc",bgsF)
bgsFD = gsub("mappingF/Fc","nearGenic/hybridization/FDc",bgsF)
cmdRun( paste0("grep '^A' ",bgsF," | sed 's/^A/Chr/g' > ",bgsFA ) )
cmdRun( paste0("grep '^D' ",bgsF," | sed 's/^D/Chr/g' > ",bgsFD ) )
bgs = c(bgsA,bgsD,bgsFA,bgsFD)
# get factor for depth
list.files("nearGenic/hybridization")
readN = filelines(beds, threads = getOption("threads", 1L))
effSize = sapply(bgs,getEffectiveRefSize.bg)
print(data.frame(readN, effSize))
# deption correction: effRefSize/totalReadsreadN
sbgs=gsub(".*/","",gsub(".bg$",".depth.bg",bgs))
cmd = paste0("awk '{print $1,$2  ,$3 ,$4* ",effSize/readN," }' OFS='\t' ",bgs," > nearGenic/hybridization/",sbgs)
cmdRun(cmd,threads = getOption("threads", 1L))
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/nearGenic/hybridization")
# quantiel normalize
qsbgs=bgQuantileNorm(sbgs)
print(qsbgs)
# examine SD, see what correction and qnorm did
print(data.frame(sapply(paste0("../../",bgs),getSD), sapply(sbgs,getSD), sapply(qsbgs,getSD) ))
# now prepare normlized w20 coverge files for L-H
l=grep("L",qsbgs,value=TRUE)
h=grep("H",qsbgs,value=TRUE)
# check to see if light and heavy are paired properly
print(data.frame(l,h))
# calculate difference between light and heavy
dbgs=bgOps(l,"difference",h,pattern="L",replacement="D")
print(dbgs)
# convert bedGraph files to bigWig
dbgsA = grep("^FA|^A",c(l,h,dbgs), value=TRUE)
dbgsD = grep("^FD|^D",c(l,h,dbgs), value=TRUE)
abwsA=bedGraphToBigWig(dbgsA,"../../mappingA_new/chr.size.txt")
abwsD=bedGraphToBigWig(dbgsD,"../../mappingD/chr.size.txt")
#or
abwsA=list.files(pattern="^(A|FA).*bw")
abwsD=list.files(pattern="^(D|FD).*bw")
names(abwsA)=gsub("_.*","",abwsA)
names(abwsD)=gsub("_.*","",abwsD)
grlA= prepWindows(list(all=rownames(res)),"A2")
grlD= prepWindows(list(all=rownames(res)),"D5")
# do A2 and At line up
plotBWsWindows(abwsA, grlA, compare="bw")
plotBWsWindows(abwsD, grlD, compare="bw") # expected
dev.off()


# pick only Dfference
bwA ="A6Dn_q20.genic_chr.size_w20.genic.depth_qnorm.bw"
bwAt="FAcD_q20.genic_chr.size_w20.genic.depth_qnorm.bw"
bwD ="DcD_q20.genic_chr.size_w20.genic.depth_qnorm.bw"
bwDt="FDcD_q20.genic_chr.size_w20.genic.depth_qnorm.bw"

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
plotBWsWindows(bws,list(At=grl.A$DoA,Dt=grl.D$DoA), compare="paired",paste0("At<Dt: ",length(DoA)))
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
plotBWsWindows(bws,list(At=grl.A$DoA,Dt=grl.D$DoA), compare="paired",paste0("A2<D5: ",length(DoA)))
plotBWsWindows(bws,list(At=grl.A$eqAD,Dt=grl.D$eqAD), compare="paired",paste0("A2=D5: ",length(eqAD)))
dev.off()

# A2 vs At, D5 vs Dt

bwsA=c(bwA,bwAt)
bwsD=c(bwD,bwDt)
names(bwsA)=c("A2","At")
names(bwsD)=c("D5","Dt")

geneL = split(rownames(res),res$category)
grl.A=prepWindows(geneL,"A2")
grl.D=prepWindows(geneL,"D5")
pdf("plotHybridization.ct.pdf")
plotBWsWindows(bwsA, grl.A, compare="bw")
plotBWsWindows(bwsD, grl.D, compare="bw")
dev.off()

geneL=split(rownames(res), res$dominance.F1)
geneL$all = rownames(res)
grl.A=prepWindows(geneL,"A2")
grl.D=prepWindows(geneL,"D5")
pdf("plotHybridization.F1dominance.pdf")
plotBWsWindows(bwsA, grl.A, compare="bw")
plotBWsWindows(bwsD, grl.D, compare="bw")
dev.off()

geneL=split(rownames(res), res$Hr.reg)
grl.A=prepWindows(geneL,"A2")
grl.D=prepWindows(geneL,"D5")
pdf("plotHybridization.Hr.pdf")
plotBWsWindows(bwsA, grl.A, compare="bw")
plotBWsWindows(bwsD, grl.D, compare="bw")
dev.off()


# genes changing biased direction
length(which(res$A*res$B<0 & !is.na(res$A*res$B) & res$Hr.reg!="Hr=0"))
# 155

## generate log2ratio At/A2, Dt/D5
cd ../nearGenic/hybridization/
bigwigCompare -b1 FDcD_q20.genic_chr.size_w20.genic.depth_qnorm.bw -b2 DcD_q20.genic_chr.size_w20.genic.depth_qnorm.bw --binSize 20 -o FDvsD5.log2ratio.bw
bigwigCompare -b1 FAcD_q20.genic_chr.size_w20.genic.depth_qnorm.bw -b2 A6Dn_q20.genic_chr.size_w20.genic.depth_qnorm.bw --binSize 20 -o FAvsA2.log2ratio.bw
#R script
bwDr = "FDvsD5.log2ratio.bw"
bwAr = "FAvsA2.log2ratio.bw"
pdf("plotHybridization.Hr.log2ratio.pdf")
plotBWsWindows(bwAr, grl.A, compare="bw")
plotBWsWindows(bwDr, grl.D, compare="bw")
dev.off()

# upstream 500bp area of D>0
i=2
genome = genomes[i]
dir=dirs[i]
txdb=loadDb(txdbs[i])
TSS = getTSS(txdb)
# rpkm
rpkm = get(genome)[["rpkm"]]
# ortholog expression
exp=rpkm[rownames(rpkm) %in% ogL[[genome]]$index,]
exp$id = ogL[[genome]]$id[match(rownames(exp), ogL[[genome]]$index)]
up=getTSS(txdb, be=500, af=0, include=exp$id)
# import bws
grA=import(bwAr); seqlevels(grA) =gsub("Chr","A", seqlevels(grA))
grD=import(bwDr); seqlevels(grD) =gsub("Chr","D", seqlevels(grD))
# get overlap
ov = findOverlaps(up, grA)
mA = aggregate(subjectHits(ov),by=list(queryHits(ov)), function(x){scores=grA[x]$score; mean(scores)})
seA = aggregate(subjectHits(ov),by=list(queryHits(ov)), function(x){scores=grA[x]$score; sd(scores)/sqrt(length(scores))})
exp$mA = mA$x[match(exp$id,up$gene_id[unique(queryHits(ov))])]
exp$seA = seA$x[match(exp$id,up$gene_id[unique(queryHits(ov))])]
ov = findOverlaps(up, grD)
mD = aggregate(subjectHits(ov),by=list(queryHits(ov)), function(x){scores=grD[x]$score; mean(scores)})
seD = aggregate(subjectHits(ov),by=list(queryHits(ov)), function(x){scores=grD[x]$score; sd(scores)/sqrt(length(scores))})
exp$mD = mD$x[match(exp$id,up$gene_id[unique(queryHits(ov))])]
exp$seD = seD$x[match(exp$id,up$gene_id[unique(queryHits(ov))])]
#
res$mA = exp[paste0(rownames(res),".A"),"mA"]
res$mD = exp[paste0(rownames(res),".D"),"mD"]
res$seA = exp[paste0(rownames(res),".A"),"seA"]
res$seD = exp[paste0(rownames(res),".D"),"seD"]
save(res,file="log2ratio.bw.rdata")

ggplot(res, aes(x=mD, y=mA,col=Hr.reg)) + geom_point()


boxplot(res[,c("mA","mD")])
library(ggplot2)
# Basic box plot
res=as.data.frame(res)
res1=res[,1:27];res1$genome="At/A2"; res1$m = res$mA; res1$se=res$seA
res2=res[,1:27];res2$genome="Dt/D5"; res2$m = res$mD; res2$se=res$seD
dat <- rbind(res1,res2)
di=dat$cisNtrans
di[dat$A>0&!is.na(dat$A)&di=="A!=0"]="A2>D5"
di[dat$A<0&!is.na(dat$A)&di=="A!=0"]="A2<D5"

pdf("plotDNAratioByExp.pdf")
p <- ggplot(dat, aes(x=genome,y=m, col=genome))+ geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=T)+ theme_minimal() + labs(title="Hr",x="Expression pattern", y = "DNS log2ratio")
print(p)
p <- ggplot(dat, aes(x=Hr.reg, y=m,col=genome))+ geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=T)+ theme_minimal() + labs(title="Hr",x="Expression pattern", y = "DNS log2ratio")
print(p)
p <- ggplot(dat, aes(x=Hr.reg, y=m,col=genome))+ geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=T)+ theme_minimal() + labs(title="Hr, subseted by parental state",x="Expression pattern", y = "DNS log2ratio")+ facet_grid(. ~ di) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
print(p)
p <- ggplot(dat, aes(x=genome, y=m,col=di))+ geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=T)+ theme_minimal() + labs(title="Hr, subseted by parental state",x="Expression pattern", y = "DNS log2ratio")+ facet_grid(. ~ Hr.reg) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
> print(p)
p <- ggplot(dat, aes(x=category, y=m,col=genome))+ geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=T)+ theme_minimal() + labs(title="Cis-trans regulation",x="Expression pattern", y = "DNS log2ratio") + theme(axis.text.x = element_text(angle = 60, hjust = 1))
print(p)
p <- ggplot(dat, aes(x=dominance.F1, y=m,col=genome))+ geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=T)+ theme_minimal() + labs(title="Expression additivity",x="Expression pattern", y = "DNS log2ratio") + theme(axis.text.x = element_text(angle = 60, hjust = 1))
print(p)
dev.off()
fit<- lm(m~Hr.reg+di+genome,data=dat)
summary(a1<-aov(fit))
# Df Sum Sq Mean Sq F value   Pr(>F)
# Hr.reg          2      1    0.48   0.391    0.676
# di              2     30   15.03  12.193 5.08e-06 ***
# genome          1     93   92.85  75.305  < 2e-16 ***
# Residuals   40691  50171    1.23
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 27 observations deleted due to missingness
TukeyHSD(x=a1,  conf.level=0.95)

fit<- lm(m~Hr.reg:di+genome,data=dat)
summary(a2<-aov(fit))
TukeyHSD(x=a2, conf.level=0.95)
