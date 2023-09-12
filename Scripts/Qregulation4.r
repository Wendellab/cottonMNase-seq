# This script aims to examine duplicated expression patterns in association with:
# 1. promoter accessibility by aggregation profiles
# 2. identified ACRs in promoter regions
# 3. DA in promoter regions
# By Guanjing Hu


# module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core
# module load r/3.5.0-py2-ufvuwmm
options(scipen=999)
library(GenomicRanges)
#library(RColorBrewer)
#library(ggbio)
library(GenomicFeatures)
library(genomation)

setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/")

# incorrect for - strand if be!=af
getTSS0=function(txdb, include=NULL, exclude.pattern=NULL, be=1000, af=1000){
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
getTSS=function(txdb, include=NULL, exclude.pattern=NULL, be=1000, af=1000){
    tss =  genes(txdb)
    if(!is.null(exclude.pattern)) {tss = tss[grep(exclude.pattern,tss$gene_id,invert=T)]}
    if(!is.null(include)) {tss = tss[tss$gene_id %in% include]}
    # +
    end(tss[strand(tss)=="+",])  =start(tss[strand(tss)=="+",])+af
    start(tss[strand(tss)=="+",])  =start(tss[strand(tss)=="+",])-be
    # -
    start(tss[strand(tss)=="-",])=end(tss[strand(tss)=="-",])-af
    end(tss[strand(tss)=="-",])=end(tss[strand(tss)=="-",])+be
    # remove duplicated TSSes ie alternative transcripts
    # this keeps the first instance and removes duplicates
    tss=tss[!duplicated(tss),]
    seqlevels(tss)=unique(as.character(seqnames(tss)))
    return(tss)
}
# extract target data over pre-defined windowns, multiple file
plotMetaPair = function(bw, windowA, windowD, overlay=TRUE, main="Meta-profile", ylim=NULL)
{
    sml = list(A = ScoreMatrix(target = bw, windows = windowA, type="bigWig",strand.aware=TRUE), Dt = ScoreMatrix(target = bw, windows = windowD, type="bigWig",strand.aware=TRUE) )
    names(sml) =c("A","D")
    sml= new("ScoreMatrixList",sml)
    if(overlay){
        plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="differential nuclease sensitivity score (tpm)",dispersion = "se", xlab="bases around TSS", main=main, ylim=ylim)
        abline(v=0,lty=2)
    }else{
        par(mfrow=c(1,length(sml)))
        plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="differential nuclease sensitivity score (tpm)",dispersion = "se", xlab="bases around TSS", main=main, overlay=FALSE, ylim=ylim)
    }
}

plotBWsWindows = function(bws, windowList, compare=c("bw","window","paired"), xlab="bases around TSS", main="", ylim = NULL, line.col=NULL, dispersion.col = NULL){
    nbw=length(bws)
    if(compare=="window"){
        for(b in 1:nbw)
        {
            sml=list()
            for(w in names(windowList)){
                sml[[w]] =  ScoreMatrix(target = bws[b], windows = windowList[[w]], type="bigWig",strand.aware=TRUE)
            }
            sml= new("ScoreMatrixList",sml)
            plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="",dispersion = "se", xlab=xlab, main=paste(main, names(bws)[b]),line.col=line.col, dispersion.col = dispersion.col, ylim=ylim)
            abline(v=0,lty=2)
        }
    }
    if(compare=="bw"){
        for(w in names(windowList)){
            sml=ScoreMatrixList(target = bws, windowList[[w]], type="bigWig",strand.aware=TRUE)
            names(sml) = names(bws)
            plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="",dispersion = "se", xlab=xlab, main=paste(main,w,length(windowList[[w]])), ,line.col=line.col, dispersion.col = dispersion.col, ylim=ylim)
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
            plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(bws), ylab="",dispersion = "se", xlab=xlab, main=main,,line.col=line.col, dispersion.col = dispersion.col, ylim=ylim)
            abline(v=0,lty=2)
        }
    }
}

makeTransparent = function(..., alpha=0.4) {
    if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
    
    alpha = floor(255*alpha)
    newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
    
    .makeTransparent = function(col, alpha) {
        rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
    }
    
    newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
    
    return(newColor)
}

ogQ<-read.table("orthohomoeolog052421.txt", sep="\t", header=FALSE,stringsAsFactor = FALSE)
ogQ$V3<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",ogQ$V3)))
ogQ$V4<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",ogQ$V4)))



#################
## Compare ref ##
#################
load('~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/4.RNAseqAnalysis/IndivGNT2021rRNA/regPattern.rdata')
rI<-res  # Gorai.002G001200
adI<-dominance.AD1
fI<-dominance.F1
load('~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/4.RNAseqAnalysis/AD1GNT2021rRNA/regPattern.rdata')
rM<-res  # Gorai.002G001200.1
adM<-dominance.AD1
fM<-dominance.F1
load('~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/4.RNAseqAnalysis/Dref2021/regPattern.rdata')
rD<-res  # Gorai.002G001200.1
adD<-dominance.AD1
fD<-dominance.F1

# table(paste0(rownames(rI),".1") %in% rownames(rM))
rownames(rI)<-paste0(rownames(rI),".1")
rownames(adI)<-paste0(rownames(adI),".1")
rownames(fI)<-paste0(rownames(fI),".1")

# annotate direction for dominance
addDirection0<-function(D){
    
    D$TvsA.reg = ifelse(D$TvsA.reg=="T=A","T=A",ifelse(D$TvsA>0,"T>A","T<A"))
    D$TvsD.reg = ifelse(D$TvsD.reg=="T=D","T=D",ifelse(D$TvsD>0,"T>D","T<D"))
    return(D)
}
adI=addDirection0(adI)
fI=addDirection0(fI)
adM=addDirection0(adM)
fM=addDirection0(fM)
adD=addDirection0(adD)
fD=addDirection0(fD)

# annotate direction
addDirection<-function(res){
    res$A0 = ifelse(res$cisNtrans=="A=0","A=0",ifelse(res$A>0,"A>0","A<0"))
    res$B0 = ifelse(res$cis=="B=0","B=0",ifelse(res$B>0,"B>0","B<0"))
    res$Bp0 = ifelse(is.na(res$Bp.padj),"Bp=0",ifelse(res$Bp.padj>0.05,"Bp=0",ifelse(res$Bp>0,"Bp>0","Bp<0")))
    res$Hr0 = ifelse(res$Hr.reg=="Hr=0","Hr=0",ifelse(res$Hr>0,"Hr>0","Hr<0"))
    res$Pr0 = ifelse(res$Pr.reg=="Pr=0","Pr=0",ifelse(res$Pr>0,"Pr>0","Pr<0"))
    res$Wr0 = ifelse(res$Wr.reg=="Wr=0","Wr=0",ifelse(res$Wr>0,"Wr>0","Wr<0"))
    return(res)
}
rI=addDirection(rI)
rM=addDirection(rM)
rD=addDirection(rD)

con<-data.frame(D=ogQ$V5)
rownames(con)=ogQ$V5
# F1 dominance
to=data.frame(Indiv=fI[ogQ$V5,"category"], AD1=fM[ogQ$V5,"category"], D5=fD[ogQ$V5,"category"])
to$cat = apply(to,1,function(x)ifelse(length(unique(x))==1,x[1],"-"))
apply(to,2,table)
con$F1.dominance=to$cat
to=data.frame(Indiv=fI[ogQ$V5,"TvsA.reg"], AD1=fM[ogQ$V5,"TvsA.reg"], D5=fD[ogQ$V5,"TvsA.reg"])
to$cat = apply(to,1,function(x)ifelse(length(unique(x))==1,x[1],"-"))
apply(to,2,table)
con$F1.TvsA=to$cat
to=data.frame(Indiv=fI[ogQ$V5,"TvsD.reg"], AD1=fM[ogQ$V5,"TvsD.reg"], D5=fD[ogQ$V5,"TvsD.reg"])
to$cat = apply(to,1,function(x)ifelse(length(unique(x))==1,x[1],"-"))
apply(to,2,table)
con$F1.TvsD=to$cat

# AD1 dominance
to=data.frame(Indiv=adI[ogQ$V5,"category"], AD1=adM[ogQ$V5,"category"], D5=adD[ogQ$V5,"category"])
to$cat = apply(to,1,function(x)ifelse(length(unique(x))==1,x[1],"-"))
apply(to,2,table)
con$AD1.dominance=to$cat
to=data.frame(Indiv=adI[ogQ$V5,"TvsA.reg"], AD1=adM[ogQ$V5,"TvsA.reg"], D5=adD[ogQ$V5,"TvsA.reg"])
to$cat = apply(to,1,function(x)ifelse(length(unique(x))==1,x[1],"-"))
apply(to,2,table)
con$AD1.TvsA=to$cat
to=data.frame(Indiv=adI[ogQ$V5,"TvsD.reg"], AD1=adM[ogQ$V5,"TvsD.reg"], D5=adD[ogQ$V5,"TvsD.reg"])
to$cat = apply(to,1,function(x)ifelse(length(unique(x))==1,x[1],"-"))
apply(to,2,table)
con$AD1.TvsD=to$cat


#A
to=data.frame(Indiv=rI[ogQ$V5,"A0"], AD1=rM[ogQ$V5,"A0"], D5=rD[ogQ$V5,"A0"])
to$cat = apply(to,1,function(x)ifelse(length(unique(x))==1,x[1],"-"))
apply(to,2,table)
con$A=to$cat

#B
to=data.frame(Indiv=rI[ogQ$V5,"B0"], AD1=rM[ogQ$V5,"B0"], D5=rD[ogQ$V5,"B0"])
to$cat = apply(to,1,function(x)ifelse(length(unique(x))==1,x[1],"-"))
apply(to,2,table)
con$B=to$cat

#Bp
to=data.frame(Indiv=rI[ogQ$V5,"Bp0"], AD1=rM[ogQ$V5,"Bp0"], D5=rD[ogQ$V5,"Bp0"])
to$cat = apply(to,1,function(x)ifelse(length(unique(x))==1,x[1],"-"))
apply(to,2,table)
con$Bp=to$cat

#Hr
to=data.frame(Indiv=rI[ogQ$V5,"Hr0"], AD1=rM[ogQ$V5,"Hr0"], D5=rD[ogQ$V5,"Hr0"])
to$cat = apply(to,1,function(x)ifelse(length(unique(x))==1,x[1],"-"))
apply(to,2,table)
con$Hr=to$cat

#Pr
to=data.frame(Indiv=rI[ogQ$V5,"Pr0"], AD1=rM[ogQ$V5,"Pr0"], D5=rD[ogQ$V5,"Pr0"])
to$cat = apply(to,1,function(x)ifelse(length(unique(x))==1,x[1],"-"))
apply(to,2,table)
con$Pr=to$cat

#Wr
to=data.frame(Indiv=rI[ogQ$V5,"Wr0"], AD1=rM[ogQ$V5,"Wr0"], D5=rD[ogQ$V5,"Wr0"])
to$cat = apply(to,1,function(x)ifelse(length(unique(x))==1,x[1],"-"))
apply(to,2,table)
con$Wr=to$cat

#CisTrans
to=data.frame(Indiv=rI[ogQ$V5,"category"], AD1=rM[ogQ$V5,"category"], D5=rD[ogQ$V5,"category"])
to$cat = apply(to,1,function(x)ifelse(length(unique(x))==1,x[1],"-"))
apply(to,2,table)
con$cistrans=to$cat

write.table(cbind(con, ogQ),file="conserved.regPattern.txt",sep="\t")

con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")

########
## Bp ##
########
con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")

txdb = loadDb("refGenomes/txdb.AD1utx.sqlite")
TSS.all = getTSS(txdb, exclude.pattern="Gohir.1Z")
TSS.A = TSS.all[grep("A",seqnames(TSS.all))] #36116
TSS.D = TSS.all[grep("D",seqnames(TSS.all))] #38784
# OG 22889 groups
TSS.Aog = TSS.all[ogQ$V3]
TSS.Dog = TSS.all[ogQ$V4]
# Not in OG
TSS.Aother = TSS.A[!names(TSS.A) %in% ogQ$V3] #13227
TSS.Dother = TSS.D[!names(TSS.D) %in% ogQ$V4] #15895


# bw from DNS,mapped AD1
bw = "mappingM_UTX/McD_q20_chr.size_w20_qnorm.bw"
# bw="Qregulation4/McD_q20_chr.size_w20_qnorm_qnorm.bw"
# bw = "Qregulation_diffACRs/AD1ref/mapping/McD_q20.full.bw"
## DE based on conserved
Abias.A= as.character(con[con$Bp=="Bp>0","V3"])
Abias.D= as.character(con[con$Bp=="Bp>0","V4"])
Dbias.A= as.character(con[con$Bp=="Bp<0","V3"])
Dbias.D= as.character(con[con$Bp=="Bp<0","V4"])
unbias.A= as.character(con[con$Bp=="Bp=0","V3"])
unbias.D= as.character(con[con$Bp=="Bp=0","V4"])
TSS.Abias.A =  TSS.A[Abias.A] #2122
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A] #2400
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
pdf("Qregulation4/plotDNS.Bp.con2AD1.pdf")
plotMetaPair(bw, windowA=TSS.A, windowD=TSS.D, main=paste("All: A vs D,",length(TSS.A),length(TSS.D)))
plotMetaPair(bw, windowA=TSS.Aog, windowD=TSS.Dog, main=paste("OG: A vs D,",length(TSS.Aog),length(TSS.Dog)))
plotMetaPair(bw, windowA=TSS.Aother, windowD=TSS.Dother, main=paste("Other: A vs D,",length(TSS.Aother),length(TSS.Dother)))
plotMetaPair(bw, windowA=TSS.Abias.A, windowD=TSS.Abias.D, main=paste("A bias expression: AT > Dt,",length(TSS.Abias.A)))
plotMetaPair(bw, windowA=TSS.Dbias.A, windowD=TSS.Dbias.D, main=paste("D bias expression: At < Dt,",length(TSS.Dbias.A)))
plotMetaPair(bw, windowA=TSS.unbias.A, windowD=TSS.unbias.D, main=paste("Unbias expression: At = Dt,",length(TSS.unbias.A)))
dev.off()




## DE based AD1 reference: same results from IndivGNRr and AD1GNTr
## --------
#load("RNAseq/RanalysisIndivGNTr/diffExprPattern.rdata")->l;l
load("RNAseq/RanalysisAD1GNTr/diffExprPattern.rdata")->l;l
Bp=data.frame(Bp)
Bp$cat<-"Bp=0"
Bp$cat[Bp$padj<0.05&Bp$log2FoldChange>0&!is.na(Bp$padj)]<-"Bp>0"
Bp$cat[Bp$padj<0.05&Bp$log2FoldChange<0&!is.na(Bp$padj)]<-"Bp<0"
table(Bp$cat)
BpP<-Bp
Abias.A= ogQ$V3[match( rownames(Bp[Bp$cat=="Bp>0",]), ogQ$V5)]
Abias.D= ogQ$V4[match( rownames(Bp[Bp$cat=="Bp>0",]), ogQ$V5)]
Dbias.A= ogQ$V3[match( rownames(Bp[Bp$cat=="Bp<0",]), ogQ$V5)]
Dbias.D= ogQ$V4[match( rownames(Bp[Bp$cat=="Bp<0",]), ogQ$V5)]
unbias.A= ogQ$V3[match( rownames(Bp[Bp$cat=="Bp=0",]), ogQ$V5)]
unbias.D= ogQ$V4[match( rownames(Bp[Bp$cat=="Bp=0",]), ogQ$V5)]
TSS.Abias.A =  TSS.A[Abias.A] #2122
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A] #2400
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
pdf("Qregulation4/plotDNSdt.Bp.AD1ref.pdf")
plotMetaPair(bwA, windowA=TSS.A, windowD=TSS.D, main=paste("All: A vs D,",length(TSS.A),length(TSS.D)))
plotMetaPair(bwA, windowA=TSS.Aog, windowD=TSS.Dog, main=paste("OG: A vs D,",length(TSS.Aog),length(TSS.Dog)))
plotMetaPair(bwA, windowA=TSS.Aother, windowD=TSS.Dother, main=paste("Other: A vs D,",length(TSS.Aother),length(TSS.Dother)))
plotMetaPair(bwA, windowA=TSS.Abias.A, windowD=TSS.Abias.D, main=paste("A bias expression: AT > Dt,",length(TSS.Abias.A)))
plotMetaPair(bwA, windowA=TSS.Dbias.A, windowD=TSS.Dbias.D, main=paste("D bias expression: At < Dt,",length(TSS.Dbias.A)))
plotMetaPair(bwA, windowA=TSS.unbias.A, windowD=TSS.unbias.D, main=paste("Unbias expression: At = Dt,",length(TSS.unbias.A)))
dev.off()
pdf("Qregulation4/plotDNS.Bp.AD1ref.pdf")
plotMetaPair(bw, windowA=TSS.A, windowD=TSS.D, main=paste("All: A vs D,",length(TSS.A),length(TSS.D)))
plotMetaPair(bw, windowA=TSS.Aog, windowD=TSS.Dog, main=paste("OG: A vs D,",length(TSS.Aog),length(TSS.Dog)))
plotMetaPair(bw, windowA=TSS.Aother, windowD=TSS.Dother, main=paste("Other: A vs D,",length(TSS.Aother),length(TSS.Dother)))
plotMetaPair(bw, windowA=TSS.Abias.A, windowD=TSS.Abias.D, main=paste("A bias expression: AT > Dt,",length(TSS.Abias.A)))
plotMetaPair(bw, windowA=TSS.Dbias.A, windowD=TSS.Dbias.D, main=paste("D bias expression: At < Dt,",length(TSS.Dbias.A)))
plotMetaPair(bw, windowA=TSS.unbias.A, windowD=TSS.unbias.D, main=paste("Unbias expression: At = Dt,",length(TSS.unbias.A)))
dev.off()




load("RNAseq/RanalysisDref/diffExprPattern.rdata")->l;l
Bp=data.frame(Bp)
Bp$cat<-"Bp=0"
Bp$cat[Bp$padj<0.05&Bp$log2FoldChange>0&!is.na(Bp$padj)]<-"Bp>0"
Bp$cat[Bp$padj<0.05&Bp$log2FoldChange<0&!is.na(Bp$padj)]<-"Bp<0"
table(Bp$cat)
Abias.A= ogQ$V3[match( rownames(Bp[Bp$cat=="Bp>0",]), ogQ$V5)]
Abias.D= ogQ$V4[match( rownames(Bp[Bp$cat=="Bp>0",]), ogQ$V5)]
Dbias.A= ogQ$V3[match( rownames(Bp[Bp$cat=="Bp<0",]), ogQ$V5)]
Dbias.D= ogQ$V4[match( rownames(Bp[Bp$cat=="Bp<0",]), ogQ$V5)]
unbias.A= ogQ$V3[match( rownames(Bp[Bp$cat=="Bp=0",]), ogQ$V5)]
unbias.D= ogQ$V4[match( rownames(Bp[Bp$cat=="Bp=0",]), ogQ$V5)]
TSS.Abias.A =  TSS.all[na.omit(Abias.A)] #1955
TSS.Abias.D =  TSS.all[na.omit(Abias.D)]
TSS.Dbias.A =  TSS.all[na.omit(Dbias.A)] #1946
TSS.Dbias.D =  TSS.all[na.omit(Dbias.D)]
TSS.unbias.A =  TSS.all[na.omit(unbias.A)] #15783
TSS.unbias.D =  TSS.all[na.omit(unbias.D)]
pdf("Qregulation4/plotDNS.Bp.Dref.pdf")
plotMetaPair(bw, windowA=TSS.A, windowD=TSS.D, main="All: A vs D")
plotMetaPair(bw, windowA=TSS.Aog, windowD=TSS.Dog, main="OG: A vs D")
plotMetaPair(bw, windowA=TSS.Aother, windowD=TSS.Dother, main="Other: A vs D")
plotMetaPair(bw, windowA=TSS.Abias.A, windowD=TSS.Abias.D, main="A bias expression: AT > Dt")
plotMetaPair(bw, windowA=TSS.Dbias.A, windowD=TSS.Dbias.D, main="D bias expression: At < Dt")
plotMetaPair(bw, windowA=TSS.unbias.A, windowD=TSS.unbias.D, main="Unbias expression: At = Dt")
dev.off()

# compare results between AD1ref and A2+D5
Bp$catAD<-BpP[rownames(Bp),"cat"]
Bp$catAD[is.na(Bp$catAD)] <-"other"
table(Bp$cat,Bp$catAD)
#       Bp<0  Bp=0  Bp>0 other
# Bp<0  1682   252    12  1324
# Bp=0   704 17724   560 11917
# Bp>0    14   391  1550  1093

#######
## B ##
#######

con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")

# bw from aggregation, mapped to AD1
bw = "Qregulation_diffACRs/AD1ref/mapping/FcD_q20.full.bw"
## DE based on conserved
Abias.A= as.character(con[con$B=="B>0","V3"])
Abias.D= as.character(con[con$B=="B>0","V4"])
Dbias.A= as.character(con[con$B=="B<0","V3"])
Dbias.D= as.character(con[con$B=="B<0","V4"])
unbias.A= as.character(con[con$B=="B=0","V3"])
unbias.D= as.character(con[con$B=="B=0","V4"])
TSS.Abias.A =  TSS.A[Abias.A]
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A]
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
pdf("Qregulation4/plotDNS.B.con2AD1.pdf")
plotMetaPair(bw, windowA=TSS.A, windowD=TSS.D, main=paste("All: A vs D,",length(TSS.A),length(TSS.D)))
plotMetaPair(bw, windowA=TSS.Aog, windowD=TSS.Dog, main=paste("OG: A vs D,",length(TSS.Aog),length(TSS.Dog)))
plotMetaPair(bw, windowA=TSS.Aother, windowD=TSS.Dother, main=paste("Other: A vs D,",length(TSS.Aother),length(TSS.Dother)))
plotMetaPair(bw, windowA=TSS.Abias.A, windowD=TSS.Abias.D, main=paste("A bias expression: AT > Dt,",length(TSS.Abias.A)))
plotMetaPair(bw, windowA=TSS.Dbias.A, windowD=TSS.Dbias.D, main=paste("D bias expression: At < Dt,",length(TSS.Dbias.A)))
plotMetaPair(bw, windowA=TSS.unbias.A, windowD=TSS.unbias.D, main=paste("Unbias expression: At = Dt,",length(TSS.unbias.A)))
dev.off()

txdb = loadDb("refGenomes/txdb.F2020.sqlite")
TSS.all = getTSS(txdb, exclude.pattern="Gorai.N")
TSS.A = TSS.all[grepl("Gar",names(TSS.all))&!grepl("Contig",seqnames(TSS.all))] #41739
TSS.D = TSS.all[grep("Gorai",names(TSS.all))] #37218
#
con$D5<-gsub("[.]1$","",con$V5)
con$A2<-paste0(con$V2, ".gene")
# OG 22889 groups
TSS.Aog = TSS.all[con$A2]
TSS.Dog = TSS.all[con$D5]
# Not in OG
TSS.Aother = TSS.A[!names(TSS.A) %in% con$A2] #18850
TSS.Dother = TSS.D[!names(TSS.D) %in% con$D5] #14329
# bw from DNS, mapped to A2+D5
bw = "mappingF2020/FcD_q20_chr.size_w20_qnorm.bw"
## DE based on conserved
Abias.A= as.character(con[con$B=="B>0","A2"])
Abias.D= as.character(con[con$B=="B>0","D5"])
Dbias.A= as.character(con[con$B=="B<0","A2"])
Dbias.D= as.character(con[con$B=="B<0","D5"])
unbias.A= as.character(con[con$B=="B=0","A2"])
unbias.D= as.character(con[con$B=="B=0","D5"])
TSS.Abias.A =  TSS.A[Abias.A] #2122
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A] #2400
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
pdf("Qregulation4/plotDNS.B.con2F1.pdf")
plotMetaPair(bw, windowA=TSS.A, windowD=TSS.D, main=paste("All: A vs D,",length(TSS.A),length(TSS.D)))
plotMetaPair(bw, windowA=TSS.Aog, windowD=TSS.Dog, main=paste("OG: A vs D,",length(TSS.Aog),length(TSS.Dog)))
plotMetaPair(bw, windowA=TSS.Aother, windowD=TSS.Dother, main=paste("Other: A vs D,",length(TSS.Aother),length(TSS.Dother)))
plotMetaPair(bw, windowA=TSS.Abias.A, windowD=TSS.Abias.D, main=paste("A bias expression: AT > Dt,",length(TSS.Abias.A)))
plotMetaPair(bw, windowA=TSS.Dbias.A, windowD=TSS.Dbias.D, main=paste("D bias expression: At < Dt,",length(TSS.Dbias.A)))
plotMetaPair(bw, windowA=TSS.unbias.A, windowD=TSS.unbias.D, main=paste("Unbias expression: At = Dt,",length(TSS.unbias.A)))
dev.off()
# bw from DNS + split subgenome + qnorm before iseg, mapped to A2 and D5
# this may be more comparable to diploids results, because relative A vs D is affected by mapping quality in F1
bws = c("Qregulation4/FcDa_q20_chr.size_w20_qnorm_qnorm.bw","Qregulation4/FcDd_q20_chr.size_w20_qnorm_qnorm.bw")
names(bws)=c("A2","D5")
txdb = loadDb("refGenomes/txdb.A2WHU.sqlite")
TSS = getTSS(txdb, exclude.pattern="Gorai.N")
TSS.A = TSS[grepl("Gar",names(TSS))&!grepl("Contig",seqnames(TSS))] #41739
txdb = loadDb("refGenomes/txdb.D5.sqlite")
TSS = getTSS(txdb, exclude.pattern="Gorai.N")
TSS.D = TSS[grep("Gorai",names(TSS))] #37218
## OG 22889 groups
TSS.Aog = TSS.A[con$A2]
TSS.Dog = TSS.D[con$D5]
# Not in OG
TSS.Aother = TSS.A[!names(TSS.A) %in% con$A2] #18850
TSS.Dother = TSS.D[!names(TSS.D) %in% con$D5] #14329
## DE based on conserved
Abias.A= as.character(con[con$B=="B>0","A2"])
Abias.D= as.character(con[con$B=="B>0","D5"])
Dbias.A= as.character(con[con$B=="B<0","A2"])
Dbias.D= as.character(con[con$B=="B<0","D5"])
unbias.A= as.character(con[con$B=="B=0","A2"])
unbias.D= as.character(con[con$B=="B=0","D5"])
TSS.Abias.A =  TSS.A[Abias.A] #2122
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A] #2400
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
pdf("Qregulation4/plotDNS.B.con2A+D.pdf")
plotBWsWindows(bws,list(A2=TSS.A, D5=TSS.D), compare="paired", main=paste("All: A vs D,",length(TSS.A),length(TSS.D)) )
plotBWsWindows(bws,list(A2=TSS.Aog, D5=TSS.Dog), compare="paired", main=paste("OG: A vs D,",length(TSS.Aog),length(TSS.Dog)) )
plotBWsWindows(bws,list(A2=TSS.Aother, D5=TSS.Dother), compare="paired", main=paste("Other: A vs D,",length(TSS.Aother),length(TSS.Dother)) )
plotBWsWindows(bws, list(A2=TSS.Abias.A, D5=TSS.Abias.D), compare="paired", main=paste("A bias expression: AT > Dt,",length(TSS.Abias.A)))
plotBWsWindows(bws, list(A2=TSS.Dbias.A, D5=TSS.Dbias.D), compare="paired", main=paste("A bias expression: AT < Dt,",length(TSS.Dbias.A)))
plotBWsWindows(bws, list(A2=TSS.unbias.A, D5=TSS.unbias.D), compare="paired", main=paste("A bias expression: AT = Dt,",length(TSS.unbias.A)))
dev.off()


#######
## A ##
#######

con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")

load("Qregulation4/results.rdata")

# bw from A2 qnorm + D5 qnorm, mapped to F1
txdb = loadDb("refGenomes/txdb.F2020.sqlite")
TSS.all = getTSS(txdb, exclude.pattern="Gorai.N")
TSS.A = TSS.all[grepl("Gar",names(TSS.all))&!grepl("Contig",seqnames(TSS.all))] #41739
TSS.D = TSS.all[grep("Gorai",names(TSS.all))] #37218
#
con$D5<-gsub("[.]1$","",con$V5)
con$A2<-paste0(con$V2, ".gene")
# OG 22889 groups
TSS.Aog = TSS.all[con$A2]
TSS.Dog = TSS.all[con$D5]
# Not in OG
TSS.Aother = TSS.A[!names(TSS.A) %in% con$A2] #18850
TSS.Dother = TSS.D[!names(TSS.D) %in% con$D5] #14329
# bw from DNS, mapped to A2+D5
bw = "Qregulation4/ScD_q20_chr.size_w20_qnorm_qnorm.bw"
## DE based on conserved
Abias.A= as.character(con[con$A=="A>0","A2"])
Abias.D= as.character(con[con$A=="A>0","D5"])
Dbias.A= as.character(con[con$A=="A<0","A2"])
Dbias.D= as.character(con[con$A=="A<0","D5"])
unbias.A= as.character(con[con$A=="A=0","A2"])
unbias.D= as.character(con[con$A=="A=0","D5"])
TSS.Abias.A =  TSS.A[Abias.A] #2122
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A] #2400
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
pdf("Qregulation4/plotDNS.A.con2F1.pdf")
plotMetaPair(bw, windowA=TSS.A, windowD=TSS.D, main=paste("All: A vs D,",length(TSS.A),length(TSS.D)))
plotMetaPair(bw, windowA=TSS.Aog, windowD=TSS.Dog, main=paste("OG: A vs D,",length(TSS.Aog),length(TSS.Dog)))
plotMetaPair(bw, windowA=TSS.Aother, windowD=TSS.Dother, main=paste("Other: A vs D,",length(TSS.Aother),length(TSS.Dother)))
plotMetaPair(bw, windowA=TSS.Abias.A, windowD=TSS.Abias.D, main=paste("Parental divergence: A2 > D5,",length(TSS.Abias.A)))
plotMetaPair(bw, windowA=TSS.Dbias.A, windowD=TSS.Dbias.D, main=paste("Parental divergence: A2 < D5,",length(TSS.Dbias.A)))
plotMetaPair(bw, windowA=TSS.unbias.A, windowD=TSS.unbias.D, main=paste("Parental divergence: A2 = D5,",length(TSS.unbias.A)))
dev.off()


## parental expression divergence A2 vs D5, when A!=0, mapped to each ref
bws=c("Qregulation4/A6Dn_q20_chr.size_w20_qnorm_qnorm.bw","Qregulation4/DcD_q20_chr.size_w20_qnorm_qnorm.bw")
names(bws)=c("A2","D5")
txdb = loadDb("refGenomes/txdb.A2WHU.sqlite")
TSS = getTSS(txdb, exclude.pattern="Gorai.N")
TSS.A = TSS[grepl("Gar",names(TSS))&!grepl("Contig",seqnames(TSS))] #41739
txdb = loadDb("refGenomes/txdb.D5.sqlite")
TSS = getTSS(txdb, exclude.pattern="Gorai.N")
TSS.D = TSS[grep("Gorai",names(TSS))] #37218
## OG 22889 groups
TSS.Aog = TSS.A[con$A2]
TSS.Dog = TSS.D[con$D5]
# Not in OG
TSS.Aother = TSS.A[!names(TSS.A) %in% con$A2] #18850
TSS.Dother = TSS.D[!names(TSS.D) %in% con$D5] #14329
## DE based on conserved
Abias.A= as.character(con[con$A=="A>0","A2"])
Abias.D= as.character(con[con$A=="A>0","D5"])
Dbias.A= as.character(con[con$A=="A<0","A2"])
Dbias.D= as.character(con[con$A=="A<0","D5"])
unbias.A= as.character(con[con$A=="A=0","A2"])
unbias.D= as.character(con[con$A=="A=0","D5"])
TSS.Abias.A =  TSS.A[Abias.A] #2122
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A] #2400
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
pdf("Qregulation4/plotDNS.A.con2A+D.pdf")
plotBWsWindows(bws,list(A2=TSS.A, D5=TSS.D), compare="paired", main=paste("All: A vs D,",length(TSS.A),length(TSS.D)) )
plotBWsWindows(bws,list(A2=TSS.Aog, D5=TSS.Dog), compare="paired", main=paste("OG: A vs D,",length(TSS.Aog),length(TSS.Dog)) )
plotBWsWindows(bws,list(A2=TSS.Aother, D5=TSS.Dother), compare="paired", main=paste("Other: A vs D,",length(TSS.Aother),length(TSS.Dother)) )
plotBWsWindows(bws, list(A2=TSS.Abias.A, D5=TSS.Abias.D), compare="paired", main=paste("Parental divergence: A2 > D5,",length(TSS.Abias.A)))
plotBWsWindows(bws, list(A2=TSS.Dbias.A, D5=TSS.Dbias.D), compare="paired", main=paste("Parental divergence: A2 < D5,",length(TSS.Dbias.A)))
plotBWsWindows(bws, list(A2=TSS.unbias.A, D5=TSS.unbias.D), compare="paired", main=paste("Parental divergence: A2 = D5,",length(TSS.unbias.A)))
dev.off()


########
## Hr ##
########
con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")
con$D5<-gsub("[.]1$","",con$V5)
con$A2<-paste0(con$V2, ".gene")

txdb = loadDb("refGenomes/txdb.A2WHU.sqlite")
TSS = getTSS(txdb, exclude.pattern="Gorai.N")
TSS.A = TSS[grepl("Gar",names(TSS))&!grepl("Contig",seqnames(TSS))] #41739
txdb = loadDb("refGenomes/txdb.D5.sqlite")
TSS = getTSS(txdb, exclude.pattern="Gorai.N")
TSS.D = TSS[grep("Gorai",names(TSS))] #37218
## OG 22889 groups
TSS.Aog = TSS.A[con$A2]
TSS.Dog = TSS.D[con$D5]
# Not in OG
TSS.Aother = TSS.A[!names(TSS.A) %in% con$A2] #18850
TSS.Dother = TSS.D[!names(TSS.D) %in% con$D5] #14329


## DE based on conserved
up.A= as.character(con[con$Hr=="Hr>0","A2"])
up.D= as.character(con[con$Hr=="Hr>0","D5"])
down.A= as.character(con[con$Hr=="Hr<0","A2"])
down.D= as.character(con[con$Hr=="Hr<0","D5"])
uc.A= as.character(con[con$Hr=="Hr=0","A2"])
uc.D= as.character(con[con$Hr=="Hr=0","D5"])
TSS.up.A =  TSS.A[up.A] #2122
TSS.up.D =  TSS.D[up.D]
TSS.down.A =  TSS.A[down.A] #2400
TSS.down.D =  TSS.D[down.D]
TSS.uc.A =  TSS.A[uc.A]
TSS.uc.D =  TSS.D[uc.D]


## parental expression divergence A2 vs D5, when A!=0, mapped to each ref
#bws=c("Qregulation4/A6Dn_q20_chr.size_w20_qnorm_qnorm.bw","Qregulation4/DcD_q20_chr.size_w20_qnorm_qnorm.bw", "Qregulation4/FcDa_q20_chr.size_w20_qnorm_qnorm.bw", "Qregulation4/FcDd_q20_chr.size_w20_qnorm_qnorm.bw")
#names(bws)=c("A2","D5","At","Dt")

bws=c("Qregulation_diffACRs/Fref/bwRPGC/A6Dn_q20.full.bw", "Qregulation_diffACRs/Fref/bwRPGC/DcD_q20.full.bw", "Qregulation_diffACRs/Fref/bwRPGC/FcDa_q20.full.bw", "Qregulation_diffACRs/Fref/bwRPGC/FcDd_q20.full.bw", "Qregulation_diffACRs/Fref/bwRPGC/McDa_q20.full.bw", "Qregulation_diffACRs/Fref/bwRPGC/McDd_q20.full.bw")
names(bws)=c("A2","D5","F1.At","F1.Dt","AD1.At","AD1.Dt")
color=c("purple","brown","blue")


pdf("Qregulation4/plotDNS_RPGC.Hr.con2A+D.pdf")

plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="All ", ylim=c(0,1))
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="All ", ylim=c(0,1) )

plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.Aog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG ", ylim=c(0,1))
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.Dog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG ", ylim=c(0,1))

plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.Aother), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Other ", ylim=c(0,1))
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.Dother), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Other ", ylim=c(0,1))

plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.up.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr>0 ", ylim=c(0,1))
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.up.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr>0 ", ylim=c(0,1))

plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.down.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr<0 ", ylim=c(0,1))
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.down.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr<0 ", ylim=c(0,1))

plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.uc.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr=0 ", ylim=c(0,1))
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.uc.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr=0 ", ylim=c(0,1))

## DE based on conserved
up.A= as.character(con[con$Wr=="Wr>0","A2"])
up.D= as.character(con[con$Wr=="Wr>0","D5"])
down.A= as.character(con[con$Wr=="Wr<0","A2"])
down.D= as.character(con[con$Wr=="Wr<0","D5"])
uc.A= as.character(con[con$Wr=="Wr=0","A2"])
uc.D= as.character(con[con$Wr=="Wr=0","D5"])
TSS.up.A =  TSS.A[up.A] #2122
TSS.up.D =  TSS.D[up.D]
TSS.down.A =  TSS.A[down.A] #2400
TSS.down.D =  TSS.D[down.D]
TSS.uc.A =  TSS.A[uc.A]
TSS.uc.D =  TSS.D[uc.D]

plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.up.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr>0 ", ylim=c(0,1))
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.up.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr>0 ", ylim=c(0,1))

plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.down.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr<0 ", ylim=c(0,1))
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.down.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr<0 ", ylim=c(0,1))

plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.uc.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr=0 ", ylim=c(0,1))
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.uc.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr=0 ", ylim=c(0,1))

dev.off()


########
## Wr ##
########
con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")

txdb = loadDb("refGenomes/txdb.AD1utx.sqlite")
TSS.all = getTSS(txdb, exclude.pattern="Gohir.1Z")
TSS.A = TSS.all[grep("A",seqnames(TSS.all))] #36116
TSS.D = TSS.all[grep("D",seqnames(TSS.all))] #38784
con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")
con$D5<-gsub("[.]1$","",con$V5)
con$A2<-paste0(con$V2, ".gene")

# OG 22889 groups
TSS.Aog = TSS.all[con$V3]
TSS.Dog = TSS.all[con$V4]
# Not in OG
TSS.Aother = TSS.A[!names(TSS.A) %in% con$V3] #13227
TSS.Dother = TSS.D[!names(TSS.D) %in% con$V4] #15895

# bw from aggregation, mapped to AD1
bws = c("Qregulation_diffACRs/AD1ref/mapping/bwRPGC/ScD_q20.full.bw", "Qregulation_diffACRs/AD1ref/mapping/bwRPGC/McD_q20.full.bw", "Qregulation_diffACRs/AD1ref/mapping/bwRPGC/FcD_q20.full.bw")
bws = c("Qregulation4/ScD_q20.full.bw", "Qregulation4/McD_q20.full.bw", "Qregulation4/FcD_q20.full.bw")
names(bws)=c("Diploids","AD1","F1")
color=c("purple","blue","brown")


pdf("Qregulation4/plotDNS_RPGC.AD1ref.Wr.pdf")
plotBWsWindows(bws, list(A=TSS.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="All ", ylim=c(-0.1,1) )
plotBWsWindows(bws, list(D=TSS.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="All ", ylim=c(-0.1,1) )

plotBWsWindows(bws, list(A=TSS.Aog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.Dog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG ", ylim=c(-0.1,1))

## DE based on conserved
up.A= as.character(con[con$Wr=="Wr>0","V3"])
up.D= as.character(con[con$Wr=="Wr>0","V4"])
down.A= as.character(con[con$Wr=="Wr<0","V3"])
down.D= as.character(con[con$Wr=="Wr<0","V4"])
uc.A= as.character(con[con$Wr=="Wr=0","V3"])
uc.D= as.character(con[con$Wr=="Wr=0","V4"])
TSS.up.A =  TSS.A[up.A] #2122
TSS.up.D =  TSS.D[up.D]
TSS.down.A =  TSS.A[down.A] #2400
TSS.down.D =  TSS.D[down.D]
TSS.uc.A =  TSS.A[uc.A]
TSS.uc.D =  TSS.D[uc.D]

plotBWsWindows(bws, list(A=TSS.Aother), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Other ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.Dother), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Other ", ylim=c(-0.1,1))

plotBWsWindows(bws, list(A=TSS.up.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr>0 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.up.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr>0 ", ylim=c(-0.1,1))

plotBWsWindows(bws, list(A=TSS.down.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr<0 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.down.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr<0 ", ylim=c(-0.1,1))

plotBWsWindows(bws, list(A=TSS.uc.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr=0 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.uc.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr=0 ", ylim=c(-0.1,1))

dev.off()

###############
## dominance ## 2023
###############
con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")

txdb = loadDb("refGenomes/txdb.AD1utx.sqlite")
TSS.all = getTSS(txdb, exclude.pattern="Gohir.1Z")
TSS.A = TSS.all[grep("A",seqnames(TSS.all))] #36116
TSS.D = TSS.all[grep("D",seqnames(TSS.all))] #38784
con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")
con$D5<-gsub("[.]1$","",con$V5)
con$A2<-paste0(con$V2, ".gene")

# OG 22889 groups
TSS.Aog = TSS.all[con$V3]
TSS.Dog = TSS.all[con$V4]
# Not in OG
TSS.Aother = TSS.A[!names(TSS.A) %in% con$V3] #13227
TSS.Dother = TSS.D[!names(TSS.D) %in% con$V4] #15895

# bw from aggregation, mapped to AD1
# bws = c("Qregulation4/ScD_q20.full.bw", "Qregulation4/McD_q20.full.bw", "Qregulation4/FcD_q20.full.bw")
# names(bws)=c("Diploids","AD1","F1")
# color=c("purple","blue","brown")
bws = c("Qregulation4/ScD_q20.full.bw", "Qregulation4/McD_q20.full.bw")
names(bws)=c("Diploids","AD1")
color=c("purple","blue")
pdf("Qregulation4/plotDNS_RPGC.AD1ref.AD1dominance.pdf")

## DE based on conserved
TSS.ad1.Adominance.1.A =  TSS.A[as.character(con[con$AD1.dominance=="A-dominant"&con$A=="A>0","V3"])] #506
TSS.ad1.Adominance.1.D =  TSS.D[as.character(con[con$AD1.dominance=="A-dominant"&con$A=="A>0","V4"])]
TSS.ad1.Adominance.2.A =  TSS.A[as.character(con[con$AD1.dominance=="A-dominant"&con$A=="A<0","V3"])] #506
TSS.ad1.Adominance.2.D =  TSS.D[as.character(con[con$AD1.dominance=="A-dominant"&con$A=="A<0","V4"])]
TSS.ad1.Ddominance.1.A =  TSS.A[as.character(con[con$AD1.dominance=="D-dominant"&con$A=="A>0","V3"])] #1006
TSS.ad1.Ddominance.1.D =  TSS.D[as.character(con[con$AD1.dominance=="D-dominant"&con$A=="A>0","V4"])]
TSS.ad1.Ddominance.2.A =  TSS.A[as.character(con[con$AD1.dominance=="D-dominant"&con$A=="A<0","V3"])] #1006
TSS.ad1.Ddominance.2.D =  TSS.D[as.character(con[con$AD1.dominance=="D-dominant"&con$A=="A<0","V4"])]
TSS.ad1.Additivity.A =  TSS.A[as.character(con[con$AD1.dominance=="Additivity","V3"])] #11900
TSS.ad1.Additivity.D =  TSS.D[as.character(con[con$AD1.dominance=="Additivity","V4"])]
TSS.ad1.Tup.1.A =  TSS.A[as.character(con[con$AD1.dominance=="Transgressive Up"&con$A=="A>0","V3"])] #11900
TSS.ad1.Tup.1.D =  TSS.D[as.character(con[con$AD1.dominance=="Transgressive Up"&con$A=="A>0","V4"])]
TSS.ad1.Tup.2.A =  TSS.A[as.character(con[con$AD1.dominance=="Transgressive Up"&con$A=="A<0","V3"])] #11900
TSS.ad1.Tup.2.D =  TSS.D[as.character(con[con$AD1.dominance=="Transgressive Up"&con$A=="A<0","V4"])]
TSS.ad1.Tdown.1.A =  TSS.A[as.character(con[con$AD1.dominance=="Transgressive Down"&con$A=="A>0","V3"])] #11900
TSS.ad1.Tdown.1.D =  TSS.D[as.character(con[con$AD1.dominance=="Transgressive Down"&con$A=="A>0","V4"])]
TSS.ad1.Tdown.2.A =  TSS.A[as.character(con[con$AD1.dominance=="Transgressive Down"&con$A=="A<0","V3"])] #11900
TSS.ad1.Tdown.2.D =  TSS.D[as.character(con[con$AD1.dominance=="Transgressive Down"&con$A=="A<0","V4"])]


plotBWsWindows(bws, list(A=TSS.ad1.Adominance.1.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:A-dominant A2>D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Adominance.1.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:A-dominant A2>D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(A=TSS.ad1.Adominance.2.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:A-dominant A2<D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Adominance.2.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:A-dominant A2<d5 ", ylim=c(-0.1,1))

plotBWsWindows(bws, list(A=TSS.ad1.Ddominance.1.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:D-dominant A2>D5", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Ddominance.1.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:D-dominant A2>D5", ylim=c(-0.1,1))
plotBWsWindows(bws, list(A=TSS.ad1.Ddominance.2.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:D-dominant A2<D5", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Ddominance.2.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:D-dominant A2<D5", ylim=c(-0.1,1))

plotBWsWindows(bws, list(A=TSS.Aog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.Dog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(A=TSS.ad1.Additivity.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:Additivity ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Additivity.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:Additivity ", ylim=c(-0.1,1))

plotBWsWindows(bws, list(A=TSS.ad1.Tup.1.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:Transgressive Up A2>D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Tup.1.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:Transgressive Up A2>D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(A=TSS.ad1.Tup.2.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:Transgressive Up A2<D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Tup.2.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:Transgressive Up A2<D5 ", ylim=c(-0.1,1))

plotBWsWindows(bws, list(A=TSS.ad1.Tdown.1.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:Transgressive Down A2>D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Tdown.1.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:Transgressive Down A2>D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(A=TSS.ad1.Tdown.2.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:Transgressive Down A2<D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Tdown.2.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="AD1:Transgressive Down A2<D5 ", ylim=c(-0.1,1))

dev.off()

bws = c("Qregulation4/ScD_q20.full.bw", "Qregulation4/FcD_q20.full.bw")
names(bws)=c("Diploids","F1")
color=c("purple","brown")
pdf("Qregulation4/plotDNS_RPGC.AD1ref.F1dominance.pdf")

## DE based on conserved
TSS.ad1.Adominance.1.A =  TSS.A[as.character(con[con$F1.dominance=="A-dominant"&con$A=="A>0","V3"])] #506
TSS.ad1.Adominance.1.D =  TSS.D[as.character(con[con$F1.dominance=="A-dominant"&con$A=="A>0","V4"])]
TSS.ad1.Adominance.2.A =  TSS.A[as.character(con[con$F1.dominance=="A-dominant"&con$A=="A<0","V3"])] #506
TSS.ad1.Adominance.2.D =  TSS.D[as.character(con[con$F1.dominance=="A-dominant"&con$A=="A<0","V4"])]
TSS.ad1.Ddominance.1.A =  TSS.A[as.character(con[con$F1.dominance=="D-dominant"&con$A=="A>0","V3"])] #1006
TSS.ad1.Ddominance.1.D =  TSS.D[as.character(con[con$F1.dominance=="D-dominant"&con$A=="A>0","V4"])]
TSS.ad1.Ddominance.2.A =  TSS.A[as.character(con[con$F1.dominance=="D-dominant"&con$A=="A<0","V3"])] #1006
TSS.ad1.Ddominance.2.D =  TSS.D[as.character(con[con$F1.dominance=="D-dominant"&con$A=="A<0","V4"])]
TSS.ad1.Additivity.A =  TSS.A[as.character(con[con$F1.dominance=="Additivity","V3"])] #11900
TSS.ad1.Additivity.D =  TSS.D[as.character(con[con$F1.dominance=="Additivity","V4"])]
TSS.ad1.Tup.1.A =  TSS.A[as.character(con[con$F1.dominance=="Transgressive Up"&con$A=="A>0","V3"])] #11900
TSS.ad1.Tup.1.D =  TSS.D[as.character(con[con$F1.dominance=="Transgressive Up"&con$A=="A>0","V4"])]
TSS.ad1.Tup.2.A =  TSS.A[as.character(con[con$F1.dominance=="Transgressive Up"&con$A=="A<0","V3"])] #11900
TSS.ad1.Tup.2.D =  TSS.D[as.character(con[con$F1.dominance=="Transgressive Up"&con$A=="A<0","V4"])]
TSS.ad1.Tdown.1.A =  TSS.A[as.character(con[con$F1.dominance=="Transgressive Down"&con$A=="A>0","V3"])] #11900
TSS.ad1.Tdown.1.D =  TSS.D[as.character(con[con$F1.dominance=="Transgressive Down"&con$A=="A>0","V4"])]
TSS.ad1.Tdown.2.A =  TSS.A[as.character(con[con$F1.dominance=="Transgressive Down"&con$A=="A<0","V3"])] #11900
TSS.ad1.Tdown.2.D =  TSS.D[as.character(con[con$F1.dominance=="Transgressive Down"&con$A=="A<0","V4"])]


plotBWsWindows(bws, list(A=TSS.ad1.Adominance.1.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:A-dominant A2>D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Adominance.1.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:A-dominant A2>D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(A=TSS.ad1.Adominance.2.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:A-dominant A2<D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Adominance.2.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:A-dominant A2<d5 ", ylim=c(-0.1,1))

plotBWsWindows(bws, list(A=TSS.ad1.Ddominance.1.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:D-dominant A2>D5", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Ddominance.1.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:D-dominant A2>D5", ylim=c(-0.1,1))
plotBWsWindows(bws, list(A=TSS.ad1.Ddominance.2.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:D-dominant A2<D5", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Ddominance.2.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:D-dominant A2<D5", ylim=c(-0.1,1))

plotBWsWindows(bws, list(A=TSS.Aog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.Dog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(A=TSS.ad1.Additivity.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:Additivity ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Additivity.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:Additivity ", ylim=c(-0.1,1))

plotBWsWindows(bws, list(A=TSS.ad1.Tup.1.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:Transgressive Up A2>D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Tup.1.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:Transgressive Up A2>D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(A=TSS.ad1.Tup.2.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:Transgressive Up A2<D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Tup.2.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:Transgressive Up A2<D5 ", ylim=c(-0.1,1))

plotBWsWindows(bws, list(A=TSS.ad1.Tdown.1.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:Transgressive Down A2>D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Tdown.1.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:Transgressive Down A2>D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(A=TSS.ad1.Tdown.2.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:Transgressive Down A2<D5 ", ylim=c(-0.1,1))
plotBWsWindows(bws, list(D=TSS.ad1.Tdown.2.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="F1:Transgressive Down A2<D5 ", ylim=c(-0.1,1))

dev.off()

---2023

pdf("Qregulation4/plotDNS_RPGC.AD1ref.ALL.pdf") #FigS7
#all
plotBWsWindows(bws, list(A=TSS.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="All " )
plotBWsWindows(bws, list(D=TSS.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="All ")
#og
plotBWsWindows(bws, list(A=TSS.Aog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG ")
plotBWsWindows(bws, list(D=TSS.Dog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG ")
#other
plotBWsWindows(bws, list(A=TSS.Aother), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Other ")
plotBWsWindows(bws, list(D=TSS.Dother), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Other ")
## A
Abias.A= as.character(con[con$A=="A>0","V3"])
Abias.D= as.character(con[con$A=="A>0","V4"])
Dbias.A= as.character(con[con$A=="A<0","V3"])
Dbias.D= as.character(con[con$A=="A<0","V4"])
unbias.A= as.character(con[con$A=="A=0","V3"])
unbias.D= as.character(con[con$A=="A=0","V4"])
TSS.Abias.A =  TSS.A[Abias.A] #1507
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A] #1416
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
#
plotBWsWindows(bws, list(A=TSS.Abias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A>0 " )
plotBWsWindows(bws, list(D=TSS.Abias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A>0 " )
#
plotBWsWindows(bws, list(A=TSS.Dbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A<0 " )
plotBWsWindows(bws, list(D=TSS.Dbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A<0 " )
#
plotBWsWindows(bws, list(A=TSS.unbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A=0 " )
plotBWsWindows(bws, list(D=TSS.unbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A=0 " )
## B
Abias.A= as.character(con[con$B=="B>0","V3"])
Abias.D= as.character(con[con$B=="B>0","V4"])
Dbias.A= as.character(con[con$B=="B<0","V3"])
Dbias.D= as.character(con[con$B=="B<0","V4"])
unbias.A= as.character(con[con$B=="B=0","V3"])
unbias.D= as.character(con[con$B=="B=0","V4"])
TSS.Abias.A =  TSS.A[Abias.A] #807
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A] #676
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
#
plotBWsWindows(bws, list(A=TSS.Abias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B>0 " )
plotBWsWindows(bws, list(D=TSS.Abias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B>0 " )
#
plotBWsWindows(bws, list(A=TSS.Dbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B<0 " )
plotBWsWindows(bws, list(D=TSS.Dbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B<0 " )
#
plotBWsWindows(bws, list(A=TSS.unbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B=0 " )
plotBWsWindows(bws, list(D=TSS.unbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B=0 " )
## Bp
Abias.A= as.character(con[con$Bp=="Bp>0","V3"])
Abias.D= as.character(con[con$Bp=="Bp>0","V4"])
Dbias.A= as.character(con[con$Bp=="Bp<0","V3"])
Dbias.D= as.character(con[con$Bp=="Bp<0","V4"])
unbias.A= as.character(con[con$Bp=="Bp=0","V3"])
unbias.D= as.character(con[con$Bp=="Bp=0","V4"])
TSS.Abias.A =  TSS.A[Abias.A] #1682
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A] #1550
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
#
plotBWsWindows(bws, list(A=TSS.Abias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Bp>0 " )
plotBWsWindows(bws, list(D=TSS.Abias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Bp>0 " )
#
plotBWsWindows(bws, list(A=TSS.Dbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Bp<0 " )
plotBWsWindows(bws, list(D=TSS.Dbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Bp<0 " )
#
plotBWsWindows(bws, list(A=TSS.unbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Bp=0 " )
plotBWsWindows(bws, list(D=TSS.unbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Bp=0 " )
## Hr
up.A= as.character(con[con$Hr=="Hr>0","V3"])
up.D= as.character(con[con$Hr=="Hr>0","V4"])
down.A= as.character(con[con$Hr=="Hr<0","V3"])
down.D= as.character(con[con$Hr=="Hr<0","V4"])
uc.A= as.character(con[con$Hr=="Hr=0","V3"])
uc.D= as.character(con[con$Hr=="Hr=0","V4"])
TSS.up.A =  TSS.A[up.A] #80
TSS.up.D =  TSS.D[up.D]
TSS.down.A =  TSS.A[down.A] #62
TSS.down.D =  TSS.D[down.D]
TSS.uc.A =  TSS.A[uc.A]
TSS.uc.D =  TSS.D[uc.D]
plotBWsWindows(bws, list(A=TSS.up.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr>0 ")
plotBWsWindows(bws, list(D=TSS.up.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr>0 ")
#
plotBWsWindows(bws, list(A=TSS.down.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr<0 ")
plotBWsWindows(bws, list(D=TSS.down.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr<0 ")
#
plotBWsWindows(bws, list(A=TSS.uc.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr=0 ")
plotBWsWindows(bws, list(D=TSS.uc.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr=0 ")
## Wr
up.A= as.character(con[con$Wr=="Wr>0","V3"])
up.D= as.character(con[con$Wr=="Wr>0","V4"])
down.A= as.character(con[con$Wr=="Wr<0","V3"])
down.D= as.character(con[con$Wr=="Wr<0","V4"])
uc.A= as.character(con[con$Wr=="Wr=0","V3"])
uc.D= as.character(con[con$Wr=="Wr=0","V4"])
TSS.up.A =  TSS.A[up.A] #283
TSS.up.D =  TSS.D[up.D]
TSS.down.A =  TSS.A[down.A] #396
TSS.down.D =  TSS.D[down.D]
TSS.uc.A =  TSS.A[uc.A]
TSS.uc.D =  TSS.D[uc.D]
plotBWsWindows(bws, list(A=TSS.up.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr>0 ")
plotBWsWindows(bws, list(D=TSS.up.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr>0 ")
#
plotBWsWindows(bws, list(A=TSS.down.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr<0 ")
plotBWsWindows(bws, list(D=TSS.down.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr<0 ")
#
plotBWsWindows(bws, list(A=TSS.uc.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr=0 ")
plotBWsWindows(bws, list(D=TSS.uc.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr=0 ")
#
dev.off()





###############
## A2+D5 ref ##
###############
con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")
con$D5<-gsub("[.]1$","",con$V5)
con$A2<-paste0(con$V2, ".gene")

txdb = loadDb("refGenomes/txdb.A2WHU.sqlite")
TSS = getTSS(txdb, exclude.pattern="Gorai.N")
TSS.A = TSS[grepl("Gar",names(TSS))&!grepl("Contig",seqnames(TSS))] #41739
txdb = loadDb("refGenomes/txdb.D5.sqlite")
TSS = getTSS(txdb, exclude.pattern="Gorai.N")
TSS.D = TSS[grep("Gorai",names(TSS))] #37218
## OG 22889 groups
TSS.Aog = TSS.A[con$A2]
TSS.Dog = TSS.D[con$D5]
# Not in OG
TSS.Aother = TSS.A[!names(TSS.A) %in% con$A2] #18850
TSS.Dother = TSS.D[!names(TSS.D) %in% con$D5] #14329

bws=c("Qregulation_diffACRs/Fref/bwRPGC/A6Dn_q20.full.bw", "Qregulation_diffACRs/Fref/bwRPGC/DcD_q20.full.bw", "Qregulation_diffACRs/Fref/bwRPGC/FcDa_q20.full.bw", "Qregulation_diffACRs/Fref/bwRPGC/FcDd_q20.full.bw", "Qregulation_diffACRs/Fref/bwRPGC/McDa_q20.full.bw", "Qregulation_diffACRs/Fref/bwRPGC/McDd_q20.full.bw")
names(bws)=c("A2","D5","F1.At","F1.Dt","AD1.At","AD1.Dt")
color=c("purple","brown","blue")

pdf("Qregulation4/plotDNS_RPGC.A2D5ref.og.pdf")
#
plotBWsWindows(bws[c("A2","F1.At")], list(A=TSS.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="All ", ylim=c(-0.1,0.5))
plotBWsWindows(bws[c("D5","F1.Dt")], list(D=TSS.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="All ", ylim=c(-0.1,0.5) )
#
plotBWsWindows(bws[c("A2","F1.At")], list(A=TSS.Aog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG ", ylim=c(-0.1,0.5))
plotBWsWindows(bws[c("D5","F1.Dt")], list(D=TSS.Dog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG ", ylim=c(-0.1,0.5))
#
plotBWsWindows(bws[c("A2","F1.At")], list(A=TSS.Aother), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Other ", ylim=c(-0.1,0.5))
plotBWsWindows(bws[c("D5","F1.Dt")], list(D=TSS.Dother), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Other ", ylim=c(-0.1,0.5))
#
dev.off()

## A
Abias.A= as.character(con[con$A=="A>0","A2"])
Abias.D= as.character(con[con$A=="A>0","D5"])
Dbias.A= as.character(con[con$A=="A<0","A2"])
Dbias.D= as.character(con[con$A=="A<0","D5"])
unbias.A= as.character(con[con$A=="A=0","A2"])
unbias.D= as.character(con[con$A=="A=0","D5"])
TSS.Abias.A =  TSS.A[Abias.A] #1507
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A] #1416
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
pdf("Qregulation4/plotDNS_RPGC.A2D5ref.A.pdf")
#
plotBWsWindows(bws[c("A2","F1.At")], list(A=TSS.Abias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A>0 ", ylim=c(-0.1,0.5))
plotBWsWindows(bws[c("D5","F1.Dt")], list(D=TSS.Abias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A>0 ", ylim=c(-0.1,0.5) )
#
plotBWsWindows(bws[c("A2","F1.At")], list(A=TSS.Dbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A<0 ", ylim=c(-0.1,0.5))
plotBWsWindows(bws[c("D5","F1.Dt")], list(D=TSS.Dbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A<0 ", ylim=c(-0.1,0.5) )
#
plotBWsWindows(bws[c("A2","F1.At")], list(A=TSS.unbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A=0 ", ylim=c(-0.1,0.5))
plotBWsWindows(bws[c("D5","F1.Dt")], list(D=TSS.unbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A=0 ", ylim=c(-0.1,0.5) )
dev.off()

## B
Abias.A= as.character(con[con$B=="B>0","A2"])
Abias.D= as.character(con[con$B=="B>0","D5"])
Dbias.A= as.character(con[con$B=="B<0","A2"])
Dbias.D= as.character(con[con$B=="B<0","D5"])
unbias.A= as.character(con[con$B=="B=0","A2"])
unbias.D= as.character(con[con$B=="B=0","D5"])
TSS.Abias.A =  TSS.A[Abias.A] #1507
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A] #1416
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
pdf("Qregulation4/plotDNS_RPGC.A2D5ref.B.pdf")
#
plotBWsWindows(bws[c("A2","F1.At")], list(A=TSS.Abias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B>0 ", ylim=c(-0.1,0.5))
plotBWsWindows(bws[c("D5","F1.Dt")], list(D=TSS.Abias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B>0 ", ylim=c(-0.1,0.5) )
#
plotBWsWindows(bws[c("A2","F1.At")], list(A=TSS.Dbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B<0 ", ylim=c(-0.1,0.5))
plotBWsWindows(bws[c("D5","F1.Dt")], list(D=TSS.Dbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B<0 ", ylim=c(-0.1,0.5) )
#
plotBWsWindows(bws[c("A2","F1.At")], list(A=TSS.unbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B=0 ", ylim=c(-0.1,0.5))
plotBWsWindows(bws[c("D5","F1.Dt")], list(D=TSS.unbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B=0 ", ylim=c(-0.1,0.5) )
dev.off()

## Hr
up.A= as.character(con[con$Hr=="Hr>0","A2"])
up.D= as.character(con[con$Hr=="Hr>0","D5"])
down.A= as.character(con[con$Hr=="Hr<0","A2"])
down.D= as.character(con[con$Hr=="Hr<0","D5"])
uc.A= as.character(con[con$Hr=="Hr=0","A2"])
uc.D= as.character(con[con$Hr=="Hr=0","D5"])
TSS.up.A =  TSS.A[up.A] #2122
TSS.up.D =  TSS.D[up.D]
TSS.down.A =  TSS.A[down.A] #2400
TSS.down.D =  TSS.D[down.D]
TSS.uc.A =  TSS.A[uc.A]
TSS.uc.D =  TSS.D[uc.D]
pdf("Qregulation4/plotDNS_RPGC.A2D5ref.Hr.pdf")
plotBWsWindows(bws[c("A2","F1.At")], list(A=TSS.up.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr>0 ", ylim=c(-0.1,0.5))
plotBWsWindows(bws[c("D5","F1.Dt")], list(D=TSS.up.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr>0 ", ylim=c(-0.1,0.5))
#
plotBWsWindows(bws[c("A2","F1.At")], list(A=TSS.down.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr<0 ", ylim=c(-0.1,0.5))
plotBWsWindows(bws[c("D5","F1.Dt")], list(D=TSS.down.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr<0 ", ylim=c(-0.1,0.5))
#
plotBWsWindows(bws[c("A2","F1.At")], list(A=TSS.uc.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr=0 ", ylim=c(-0.1,0.5))
plotBWsWindows(bws[c("D5","F1.Dt")], list(D=TSS.uc.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr=0 ", ylim=c(-0.1,0.5))
dev.off()

---
pdf("Qregulation4/plotDNS_RPGC.A2D5ref.ALL.pdf") #FigS6
#
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="All ")
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="All ")
#
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.Aog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG " )
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.Dog), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="OG " )
#
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.Aother), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Other " )
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.Dother), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Other " )

## A
Abias.A= as.character(con[con$A=="A>0","A2"])
Abias.D= as.character(con[con$A=="A>0","D5"])
Dbias.A= as.character(con[con$A=="A<0","A2"])
Dbias.D= as.character(con[con$A=="A<0","D5"])
unbias.A= as.character(con[con$A=="A=0","A2"])
unbias.D= as.character(con[con$A=="A=0","D5"])
TSS.Abias.A =  TSS.A[Abias.A] #1507
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A] #1416
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
#
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.Abias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A>0 " )
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.Abias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A>0 " )
#
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.Dbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A<0 " )
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.Dbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A<0 " )
#
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.unbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A=0 " )
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.unbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="A=0 " )


## B
Abias.A= as.character(con[con$B=="B>0","A2"])
Abias.D= as.character(con[con$B=="B>0","D5"])
Dbias.A= as.character(con[con$B=="B<0","A2"])
Dbias.D= as.character(con[con$B=="B<0","D5"])
unbias.A= as.character(con[con$B=="B=0","A2"])
unbias.D= as.character(con[con$B=="B=0","D5"])
TSS.Abias.A =  TSS.A[Abias.A] #1507
TSS.Abias.D =  TSS.D[Abias.D]
TSS.Dbias.A =  TSS.A[Dbias.A] #1416
TSS.Dbias.D =  TSS.D[Dbias.D]
TSS.unbias.A =  TSS.A[unbias.A]
TSS.unbias.D =  TSS.D[unbias.D]
#
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.Abias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B>0 " )
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.Abias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B>0 " )
#
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.Dbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B<0 " )
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.Dbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B<0 " )
#
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.unbias.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B=0 " )
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.unbias.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="B=0 " )

## Hr
up.A= as.character(con[con$Hr=="Hr>0","A2"])
up.D= as.character(con[con$Hr=="Hr>0","D5"])
down.A= as.character(con[con$Hr=="Hr<0","A2"])
down.D= as.character(con[con$Hr=="Hr<0","D5"])
uc.A= as.character(con[con$Hr=="Hr=0","A2"])
uc.D= as.character(con[con$Hr=="Hr=0","D5"])
TSS.up.A =  TSS.A[up.A] #2122
TSS.up.D =  TSS.D[up.D]
TSS.down.A =  TSS.A[down.A] #2400
TSS.down.D =  TSS.D[down.D]
TSS.uc.A =  TSS.A[uc.A]
TSS.uc.D =  TSS.D[uc.D]
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.up.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr>0 ")
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.up.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr>0 ")
#
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.down.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr<0 ")
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.down.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr<0 ")
#
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.uc.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr=0 ")
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.uc.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Hr=0 ")

## Wr
up.A= as.character(con[con$Wr=="Wr>0","A2"])
up.D= as.character(con[con$Wr=="Wr>0","D5"])
down.A= as.character(con[con$Wr=="Wr<0","A2"])
down.D= as.character(con[con$Wr=="Wr<0","D5"])
uc.A= as.character(con[con$Wr=="Wr=0","A2"])
uc.D= as.character(con[con$Wr=="Wr=0","D5"])
TSS.up.A =  TSS.A[up.A] #2122
TSS.up.D =  TSS.D[up.D]
TSS.down.A =  TSS.A[down.A] #2400
TSS.down.D =  TSS.D[down.D]
TSS.uc.A =  TSS.A[uc.A]
TSS.uc.D =  TSS.D[uc.D]
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.up.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr>0 ")
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.up.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr>0 ")
#
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.down.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr<0 ")
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.down.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr<0 ")
#
plotBWsWindows(bws[c("A2","F1.At","AD1.At")], list(A=TSS.uc.A), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr=0 ")
plotBWsWindows(bws[c("D5","F1.Dt","AD1.Dt")], list(D=TSS.uc.D), compare="bw", line.col=color, dispersion.col=makeTransparent(color), main="Wr=0 ")

dev.off()


###############################
## identified ACR comparison ##
###############################
# load p and g ACRs
#setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions")
#load('ACRs_compare/ACRpg.rdata')
load("Qregulation3/ACRpg.rdata")
# load duplicated gene patterns
con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")

# Func
getGrInfo = function(gr){
    anno =gr
    catLevel = c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Downstream (<1kb)", "Downstream (1-2kb)", "Downstream (2-3kb)", "Exon", "Intron", "Distal Intergenic" , "5' UTR","3' UTR")
    category = factor(anno$annotation2,levels = catLevel)
    df =data.frame(Feature=catLevel, RegionSize = tapply(width(anno),category,sum), Number = tapply(rep(1,length(category)),category,sum) , geneN=tapply(anno$geneId,category,function(x)length(unique(x))))
    rownames(df)=NULL
    df[is.na(df)]=0
    df$NumberPerGene = df$Number/df$geneN
    df$RegionSizePerGene = df$RegionSize/df$geneN
    return(df)
}
dflist2df = function(df.list,label="Sample")
{
    samples = names(df.list)
    df = df.list[[1]]
    df[,label]=samples[1]
    for(i in 2:length(samples)){
        x = df.list[[i]]
        x[,label] = samples[i]
        df = rbind(df, x)
    }
    # levels(df$Feature)=c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Exon", "Intron", "Downstream (<=300)","Distal Intergenic" )
    # df$Sample=factor(df$Sample, levels=c("A2","D5","F1.At","F1.Dt","AD1.At","AD1.Dt"))
    return(df)
}
library(ggplot2)
library(gridExtra)
library(IRanges)

# get promoter ACR size per AD1 gene
# At
gr<-anL$AD1a@anno
gr=gr[gr$annotation2=="Promoter (<=1kb)"]
da<-data.frame(id=gr$geneId,feature=gr$annotation2, size=width(gr))
sa<-aggregate(da$size, list(da$id),sum)
names(sa)=c("id","AD1a")
sa$V1<- con$V1[match(sa$id,con$V3)]
# Dt
gr<-anL$AD1d@anno
gr=gr[gr$annotation2=="Promoter (<=1kb)"]
dd<-data.frame(id=gr$geneId,feature=gr$annotation2, size=width(gr))
sd<-aggregate(dd$size, list(dd$id),sum)
names(sd)=c("id","AD1d")
sd$V1<- con$V1[match(sd$id,con$V4)]
# assemble by OG
cc<-con
ccc<-merge(merge(cc,sa[,c("V1","AD1a")],allx=T),sd[,c("V1","AD1d")],all.x=T)
ccc$AD1a[is.na(ccc$AD1a)]=0
ccc$AD1d[is.na(ccc$AD1d)]=0
# examine ACR sizes, At always higher than Dt regardless of bias
boxplot(ccc$AD1a[ccc$Bp=="Bp<0"], ccc$AD1d[ccc$Bp=="Bp<0"], main="Bp<0")
boxplot(ccc$AD1a[ccc$Bp=="Bp>0"], ccc$AD1d[ccc$Bp=="Bp>0"], main="Bp>0")
boxplot(ccc$AD1a[ccc$Bp=="Bp=0"], ccc$AD1d[ccc$Bp=="Bp=0"], main="Bp=0")
# No clear pattern from scatterplot either
plot(ccc$AD1a,ccc$AD1d,col=factor(ccc$Bp))

# get average size (totalSize/geneNumber)

# AD1
gr<-anL$AD1a@anno
gr$test <-  con$Bp[match(gr$geneId,con$V3)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$AD1d@anno
gr$test <-  con$Bp[match(gr$geneId,con$V4)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pAD1bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1") + theme_bw()
pAD1bar
pAD1box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1") + theme_bw()
pAD1box


# F1
gr<-anL$F1a@anno
gr$test <-  con$B[match(gr$transcriptId,con$V2)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$F1d@anno
gr$test <-  con$B[match(paste0(gr$geneId,".1"),con$V5)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pF1bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1") + theme_bw()
pF1box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1") + theme_bw()
pF1bar
pF1box

# A2 D5
gr<-anL$A2@anno
gr$test <-  con$A[match(gr$transcriptId,con$V2)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$D5@anno
gr$test <-  con$A[match(paste0(gr$geneId,".1"),con$V5)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pA2D5bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("A2 vs D5") + theme_bw()
pA2D5bar
pA2D5box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("A2 vs D5") + theme_bw()
pA2D5box

pdf("Qregulation4/acr2DGE.pdf")
pA2D5bar
pF1bar
pAD1bar
grid.arrange(pA2D5bar+theme(legend.position="none"),pF1bar+theme(legend.position="none"),pAD1bar+theme(legend.position="none"),nrow=3,ncol=1)
pA2D5box
pF1box
pAD1box
grid.arrange(pA2D5box+theme(legend.position="none"),pF1box+theme(legend.position="none"),pAD1box+theme(legend.position="none"),nrow=3,ncol=1)
dev.off()


##########################
## DA result comparison ##
##########################
library(ChIPseeker)
library(ggplot2)
library(gridExtra)
# load duplicated gene patterns
con<- read.table("Qregulation4/conserved.regPattern.txt",header=T, sep="\t")
# Func
getGrInfo = function(gr){
    anno =gr
    catLevel = c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Downstream (<1kb)", "Downstream (1-2kb)", "Downstream (2-3kb)", "Exon", "Intron", "Distal Intergenic" , "5' UTR","3' UTR")
    category = factor(anno$annotation2,levels = catLevel)
    df =data.frame(Feature=catLevel, RegionSize = tapply(width(anno),category,sum), Number = tapply(rep(1,length(category)),category,sum) , geneN=tapply(anno$geneId,category,function(x)length(unique(x))))
    rownames(df)=NULL
    df[is.na(df)]=0
    df$NumberPerGene = df$Number/df$geneN
    df$RegionSizePerGene = df$RegionSize/df$geneN
    return(df)
}
dflist2df = function(df.list,label="Sample")
{
    samples = names(df.list)
    df = df.list[[1]]
    df[,label]=samples[1]
    if(length(samples)>1){
        for(i in 2:length(samples)){
            x = df.list[[i]]
            x[,label] = samples[i]
            df = rbind(df, x)
        }
    }
    # levels(df$Feature)=c( "Promoter (<=1kb)", "Promoter (1-2kb)", "Promoter (2-3kb)", "Exon", "Intron", "Downstream (<=300)","Distal Intergenic" )
    # df$Sample=factor(df$Sample, levels=c("A2","D5","F1.At","F1.Dt","AD1.At","AD1.Dt"))
    return(df)
}

# load combined DA
#setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions")
#load('~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/DAanalysis/AD1ref/merged.annotation.rdata')
load("Qregulation_diffACRs/AD1ref/merged.annotation.rdata")
names(anL)
# "resDiA" "resFA"  "resMA"  "resHrA" "resWrA" "resPrA" "resDiD" "resFD"  "resMD"  "resHrD" "resWrD" "resPrD"
#
# AD1
gr<-anL$resMA@anno
gr$test <-  con$Bp[match(gr$geneId,con$V3)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$resMD@anno
gr$test <-  con$Bp[match(gr$geneId,con$V4)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pAD1bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1 ref: AD1") + theme_bw()
pAD1bar
pAD1box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1 ref: AD1") + theme_bw()
pAD1box

# F1
gr<-anL$resFA@anno
gr$test <-  con$B[match(gr$geneId,con$V3)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$resFD@anno
gr$test <-  con$B[match(gr$geneId,con$V4)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pF1bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1 ref: F1") + theme_bw()
pF1box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1 ref: F1") + theme_bw()
pF1bar
pF1box

# A2 D5
gr<-anL$resDiA@anno
gr$test <-  con$A[match(gr$geneId,con$V3)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$resDiD@anno
gr$test <-  con$A[match(gr$geneId,con$V4)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pA2D5bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1 ref: A2 vs D5") + theme_bw()
pA2D5bar
pA2D5box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("AD1 ref: A2 vs D5") + theme_bw()
pA2D5box

pdf("Qregulation4/DAad1ref2DGE.pdf")
pA2D5bar
pF1bar
pAD1bar
grid.arrange(pA2D5bar+theme(legend.position="none"),pF1bar+theme(legend.position="none"),pAD1bar+theme(legend.position="none"),nrow=3,ncol=1)
pA2D5box
pF1box
pAD1box
grid.arrange(pA2D5box+theme(legend.position="none"),pF1box+theme(legend.position="none"),pAD1box+theme(legend.position="none"),nrow=3,ncol=1)
dev.off()

# Hr: DA up or down has no association with expression Hr
gr<-anL$resHrA@anno
gr$test <-  con$Hr[match(gr$geneId,con$V3)]
table(gr$test,gr$direction)
#       down up
# -       1  1
# Hr=0   10 10
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="A"
dfA<-df
gr<-anL$resHrD@anno
gr$test <-  con$Hr[match(gr$geneId,con$V4)]
table(gr$test) # NA
table(gr$direction)
# down   up
#   4    4
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
df$res="resHr"
ddf<-df
# Pr
gr<-anL$resPrA@anno
gr$test <-  con$Pr[match(gr$geneId,con$V3)]
table(gr$test,gr$direction)
#       down mixed up
# -       1     0 14
# Pr<0    1     0  6
# Pr=0   11     0 95
# Pr>0    0     0  3
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="A"
dfA<-df
gr<-anL$resPrD@anno
gr$test <-  con$Pr[match(gr$geneId,con$V4)]
table(gr$test,gr$direction)
#      down up
# -       0  2
# Pr=0    5 16
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
df$res="resPr"
ddf<-rbind(ddf,df)
# Wr
gr<-anL$resWrA@anno
gr$test <-  con$Wr[match(gr$geneId,con$V3)]
table(gr$test,gr$direction)
#      down up
# -       0 11
# Wr<0    0  5
# Wr=0    6 57
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="A"
dfA<-df
gr<-anL$resWrD@anno
gr$test <-  con$Wr[match(gr$geneId,con$V4)]
table(gr$test,gr$direction)
#       down up
# -       1  1
# Wr=0    6 18
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
df$res="resWr"
ddf<-rbind(ddf,df)

df<-ddf[ddf$Feature=="Promoter (<=1kb)"&ddf$direction!="mixed",]
p1<-ggplot(df, aes_string(x = "direction", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~res) + ylab("DA region per gene (bp)") + xlab("") + ggtitle("AD1 ref") + theme_bw()
p2<-ggplot(df, aes_string(x = "direction", fill = "genome", y = "RegionSize")) + geom_bar(position="dodge",stat="identity") + facet_grid(~res) + ylab("DA region (bp)") + xlab("") + ggtitle("AD1 ref") + theme_bw()
pdf("Qregulation4/DAad1ref2p1k.pdf")
grid.arrange(p1+theme(legend.position="none"),p2+theme(legend.position="none"),nrow=3,ncol=1)
p1
p2
dev.off()





# load combined DA
#setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions")
#load('~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/DAanalysis/F1ref/merged.annotation.rdata')
library(ChIPseeker)
load("Qregulation_diffACRs/Fref/merged.annotation.rdata")
names(anL)
# "resDiA" "resFA"  "resMA"  "resHrA" "resWrA" "resPrA" "resDiD" "resFD"  "resMD"  "resHrD" "resWrD" "resPrD"
#
# AD1
gr<-anL$resMA@anno
gr$test <-  con$Bp[match(gr$transcriptId,con$V2)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$resMD@anno
gr$test <-  con$Bp[match(gr$geneId,gsub("[.]1$","",con$V5))]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pAD1bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1 ref: AD1") + theme_bw()
pAD1bar
pAD1box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1 ref: AD1") + theme_bw()
pAD1box

# F1
gr<-anL$resFA@anno
gr$test <-  con$B[match(gr$transcriptId,con$V2)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$resFD@anno
gr$test <-  con$B[match(gr$geneId,gsub("[.]1$","",con$V5))]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pF1bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1 ref: F1") + theme_bw()
pF1box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1 ref: F1") + theme_bw()
pF1bar
pF1box

# A2 D5
gr<-anL$resDiA@anno
gr$test <-  con$A[match(gr$transcriptId,con$V2)]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="A"
dfA<-df
gr<-anL$resDiD@anno
gr$test <-  con$A[match(gr$geneId,gsub("[.]1$","",con$V5))]
grl<-split(gr, factor(seqnames(gr)))
df<-dflist2df(lapply(grl, function(gr)dflist2df(tapply(gr,gr$test,getGrInfo))),"chr")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
#
pA2D5bar<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1 ref: A2 vs D5") + theme_bw()
pA2D5bar
pA2D5box<-ggplot(df[grep("Promoter",df$Feature),], aes_string(x = "Sample", fill = "genome", y = "RegionSizePerGene")) + geom_boxplot() + facet_grid(~Feature) + ylab("Peak region per gene (bp)") + xlab("") + ggtitle("F1 ref: A2 vs D5") + theme_bw()
pA2D5box

pdf("Qregulation4/DAf1ref2DGE.pdf")
pA2D5bar
pF1bar
pAD1bar
grid.arrange(pA2D5bar+theme(legend.position="none"),pF1bar+theme(legend.position="none"),pAD1bar+theme(legend.position="none"),nrow=3,ncol=1)
pA2D5box
pF1box
pAD1box
grid.arrange(pA2D5box+theme(legend.position="none"),pF1box+theme(legend.position="none"),pAD1box+theme(legend.position="none"),nrow=3,ncol=1)
dev.off()


# Hr: DA up or down has no association with expression Hr
gr<-anL$resHrA@anno
gr$test <-  con$Hr[match(gr$transcriptId,con$V2)]
table(gr$test) # Hr=0 1
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="A"
dfA<-df
gr<-anL$resHrD@anno
gr$test <-  con$Hr[match(gr$geneId,gsub("[.]1$","",con$V5))]
table(gr$test) # Hr=0 7
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
df$res="resHr"
ddf<-df
# Pr
gr<-anL$resPrA@anno
gr$test <-  con$Pr[match(gr$transcriptId,con$V2)]
table(gr$test,gr$direction)
#      down mixed  up
# -       0     0  15
# Pr<0    0     0   5
# Pr=0    8     0 125
# Pr>0    0     0   3
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="A"
dfA<-df
gr<-anL$resPrD@anno
gr$test <-  con$Pr[match(gr$geneId,gsub("[.]1$","",con$V5))]
table(gr$test,gr$direction)
#      down up
# -       1  5
# Pr<0    1  1
# Pr=0   17 33
# Pr>0    1  0
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
df$res="resPr"
ddf<-rbind(ddf,df)
# Wr
gr<-anL$resWrA@anno
gr$test <-  con$Wr[match(gr$transcriptId,con$V2)]
table(gr$test,gr$direction)
#      down up
# -       1 14
# Wr<0    0  2
# Wr=0   11 96
# Wr>0    0  3
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="A"
dfA<-df
gr<-anL$resWrD@anno
gr$test <-  con$Wr[match(gr$geneId,gsub("[.]1$","",con$V5))]
table(gr$test,gr$direction)
#       down up
# -       0  9
# Wr<0    0  2
# Wr=0   10 28
# Wr>0    0  1
df<-dflist2df(tapply(gr,gr$direction,getGrInfo),"direction")
df$genome="D"
dfD<-df
df<-rbind(dfA,dfD)
df$res="resWr"
ddf<-rbind(ddf,df)

df<-ddf[ddf$Feature=="Promoter (<=1kb)"&ddf$direction!="mixed",]
p1<-ggplot(df, aes_string(x = "direction", fill = "genome", y = "RegionSizePerGene")) + geom_bar(position="dodge",stat="identity") + facet_grid(~res) + ylab("DA region per gene (bp)") + xlab("") + ggtitle("F1 ref") + theme_bw()
p2<-ggplot(df, aes_string(x = "direction", fill = "genome", y = "RegionSize")) + geom_bar(position="dodge",stat="identity") + facet_grid(~res) + ylab("DA region (bp)") + xlab("") + ggtitle("F1 ref") + theme_bw()
pdf("Qregulation4/DAf1ref2p1k.pdf")
grid.arrange(p1+theme(legend.position="none"),p2+theme(legend.position="none"),nrow=3,ncol=1)
p1
p2
dev.off()

