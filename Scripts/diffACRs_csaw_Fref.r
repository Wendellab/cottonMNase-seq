# diffACRs_csaw_Fref.r

# csaw workflow for ACR differential accessibility (DA) analysis of hybridization
# Guanjing Hu
# huguanjing@gmail.com
# https://wendellab.github.io/cottonMNase-seq

# module load py-deeptools r/3.5.0-py2-ufvuwmm bedtools2
# cd /work/LAS/jfw-lab/hugj2006/cottonLeaf
# R
options(scipen=999)
library(csaw)
library(edgeR)
library(ggplot2)

# wrapper
mergeResults =function (ranges, tab, tol, get.best = TRUE, merge.args = list(),
    combine.args = list(), best.args = list())
{
    merged <- do.call(mergeWindows, c(list(ranges, tol = tol), merge.args))
    combined <- do.call(combineTests, c(list(merged$id, tab), combine.args))
    output <- DataFrame(regions = I(merged$region), combined = I(combined))
    if (get.best) {
        output$best <- do.call(getBestTest, c(list(merged$id,tab), best.args))
    }
    metadata(output)$id <- merged$id
    output
}
overlapResults = function (ranges, tab, regions, get.best = TRUE, overlap.args = list(),    combine.args = list(), best.args = list())
{
    olap <- do.call(findOverlaps, c(list(query = regions, subject = ranges),
        overlap.args))
    combined <- do.call(combineOverlaps, c(list(olap, tab), combine.args))
    output <- DataFrame(regions = I(regions), combined = I(combined))
    if (get.best) {
        output$best <- do.call(getBestOverlaps, c(list(olap,tab), best.args))
    }
    metadata(output)$overlaps <- olap
    output
}
# merged <- mergeResults(filtered.data, resDi$table, tol=100, merge.args=list(max.width=5000))
sumTab<-function(resDi,resF,resH,med="")
{
    d<-summary(decideTests(resDi)) # FDR <0.05
    f<-summary(decideTests(resF)) # FDR <0.05
    h<-summary(decideTests(resH))
    tbl<-as.data.frame(cbind(d,f,h))
    tbl$med<-med
    return(tbl)
}
# generate MA plots
plotRes<-function(resDi,resF,resM,resHr,resPr,n=10000,out="out.pdf"){
    require("ggplot2")
    require("gridExtra")
    resDi$table$sig = factor(decideTests(resDi), levels=c(-1,0,1))
    resF$table$sig  = factor(decideTests(resF),levels=c(-1,0,1))
    resM$table$sig  = factor(decideTests(resM),levels=c(-1,0,1))
    resHr$table$sig  = factor(decideTests(resHr),levels=c(-1,0,1))
    resPr$table$sig  = factor(decideTests(resPr),levels=c(-1,0,1))
    select = sample(1:nrow(resDi$table),n)
    cvalue<-c("red","black","red");names(cvalue)=c(-1,0,1)
    pdf(out)
    p1<-ggplot(data=data.frame(resDi$table[select,]), aes(x = logCPM, y = logFC, col = sig)) + geom_point() + scale_color_manual(values =cvalue[levels(droplevels(resDi$table[select,]$sig))]) + geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + geom_hline(yintercept = 0) + labs(col = NULL) + theme_bw()+ theme(legend.position = "none")
    p2<-ggplot(data=data.frame(resF$table[select,]), aes(x = logCPM, y = logFC, col = sig)) + geom_point() + scale_color_manual(values =cvalue[levels(droplevels(resF$table[select,]$sig))]) + geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + geom_hline(yintercept = 0) + labs(col = NULL) + theme_bw()+ theme(legend.position = "none")
    p3<-ggplot(data=data.frame(resM$table[select,]), aes(x = logCPM, y = logFC, col = sig)) + geom_point() + scale_color_manual(values =cvalue[levels(droplevels(resM$table[select,]$sig))]) + geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + geom_hline(yintercept = 0) + labs(col = NULL) + theme_bw()+ theme(legend.position = "none")
    p4<-ggplot(data=data.frame(resHr$table[select,]), aes(x = logCPM, y = logFC, col = sig)) + geom_point() + scale_color_manual(values =cvalue[levels(droplevels(resHr$table[select,]$sig))]) + geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + geom_hline(yintercept = 0) + labs(col = NULL) + theme_bw()+ theme(legend.position = "none")
    p5<-ggplot(data=data.frame(resPr$table[select,]), aes(x = logCPM, y = logFC, col = sig)) + geom_point() + scale_color_manual(values =cvalue[levels(droplevels(resPr$table[select,]$sig))]) + geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + geom_hline(yintercept = 0) + labs(col = NULL) + theme_bw()+ theme(legend.position = "none")
    grid.arrange(p1, p2, p3, p4,p5,ncol = 3, nrow = 3)
    dev.off()
}
########

# The input BAM files were generated by mapping all MNase-seq reads against the F1 reference genome (concatenated by A2 and D5 referemce genome used; F2020=A2+D5), properly paired restricted, mapping quality filtered, and sorted using script "runHisat2.sh".

###################################
## From sorted bam to DB results ##
###################################
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/Qregulation_diffACRs/Fref")

# specify paired-end BAMs from A2, D5, F1
bam.files=grep("A6H_",list.files(patter="sort.bam$"),invert=T,value=T)
df=data.frame(BAM=bam.files, genome=substring(bam.files,1,1),digestion=substring(bam.files,3,3))
# Chr D01 blacklist
blacklist =GRanges("D01", IRanges(23162000, 23820000))

# input identified ACRs
load("../../Qregulation3/ACRsC.rdata")
D5=ACR.D5
seqlevels(D5) = gsub("^Chr","D",seqlevels(D5))
A2=ACR.A2
seqlevels(A2) = gsub("^Chr","A",seqlevels(A2))
F1a=ACR.F1a
seqlevels(F1a) = gsub("^Chr","A",seqlevels(F1a))
F1d=ACR.F1d
seqlevels(F1d) = gsub("^Chr","D",seqlevels(F1d))

# check overlapping
length(findOverlaps(F1d, D5))/length(D5) # 0.09980345
length(findOverlaps(F1d, D5))/length(F1d) #  0.08482602
length(findOverlaps(F1a, A2))/length(A2) # 0.1145313
length(findOverlaps(F1a, A2))/length(F1a) # 0.09501611

# define concensus peakset
peaks=list(A = union(A2,F1a),D=union(D5,F1d))


# separate A and D subgenomes for analysis
for(g in c("A","D"))
{
   ## 1. input sorted and indexed bam files
   use = df[df$genome=="F"|df$genome=="M"|df$genome==g,]
   use.chr <- c(paste0(g,"0", 1:9), paste0(g,10:13))
   print(use)
   print(use.chr)
   bams=as.character(use$BAM)
   # minimal MAPQ quality 20, paired-end
   param <- readParam(pe="both", discard=blacklist, restrict=use.chr);param
   
   ## 2.1. Count reads into ACRs
   peak.counts <- regionCounts(bams, peaks[[g]], param=param)
   # filter low abundance peaks
   peak.abundances <- aveLogCPM(asDGEList(peak.counts))
   # ......only use peaks logCPM > -3; few or no peaks should be removed; modify as desired
   keep <- peak.abundances > -3
   summary(keep)
   peak.counts.filt <- peak.counts[keep, ]
   
   ## 2.2. Count reads into 150 bp windows, based on fragment size
   win.data <- windowCounts(bams, ext=150, width=150, param=param)
   win.data
   head(assay(win.data))
   # Filtering out uninteresting regions.
   # simple restriction to filter window abundances at 5% lowe end
   abundances <- aveLogCPM(asDGEList(win.data))
   summary(abundances)
   keep.simple <- abundances > quantile(abundances,0.05)
   summary(keep.simple)
   # remove low variance
   logcounts <- cpm(asDGEList(win.data), log=TRUE)
   vars <- apply(logcounts, 1, var)
   keep.var <- vars > quantile(vars,0.25)
   summary(keep.var)
   # global enrichment to keep windows at least 1.5 fold of background signal
   bin.size <- 2000L
   binned <- windowCounts(bams, bin=TRUE, width=bin.size, param=param)
   filter.stat <- filterWindows(win.data, binned, type="global")
   min.fc=1.5
   pdf(paste0("checkBackground_",g,".pdf"))
   hist(filter.stat$back.abundances, xlab="Adjusted bin log-CPM", breaks=100, main="", col="grey80", xlim=c(min(filter.stat$back.abundances), 0))
   global.bg <- filter.stat$abundances - filter.stat$filter
   abline(v=global.bg[1], col="red", lwd=2)
   abline(v=global.bg[1]+log2(min.fc), col="blue", lwd=2)
   legend("topright", lwd=2, col=c('red', 'blue'),legend=c("Background", "Threshold"))
   dev.off()
   keep <- filter.stat$filter > log2(min.fc)
   summary(keep) # two stringent, losing too many data, not suitable for DNS-MNase-seq
   # filter used
   summary(keep.simple&keep.var)
   counts.local.filt <- win.data[keep.simple&keep.var,]
  
   ## 3. Calculating normalization factors.
   #  count BAM background bins (for TMM normalization)
   binned <- windowCounts(bams, bin=TRUE, width=10000, param=param)
   # --------------NORMALIZATION--------------------
   # method 1: ACR peaks only, TMM normalization based on binned counts
   peak.counts.tmm <- peak.counts.filt
   peak.counts.tmm <- normFactors(binned, se.out=peak.counts.tmm)
   # method 2: ACR peaks only, csaw loess-normalization
   peak.counts.loess <- peak.counts.filt
   peak.counts.loess <- normOffsets(peak.counts.loess, type="loess", se.out=TRUE)
   # from vignette: "For type="loess", a numeric matrix of the same dimensions as counts, containing the log-based offsets for use in GLM fitting.
   # method 3: csaw de novo peaks by local enrichment, TMM normalization based on binned counts
   counts.local.tmm <- counts.local.filt
   counts.local.tmm <- normFactors(binned, se.out=counts.local.tmm)
   # method 4: csaw de novo peaks by local enrichment, csaw loess-normalization
   counts.local.loess <- counts.local.filt
   counts.local.loess <- normOffsets(counts.local.loess, type="loess", se.out=TRUE)
   # from vignette: "For type="loess", a numeric matrix of the same dimensions as counts, containing the log-based offsets for use in GLM fitting."
   # --------------NORMALIZATION--------------------
   save(list=c("peak.counts.tmm","peak.counts.loess","counts.local.tmm","counts.local.loess"), file=paste0("DAnorm_",g,".rdata"))
   
   ## 4. DA analysis with 4 different methods
   for(med in c("peak.counts.tmm","peak.counts.loess","counts.local.tmm","counts.local.loess"))
   {
       # set working windows for the desired analysis
       working.windows <- get(med)
       y <- asDGEList(working.windows)
       sample<-gsub("_.*","",use$BAM)
       colnames(y$counts) <- sample
       rownames(y$samples) <- sample
       use$genome=gsub("A|D","P",use$genome)
       y$samples$group <- factor(paste0(use$genome,use$digestion))
       # setup design matrix; see edgeR manual for more information
       design <- model.matrix(~0+group, data=y$sample)
       # stabilize dispersion estimates with empirical bayes
       y <- estimateDisp(y, design)
       fit <- glmQLFit(y, design, robust=TRUE)
         # Estimate dispersion, plot e coefficient of biological variation (BCV) and quasi-likelihood (QL)dispersion
         # pdf(paste0("Dispersion_",g,".pdf"))
         # plotBCV(y)
         # plotQLDisp(fit)
         # dev.off()
       
       # DA analyses
       resDi <- glmQLFTest(fit, contrast=makeContrasts(groupPL-groupPH,levels=design))
       resF  <- glmQLFTest(fit, contrast=makeContrasts(groupFL-groupFH,levels=design))
       resM  <- glmQLFTest(fit, contrast=makeContrasts(groupML-groupMH,levels=design))
       resHr  <- glmQLFTest(fit, contrast=makeContrasts((groupFL-groupFH)-(groupPL-groupPH),levels=design))
       resPr  <- glmQLFTest(fit, contrast=makeContrasts((groupML-groupMH)-(groupFL-groupFH),levels=design))
       
       # MA plots and summary table
       # summary(decideTests(resDi)) # FDR <0.05
       tbl<-as.data.frame(cbind( summary(decideTests(resDi)), summary(decideTests(resF)), summary(decideTests(resM)), summary(decideTests(resHr)), summary(decideTests(resPr)) ))
       tbl$med<-med
       tbl$genome=g
       write.table(tbl,file=paste0("sigTable.",g,".",med,".txt"),sep="\t")
       plotRes(resDi,resF,resM,resHr,resPr,out=paste0("plotMA.",g,".",med,".pdf"))
       
       # for method 3 and 4, recalculate FDR by ACRs
       if(grepl("^counts",med))
       {
           for(i in c("resDi","resF","resM","resHr","resPr"))
           {
               print(i)
               res=get(i)
               olap <- overlapResults(working.windows, regions=peaks[[g]], tab=res$table)
               cx<-as.data.frame(table(olap$combined$direction,olap$combined$FDR<0.05))
               names(cx)<-c("direction","sig","number")
               cx$med<-med
               print(cx)
               cx$res<-i
               if(i=="resDi"){acx<-cx}else{acx<-rbind(acx,cx)}
               assign(i,olap)
           }
           write.table(acx,file=paste0("sigTable2.",g,".",med,".txt"),row.names=FALSE,sep="\t")
       }
      
      # save
       save(resDi,resF,resM,resHr,resPr, working.windows, file=paste0("DAres_",g,".",med,".rdata"))
    }
}

##########################################
## Compare DB results by 4 norm methods ##
##########################################
# put 4 method results together to inspect
library(reshape2)
fl<-list.files(pattern="sigTable[.].*peak")
for(f in fl){
    x<-read.table(f,header=T, sep="\t")
    y<-cbind(c("down","up"), x[c("Down","Up"),])
    names(y)=c("direction","P","F","M","Hr","Pr","med","genome")
    if(f==fl[1]){res<-y}else{res<-rbind(res,y)}
}
fl<-list.files(pattern="sigTable2")
for(f in fl){
    acx<-read.table(f,header=T, sep="\t")
    ss=dcast(acx,direction+sig+med~res,value.var="number")
    ss[is.na(ss)]=0
    #ss=cbind(acx[acx$res=="resDi",1:3],acx[acx$res=="resF",3],acx[acx$res=="resH",3:4])
    sss=ss[ss$sig&ss$direction!="mixed",c("direction","resDi","resF","resM","resHr","resPr","med")]
    sss$genome=unlist(strsplit(f,"[.]"))[2]
    names(sss)=c("direction","P","F","M","Hr","Pr","med","genome")
    res<-rbind(res,sss)
}
# load peaks,
load("../../Qregulation3/ACRsC.rdata")
ACR.D5
D5=ACR.D5
seqlevels(D5) = gsub("^Chr","D",seqlevels(D5))
A2=ACR.A2
seqlevels(A2) = gsub("^Chr","A",seqlevels(A2))
F1a=ACR.F1a
seqlevels(F1a) = gsub("^Chr","A",seqlevels(F1a))
F1d=ACR.F1d
seqlevels(F1d) = gsub("^Chr","D",seqlevels(F1d))
# define concensus peakset
peaks=list(A = union(A2,F1a),D=union(D5,F1d))
# examine percentage of ACRs found up in Diploid and F1
exp=c(length(peaks$A),length(peaks$D))
names(exp)=c("A","D");exp
#A      D
#619384 396174
res$acrN=exp[as.character(res$genome)]
res$Fperc=res$F/res$acrN
res$Pperc=res$P/res$acrN
res$Mperc=res$M/res$acrN
write.table(res, file="methodComp.txt",row.names=FALSE,sep="\t")

library(ggplot2)
library("gridExtra")
x<-read.table("methodComp.txt",head=T,sep="\t")
ss<- theme_bw()+ theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
p1<-ggplot(data=x,aes(x=med,y=P,group=med,fill=direction))+geom_bar(stat="identity") + facet_grid(genome~.)
p2<-ggplot(data=x,aes(x=med,y=F,group=med,fill=direction))+geom_bar(stat="identity") + facet_grid(genome~.)
p3<-ggplot(data=x,aes(x=med,y=M,group=med,fill=direction))+geom_bar(stat="identity") + facet_grid(genome~.)
p4<-ggplot(data=x,aes(x=med,y=Hr,group=med,fill=direction))+geom_bar(stat="identity") + facet_grid(genome~.)
p5<-ggplot(data=x,aes(x=med,y=Pr,group=med,fill=direction))+geom_bar(stat="identity") + facet_grid(genome~.)
p6<-ggplot(data=x,aes(x=med,y=Fperc,group=med,fill=direction))+geom_bar(stat="identity") + facet_grid(genome~.)
p7<-ggplot(data=x,aes(x=med,y=Pperc,group=med,fill=direction))+geom_bar(stat="identity") + facet_grid(genome~.)
p8<-ggplot(data=x,aes(x=med,y=Mperc,group=med,fill=direction))+geom_bar(stat="identity") + facet_grid(genome~.)
pdf("methodComp.pdf")
p1
p2
p3
grid.arrange(p1+ss, p2+ss, p3+ss, ncol = 3, nrow = 1)
p4
p5
grid.arrange(p4+ss, p5+ss, ncol = 3, nrow = 1)
p6
p7
p8
grid.arrange(p6+ss, p7+ss, p8+ss, ncol = 3, nrow = 1)
dev.off()
# based on p1~3, peak.counts.tmm > counts.local.tmm > counts.local.loess >> peak.counts.loess in finding ACRs
# Hr: all methods agree on that few changes (<30 ACRs) by hybridization, somehow more in D than A genome; more decrease then increase in accessbility

###########################################################
## De novo DA windows analysis: genome-wide and promoter ##
###########################################################
options(scipen=999)
library(csaw)
library(edgeR)
library(ggplot2)
library(GenomicFeatures)

setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/Qregulation_diffACRs/Fref")
bl =GRanges("D01", IRanges(23162000, 23820000))

# input identified ACRs
load("../../Qregulation3/ACRsC.rdata"); ACR.D5
D5=ACR.D5
seqlevels(D5) = gsub("^Chr","D",seqlevels(D5))
h=findOverlaps(D5,bl); length(h) # 0 in black region
A2=ACR.A2
seqlevels(A2) = gsub("^Chr","A",seqlevels(A2))
F1a=ACR.F1a
seqlevels(F1a) = gsub("^Chr","A",seqlevels(F1a))
F1d=ACR.F1d
seqlevels(F1d) = gsub("^Chr","D",seqlevels(F1d))
h=findOverlaps(F1d,bl); length(h) # 1357 in black region
F1d=F1d[-queryHits(h)] # from 224483 to 223126

# define concensus peakset
peaks=list(A = union(A2,F1a),D=union(D5,F1d))

# input promoters: 2kb upstream of TSS
### import genome annotation
txdb <- loadDb("../../refGenomes/txdb.F2020.sqlite")
gns= genes(txdb)
pms = promoters(txdb, upstream=1000, downstream=0) # 120545
pms= pms[grep("scaffold|Contig",seqnames(pms),invert=T)]# 118590
pms.a = pms[grep("A",seqnames(pms)),]
pms.d = pms[grepl("D",seqnames(pms))&grepl("[.]1$",names(pms)),]
pms=list(A=pms.a,D=pms.d)

# specify paired-end BAMs from A2, D5, F1
bam.files=grep("A6H_",list.files(patter="sort.bam$"),invert=T,value=T)
df=data.frame(BAM=bam.files, genome=substring(bam.files,1,1),digestion=substring(bam.files,3,3))
# Chr D01 blacklist

# separate A and D subgenomes for analysis
for(g in c("A","D"))
{
    use = df[df$genome=="F"|df$genome=="M"|df$genome==g,]
    load( paste0("DAnorm_",g,".rdata") )

    # set working windows for the desired analysis
    med="counts.local.tmm"
    working.windows <- get(med)
    y <- asDGEList(working.windows)
    sample<-gsub("_.*","",use$BAM)
    colnames(y$counts) <- sample
    rownames(y$samples) <- sample
    use$genome=gsub("A|D","P",use$genome)
    y$samples$group <- factor(paste0(use$genome,use$digestion))
    # setup design matrix; see edgeR manual for more information
    design <- model.matrix(~0+group, data=y$sample)
    # stabilize dispersion estimates with empirical bayes
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
      # Estimate dispersion, plot e coefficient of biological variation (BCV) and quasi-likelihood (QL)dispersion
      # pdf(paste0("Dispersion_",g,".pdf"))
      # plotBCV(y)
      # plotQLDisp(fit)
      # dev.off()
    
    # DA analyses
    resDi <- glmQLFTest(fit, contrast=makeContrasts(groupPL-groupPH,levels=design))
    resF  <- glmQLFTest(fit, contrast=makeContrasts(groupFL-groupFH,levels=design))
    resM  <- glmQLFTest(fit, contrast=makeContrasts(groupML-groupMH,levels=design))
    resHr  <- glmQLFTest(fit, contrast=makeContrasts((groupFL-groupFH)-(groupPL-groupPH),levels=design))
    resPr  <- glmQLFTest(fit, contrast=makeContrasts((groupML-groupMH)-(groupFL-groupFH),levels=design))
    
    for(i in c("resDi","resF","resM","resHr","resPr"))
    {
        print(i)
        res=get(i)
        # merging windows, Correcting for multiple testing.
        merged <- mergeResults(working.windows, res$table, tol=200L)
        # overlap to ACRs, Correcting for multiple testing.
        olACR <- overlapResults(working.windows, tab=res$table, regions=peaks[[g]])
        # overlap to promoters, Correcting for multiple testing.
        olP <- overlapResults(working.windows, tab=res$table, regions=pms[[g]])
        save(merged, olACR, olP, file=paste0(i,"_",g,".",med,".rdata"))
    }
}


# Examine results

sigRes=function(res, label=NULL){
    require(reshape2)
    names(res) # "regions"  "combined" "best"
    res0=data.frame(Var1=c("sigN","sigBp","sigN","sigBp"),Var2=c("up","up","down","down"),value=0)
    # combined results
    sig=which(res$combined$FDR<0.05)
    if(length(sig)>0){
        sigN=tapply(width(res$regions[sig]),res$combined$direction[sig], FUN=length)
        sigBp=tapply(width(res$regions[sig]),res$combined$direction[sig], FUN=sum)
        Cres=melt(rbind(sigN,sigBp))
    }else{ Cres= res0}
    Cres$FDR="combined"
    # best results
    sig=which(res$best$FDR<0.05)
    if(length(sig)>0){
        res$best$direction = ifelse(res$best$logFC>0, "up","down")
        sigN=tapply(width(res$regions[sig]),res$best$direction[sig], FUN=length)
        sigBp=tapply(width(res$regions[sig]),res$best$direction[sig], FUN=sum)
        Bres=melt(rbind(sigN,sigBp))
    }else{ Bres= res0}
    Bres$FDR="best"
    lres=rbind(Cres,Bres)
    if(!is.na(label)){lres$label=label}
    return(lres)
}
fl<-list.files(pattern="^res")
for(f in fl){
    load(f)
    print(f)
    flag=unlist(strsplit(f,"_|[.]"))
    sr<-rbind(sigRes(merged,"SlidingW"), sigRes(olACR,"ACR"), sigRes(olP,"Promoter"))
    sr$test=flag[1]
    sr$genome=flag[2]
    # print(sumRes)
    if(f==fl[1]){sumtbl=sr}else{sumtbl=rbind(sumtbl,sr)}
}
write.table(sumtbl,file="resNnBp.Fref.txt",row.names=FALSE,sep="\t")
sumtbl

x<-read.table("resNnBp.Fref.txt",header=T,sep="\t")
library(ggplot2)
pdf("resNnBp.Fref.pdf")
ss<- theme_bw()+ theme(legend.position = "right")+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
# resHr vs resPr, separate N and Bp.
ggplot(data=x[x$test%in%c("resHr","resPr")&x$Var1=="sigN",],aes(x=test,y=value,group=test,fill=Var2))+geom_bar(stat="identity")+ facet_grid(genome~label+FDR)+ss
ggplot(data=x[x$test%in%c("resHr","resPr")&x$Var1=="sigBp",],aes(x=test,y=value,group=test,fill=Var2))+geom_bar(stat="identity")+ facet_grid(genome~label+FDR)+ss
# resDi, F, M
ggplot(data=x[x$test%in%c("resDi","resF","resM")&x$Var1=="sigN",], aes(x=test,y=value,group=test,fill=Var2))+geom_bar(stat="identity")+ facet_grid(genome~label+FDR)+ss
ggplot(data=x[x$test%in%c("resDi","resF","resM")&x$Var1=="sigBp",], aes(x=test,y=value,group=test,fill=Var2))+geom_bar(stat="identity")+ facet_grid(genome~label+FDR)+ss
dev.off()
# CONCLUSION: 1.Pr > Hr increased accessibility for ACR, genome-wid, and 1kb promoter, what about gene body and other locations, eg. TE? 2. best and combined results are quite consistent

# --- ABANDON BELOW


#################################
## Compare DB results for MSFs ##
#################################
# wrapper
mergeResults =function (ranges, tab, tol, get.best = TRUE, merge.args = list(),
    combine.args = list(), best.args = list())
{
    merged <- do.call(mergeWindows, c(list(ranges, tol = tol), merge.args))
    combined <- do.call(combineTests, c(list(merged$id, tab), combine.args))
    output <- DataFrame(regions = I(merged$region), combined = I(combined))
    if (get.best) {
        output$best <- do.call(getBestTest, c(list(merged$id,tab), best.args))
    }
    metadata(output)$id <- merged$id
    output
}
overlapResults = function (ranges, tab, regions, get.best = TRUE, overlap.args = list(),    combine.args = list(), best.args = list())
{
    olap <- do.call(findOverlaps, c(list(query = regions, subject = ranges),
        overlap.args))
    combined <- do.call(combineOverlaps, c(list(olap, tab), combine.args))
    output <- DataFrame(regions = I(regions), combined = I(combined))
    if (get.best) {
        output$best <- do.call(getBestOverlaps, c(list(olap,tab), best.args))
    }
    metadata(output)$overlaps <- olap
    output
}
# merged <- mergeResults(filtered.data, resDi$table, tol=100, merge.args=list(max.width=5000))

### how do DB results compare to ACRs? Quite Consistent
load("../../Qregulation3/ACRsC.rdata")
D5=ACR.D5
seqlevels(D5) = gsub("^Chr","D",seqlevels(D5))
A2=ACR.A2
seqlevels(A2) = gsub("^Chr","A",seqlevels(A2))
F1a=ACR.F1a
seqlevels(F1a) = gsub("^Chr","A",seqlevels(F1a))
F1d=ACR.F1d
seqlevels(F1d) = gsub("^Chr","D",seqlevels(F1d))

load("../Fref_old/DBres_D.rdata")
#-----------------------------------
# D: resDi vs MSF.D5
olap <- overlapResults(filtered.data, regions=D5, tab=resDi$table)
Dd=olap$combined
table(olap$combined$direction,olap$combined$FDR<0.05)
#         FALSE  TRUE
#down     68     9
#mixed    83     1
#up    18535 42780   69.0% MSFs supported by csaw
#-----------------------------------
# D: resF vs MSF.F1d
olap <- overlapResults(filtered.data, regions=F1d, tab=resF$table)
Fd<-olap$combined
table(olap$combined$direction,olap$combined$FDR<0.05)
#      FALSE  TRUE
#down     18     0
#mixed    55     0
#up    16397 43498, 72.2% F MSFs supported
#-----------------------------------
# D: resH vs MSF.D5
olap <- overlapResults(filtered.data, regions=D5, tab=resH$table)
table(olap$combined$direction,olap$combined$FDR<0.05)
#      FALSE  TRUE
#down  47181    23
#mixed 10751     0
#up     3520     1    only 24 ACRs affected by hybridization
Hd<-olap$combined
#-----------------------------------
# D: resH vs MSF.F1d
olap <- overlapResults(filtered.data, regions=F1d, tab=resH$table)
Hfd<-olap$combined
table(olap$combined$direction,olap$combined$FDR<0.05)
#       FALSE
# down   5192
# mixed 13301
# up    41475     zero ACRs affected by hybridization
olapD <- DataFrame(regions = I(MSF.D5),  D=I(Dd), H=I(Hd))
olapFd<- DataFrame(regions = I(MSF.F1d), F=I(Fd), H=I(Hfd))

#-----------------------------------
# combined results
df=data.frame(regions =c("DNS_diploid","DNS_F1","F1-diploid"),sizeUp=NA,noUp=NA,sizeDown=NA,noDown=NA)
rownames(df)=c("resDi","resF","resH")
for(i in c("resDi","resF","resH")){
    print(i)
    res=get(i)
    # keep p-value significant regions, combined, FDR<0.05, DB regions
    ff = (res$table$PValue<0.05)
    merged <- mergeResults(filtered.data[ff,], res$table[ff,], tol=100, merge.args=list(max.width=5000))
    gr = merged$regions[merged$combined$FDR<0.05,]
    mcols(gr) =cbind(mcols(gr),merged$combined[merged$combined$FDR<0.05,])
    df[i,"sizeUp"] =sum(width(gr[gr$direction=="up"]))/10^6
    df[i,"noUp"] =length(gr[gr$direction=="up"])
    df[i,"sizeDown"] =sum(width(gr[gr$direction=="down"]))/10^6
    df[i,"noDown"] =length(gr[gr$direction=="down"])
    assign(paste0("gr_",i),gr)
}
print(df)
#          regions   sizeUp   noUp sizeDown noDown
#resDi DNS_diploid 89.47695 265050 68.10773 217958
#resF       DNS_F1 95.67370 268178 65.60778 207936
#resH   F1-diploid 20.62405  96885 19.09055  90140
dfD=df
combinedD = list(D=gr_resDi,F=gr_resF,H=gr_resH)
save(dfD,combinedD,olapD,olapFd, file="combinedRes.rdata")

load("DBres_A.rdata")
#-----------------------------------
# A: resDi vs MSF.A2
olap <- overlapResults(filtered.data, regions=A2, tab=resDi$table)
Da<-olap$combined
table(olap$combined$direction,olap$combined$FDR<0.05)
#       FALSE  TRUE
#down    568    12
#mixed  1480     0
#up    64551 59872, 46.5% of A2 MSFs supported, probably due to lower quality of A6H and A6L
#-----------------------------------
# A: resH vs MSF.A2
olap <- overlapResults(filtered.data, regions=A2, tab=resH$table)
Ha<-olap$combined
table(olap$combined$direction,olap$combined$FDR<0.05)
#        FALSE
# down  86988
# mixed 28453
# up    11042
#-----------------------------------
# D: resF vs MSF.F1a
olap <- overlapResults(filtered.data, regions=F1a, tab=resF$table)
Fa<-olap$combined
table(olap$combined$direction,olap$combined$FDR<0.05)
#         FALSE  TRUE
# down     29     0
# mixed    55     0
# up    36526 90900  , 70.8% F1a MSFs supported
#-----------------------------------
# D: resH vs MSF.F1a
olap <- overlapResults(filtered.data, regions=F1a, tab=resH$table)
Hfa<-olap$combined
table(olap$combined$direction,olap$combined$FDR<0.05)
#       FALSE
# down  11547
# mixed 30308
# up    85655     zero ACRs affected by hybridization
olapA <- DataFrame(regions = I(MSF.A2),  D=I(Da), H=I(Ha))
olapFa<- DataFrame(regions = I(MSF.F1a), F=I(Fa), H=I(Hfa))
#-----------------------------------

# combined results
df=data.frame(regions =c("DNS_diploid","DNS_F1","F1-diploid"),sizeUp=NA,noUp=NA,sizeDown=NA,noDown=NA)
rownames(df)=c("resDi","resF","resH")
for(i in c("resDi","resF","resH")){
    print(i)
    res=get(i)
    # keep p-value significant regions, combined, FDR<0.05, DB regions
    ff = (res$table$PValue<0.05)
    merged <- mergeResults(filtered.data[ff,], res$table[ff,], tol=100, merge.args=list(max.width=5000))
    gr = merged$regions[merged$combined$FDR<0.05,]
    mcols(gr) =cbind(mcols(gr),merged$combined[merged$combined$FDR<0.05,])
    df[i,"sizeUp"] =sum(width(gr[gr$direction=="up"]))/10^6
    df[i,"noUp"] =length(gr[gr$direction=="up"])
    df[i,"sizeDown"] =sum(width(gr[gr$direction=="down"]))/10^6
    df[i,"noDown"] =length(gr[gr$direction=="down"])
    assign(paste0("gr_",i),gr)
}
print(df)
#          regions   sizeUp   noUp sizeDown noDown
#resDi DNS_diploid 180.5343 425852 108.2926 355316
#resF       DNS_F1 175.2772 461616  96.9000 332499
#resH   F1-diploid  24.2245 117443  29.1255 136414
dfA=df
combinedA = list(D=gr_resDi,F=gr_resF,H=gr_resH)

# DB results equivalent to iSeg for diploids and F1
save(dfD,combinedD,olapD,olapFd,dfA,combinedA,olapA,olapFa, file="DBres_out.rdata")


#########################
## Annotate DB results ##
#########################
load("DBres_out.rdata")
setwd("../../")
# annotateACRs function
source("FUN.acr.r")
# load chipseeker annotation results
library(ChIPseeker)
load("Qregulation3/MSFpg.rdata")->l;l # "anL" "df"
names(anL) # "A2"     "D5"     "F1a"    "F1d"    "F1.At"  "F1.Dt"  "AD1a"   "AD1d" "AD1.At" "AD1.Dt
# check correspondance, TRUE
#table(ranges(anL$A2@anno) ==ranges(olapA$regions[olapA$an$type!="Distal"]))
#table(ranges(anL$D5@anno) ==ranges(olapD$regions[olapD$an$type!="Distal"]))


## import genome annotation
txdb <- loadDb("refGenomes/txdb.A2WHU.sqlite")
gns= genes(txdb)
gns.A2= gns[grep("Contig",seqnames(gns),invert=T)] # 41739
txdb <- loadDb("refGenomes/txdb.D5.sqlite")
gns= genes(txdb)
gns.D5= gns[grep("scaffold",seqnames(gns),invert=T)]

## check DB with annotation
olap=olapA; gns=gns.A2; fine = anL$A2@anno
an= annotateACRs(olap$regions,gns, distance=2000)
an$fine="-"
an$fine[an$type!="Distal"] = fine$annotation2
olap$an = I(an[,12:20])
olapA=olap
table(ranges(anL$A2@anno) ==ranges(olapA$regions[olapA$an$type!="Distal"]))
with(olap,table(an$type,H$direction, H$FDR<0.05))
#, ,  = FALSE
#          down mixed    up
#Distal   72096 25209 10030
#Genic     5288  1064   299
#Proximal  9604  2180   713
with(olap,table(an$fine,H$direction, H$FDR<0.05))
#, ,  = FALSE
#                    down mixed    up
#-                  72096 25209 10030
#Distal Intergenic   1292   308   119
#Downstream (<1kb)    184    38    12
#Downstream (1-2kb)   109    25    10
#Downstream (2-3kb)    88    16     3
#Exon                 496    82    26
#Intron              1451   346   106
#Promoter (<=1kb)    5174  1100   318
#Promoter (1-2kb)    3972   887   288
#Promoter (2-3kb)    2126   442   130

#---------------------------------------------------
olap=olapD; gns=gns.D5; fine = anL$D5@anno
an= annotateACRs(olap$regions,gns, distance=2000)
an$fine="-"
an$fine[an$type!="Distal"] = fine$annotation2
olap$an = I(an[,12:20])
olapD=olap
table(ranges(anL$D5@anno) ==ranges(olapD$regions[olapD$an$type!="Distal"]))
with(olap,table(an$type,H$direction, H$FDR<0.05))
#, ,  = FALSE
#            down mixed    up
#  Distal   36426  8456  2783
#  Genic     4490   771   216
#  Proximal  6265  1524   521
#, ,  = TRUE
#            down mixed    up
#  Distal      12     0     0
#  Genic        7     0     1
#  Proximal     4     0     0
with(olap,table(an$fine,H$direction, H$FDR<0.05))
#, ,  = FALSE
#                      down mixed    up
#  -                  36426  8456  2783
#  Distal Intergenic    995   257    73
#  Downstream (<1kb)    104    32    10
#  Downstream (1-2kb)    58    10     1
#  Downstream (2-3kb)    41    13     7
#  Exon                 667    94    29
#  Intron               786   151    35
#  Promoter (<=1kb)    3642   824   266
#  Promoter (1-2kb)    2890   609   225
#  Promoter (2-3kb)    1572   305    91
#, ,  = TRUE
#                      down mixed    up
#  -                     12     0     0
#  Distal Intergenic      0     0     0
#  Downstream (<1kb)      0     0     0
#  Downstream (1-2kb)     0     0     0
#  Downstream (2-3kb)     0     0     0
#  Exon                   1     0     0
#  Intron                 0     0     0
#  Promoter (<=1kb)       8     0     1
#  Promoter (1-2kb)       1     0     0
#  Promoter (2-3kb)       1     0     0
save(dfD,combinedD,olapD,olapFd,dfA,combinedA,olapA,olapFa, file="DBres_out.rdata")

###############################
## Direct comparison of MSFs ##
###############################
load("Qregulation3/ACRs_MSF.rdata")->l;l
length(h<-findOverlaps(MSF.A2,MSF.F1a)) # 12516
length(h)/length(MSF.A2)
length(h)/length(MSF.F1a) #9.7% shared peaks between A2 and F1a
length(h<-findOverlaps(MSF.D5,MSF.F1d)) # 4742
length(h)/length(MSF.D5) #7.6%
length(h)/length(MSF.F1d) #7.8% shared peaks between D5 and F1d

#################
## Conclusions ##
#################
# direction peak overlapping leads to only a small number of shared ACRs, mainly due to peak calling parameters variation
# CSAW approaches should be more reliable and rather stringent, leading to a small amount of differentiation by  hybridization

q("no")