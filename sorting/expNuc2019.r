###############################################################
## Inspect chromatin profiles in association with expression ##
###############################################################
library(GenomicRanges)
#library(RColorBrewer)
#library(ggbio)
library(GenomicFeatures)
library(genomation)
options(scipen=999)

# load annotation results of DNS segment annotation
l=load("isegRes/SegAnnotation/dns/genomicAnnotation.byType.rdata")
# "sampleInfo"       "sl"               "sl.nonzero"       "Sums" "Neg_peakAnnoList" "Pos_peakAnnoList"
sampleInfo$txDB=gsub("../../","",as.character(sampleInfo$txDB))
sampleInfo$inputBG=gsub("../","isegRes/",as.character(sampleInfo$inputBG))
sampleInfo$fullBW = c("aggregation/A*n*full.bw","aggregation/combo/D*full.bw","aggregation/combo/F*full.bw","aggregation/combo/M*full.bw")
sampleInfo

# quantile expression results from individual mapping
l=load("RNAseq/Ranalysis/expression.Indiv.rdata");l

# FUN - prepare quantile BED files
writeGRtoBED = function(gr, group="quantile",out="")
{
    df <- data.frame(seqnames=seqnames(gr), starts=start(gr)-1, ends=end(gr), names=c(rep(".", length(gr))),scores=c(rep(".", length(gr))),strands=strand(gr))
    if(!is.null(group)){
        id = gr$gene_id
        groupIdx=expr[id,group]
        dfL = split(df,groupIdx)
        for(j in 1:length(dfL))
        {
            write.table(dfL[[j]], file=paste0(out,".Q",names(dfL)[j],".bed"), quote=F, sep="\t", row.names=F, col.names=F)
        }
    }else{
        write.table(df, file=paste0(out,".bed"), quote=F, sep="\t", row.names=F, col.names=F)
    }
}

# write quantile gene annotation to bed
for(i in 1:4)
{
    # genomic features
    txdb = loadDb(as.character(sampleInfo$txDB[i]))
    genome=as.character( sampleInfo$genome[i] )
    print(genome)
    genes.gr =  genes(txdb)
    length(genes.gr)
    genes.gr=genes.gr[grep("^Chr|^A|^D",seqnames(genes.gr)),]
    length(genes.gr)
    expr=get(genome)$rpkm
    if(i %in% 1:2){
        writeGRtoBED(genes.gr, group="quantile",out=paste0("expNuc/", genome))
    }
    if(i %in% 3:4){
        writeGRtoBED(genes.gr[grep("^A",seqnames(genes.gr)),], group="quantile",out=paste0("expNuc/", genome,".A"))
        writeGRtoBED(genes.gr[grep("^D",seqnames(genes.gr)),], group="quantile",out=paste0("expNuc/", genome,".D"))
    }
}
list.files("expNuc/",pattern="Q*bed")


## prep deeptools script and run

# A2, D5, F1.At, F1.Dt, AD1.At, AD1.Dt each H/L/D with expression Q0-Q4
fileN = "expNuc.quantile.sh"
cat("module load py-deeptools", file=fileN,sep="\n")
cat("# A2, D5, F1.At, F1.Dt, AD1.At, AD1.Dt each H/L/D with expression Q0-Q4", file=fileN,sep="\n",append=TRUE)
for(i in 1:4)
{
    genome=as.character( sampleInfo$genome[i] )
    print(genome)
    if(i %in% 1:2){
        cmd = paste0("#",genome,"\ncomputeMatrix scale-regions -S ",sampleInfo$fullBW[i]," -R expNuc/",genome,"*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros\nplotHeatmap -m scaleGeneBody.gz -out expNuc/plotHeatmap.quantile.",genome,".png")
        cat(cmd, file=fileN,sep="\n",append=TRUE)
    }
    if(i %in% 3:4){
        cmd = paste0("#",genome,"\ncomputeMatrix scale-regions -S ",sampleInfo$fullBW[i]," -R expNuc/",genome,".A*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros\nplotHeatmap -m scaleGeneBody.gz -out expNuc/plotHeatmap.quantile.",genome,".A.png")
        cat(cmd, file=fileN,sep="\n",append=TRUE)
        cmd = paste0("computeMatrix scale-regions -S ",sampleInfo$fullBW[i]," -R expNuc/",genome,".D*.bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros\nplotHeatmap -m scaleGeneBody.gz -out expNuc/plotHeatmap.quantile.",genome,".D.png")
        cat(cmd, file=fileN,sep="\n",append=TRUE)
    }
}
system(paste0("cat ",fileN))
system(paste0("bash ",fileN))

# At vs Dt within F1 and AD1 for H/L/D
fileN = "expNuc.AvsD.sh"
cat("module load py-deeptools", file=fileN,sep="\n")
cat("# At vs Dt within F1 and AD1 for H/L/D", file=fileN,sep="\n",append=TRUE)
for(i in 3:4)
{
    genome=as.character( sampleInfo$genome[i] )
    print(genome)
    cmd = paste0("#",genome,"\ncomputeMatrix scale-regions -S ",sampleInfo$fullBW[i]," -R refGenomes/",genome,".gene*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros\nplotHeatmap -m scaleGeneBody.gz -out expNuc/plotHeatmap.AvsD.",genome,".png\nplotProfile -m scaleGeneBody.gz --perGroup -out expNuc/plotProfileGroup.AvsD.",genome,".png --outFileNameData expNuc/plotProfile.AvsD.",genome,".tab")
    cat(cmd, file=fileN,sep="\n",append=TRUE)
}
system(paste0("cat ",fileN))
system(paste0("bash ",fileN))

# note that direct comparison of A2 and D5 profiles, and A2 vs F1_At and D5 vs F1_Dt are easier with R package genomation, conducted by aggreNvisual.r


#########################################################################################
## Inspect ortho-homoeolog quadruplets for the presence of MSFs in 1kb promoter region ##
#########################################################################################

l = load("isegRes/SegAnnotation/dns/genomicAnnotation.byType.rdata")
## collect gene containing MSFs in 1kb promoter
getGenes=function(csAnno, feature="Promoter (<=1kb)")
{
    anno = csAnno@anno
    genes = anno$geneId[anno$annotation==feature]
    return(genes)
}
geneL = lapply(Pos_peakAnnoList, getGenes)
str(geneL)
# List of 6
# $ A2    : chr [1:12485] "Ga01G0007" "Ga01G0012" "Ga01G0017" "Ga01G0018" ...
# $ D5    : chr [1:10662] "Gorai.001G001000" "Gorai.001G002500" "Gorai.001G002600" "Gorai.001G003000" ...
# $ F1.At : chr [1:8771] "Ga01G0003" "Ga01G0005" "Ga01G0005" "Ga01G0005" ...
# $ F1.Dt : chr [1:10526] "Gorai.001G000500" "Gorai.001G001300" "Gorai.001G002400" "Gorai.001G002900" ...
# $ AD1.At: chr [1:18784] "Gohir.A01G000100" "Gohir.A01G000300" "Gohir.A01G000300" "Gohir.A01G000500" ...
# $ AD1.Dt: chr [1:17881] "Gohir.D01G000400" "Gohir.D01G000700" "Gohir.D01G000800" "Gohir.D01G000900" ...
#----- more genes were detected in AD1

## load 20362 ortholog quadruplets derived from Justin's orthoMCL results version3_081618
ogQ<-read.table("orthohomoeologQuadruplets101218.txt", sep="\t", header=TRUE)
head(ogQ)

# correspondece
og6 = data.frame(A2=ogQ$A2, D5=ogQ$D5, F1.At=ogQ$A2, F1.Dt=ogQ$D5, AD1.At=ogQ$At, AD1.Dt=ogQ$Dt)
rownames(og6)=ogQ$groupID
head(og6)
# within in each quadruplets, which contains MSFs and which doesn't
paMSF = data.frame( A2=og6$A2 %in% geneL$A2, D5=og6$D5 %in% geneL$D5, F1.At=og6$F1.At %in% geneL$F1.At, F1.Dt=og6$F1.Dt %in% geneL$F1.Dt, AD1.At=og6$AD1.At %in% geneL$AD1.At, AD1.Dt=og6$AD1.Dt %in% geneL$AD1.Dt)
rownames(paMSF) = rownames(og6)
head(paMSF)
# convert to 6 digit string code
paMSFn = apply(apply(paMSF,2,as.numeric), 1, function(x)paste0(x,collapse=""))
# examine the pattern of p/a
ddply(as.data.frame(paMSF[,-1]),.(A2,D5,F1.At,F1.Dt,AD1.At,AD1.Dt),nrow)
table(paMSFn)
# gene expression, use RPKM, log2mean
l=load("RNAseq/Ranalysis/expression.rdata")
l # "A2"  "D5"  "F1"  "AD1"
expr = data.frame(A2=A2$rpkm[og6$A2,"log2mean"], D5=D5$rpkm[og6$D5,"log2mean"], F1.At=F1$rpkm[og6$F1.At,"log2mean"], F1.Dt=F1$rpkm[og6$F1.Dt,"log2mean"], AD1.At=AD1$rpkm[og6$AD1.At,"log2mean"],AD1.Dt=F1$rpkm[og6$AD1.Dt,"log2mean"])
rownames(expr) = ogQ$groupID
# save
info = "og6 - ortho-homoeolog gene group index, paMSF - presence of 1kb promoter MSFs, paMSFn - 6 digit code for presence of 1kb promoter MSFs, expr - log2 mean rpkm+1"
save(info, og6, paMSF, paMSFn, expr, file="expNuc/ogPromoterMSFs.rdata")

## heatmap
library(RColorBrewer)
library(ComplexHeatmap)
library(heatmaps)
library(circlize)
pdf("expNuc/ogPromoterMSFs.pdf")
####### simple heatmap of total expression
select=1:500
mat = as.matrix(expr[select,])
mat_scaled = t(apply(mat, 1, scale))
#heatmap.2(mat, scale = "none", Colv=FALSE, dendrogram="row",col = bluered(100), trace = "none", density.info = "none", margins=c(8,5))
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
# heatmap.2(mat, scale = "none", Colv=FALSE, dendrogram="row",col = hmcol, trace = "none", margins=c(8,5))
####### ComplexHeatmap
ha = HeatmapAnnotation(df = data.frame(ploidy =c("diploid","diploid","F1","F1","AD1","AD1"),genome=c("A","D","A","D","A","D") ), col = list(ploidy=structure(brewer.pal(3, "Set2"), names = c("diploid","F1","AD1")),genome=c("A"="pink","D"="royalblue")))
#     present absence profiles, use default redBlue color
pa = as.matrix(apply(paMSF[select,],2,as.numeric))
hm= Heatmap(pa, name="MSF", col = c("0" = "white", "1" = "blue"),  cluster_column = FALSE, top_annotation = ha, top_annotation_height = unit(4, "mm"),column_title = "Promoter MSFs", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = T)
#     expression log2rpkm profiles, use default redBlue color
he.scaled=Heatmap(mat_scaled, name = "scaled log2rpkm", cluster_column = FALSE, top_annotation_height = unit(4, "mm"),column_title = "Scaled Expression", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = TRUE)
#     expression log2rpkm profiles, use default redBlue color
he=Heatmap(mat, name = "log2rpkm",col = hmcol, cluster_column = FALSE, top_annotation = ha,  top_annotation_height = unit(4, "mm"),column_title = "Expression", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = TRUE)
he+he.scaled+hm

df = as.data.frame(table(paMSFn))
for(i in 1:nrow(df)){
    select=(paMSFn==df$paMSFn[i])
    mat = as.matrix(expr[select,])
    mat_scaled = t(apply(mat, 1, scale))
    colnames(mat_scaled)=colnames(mat)
    pa = as.matrix(apply(paMSF[select,],2,as.numeric))
    
    #     present absence profiles, use default redBlue color
    hm= Heatmap(pa, name="MSF", col = c("0" = "white", "1" = "blue"),  cluster_column = FALSE, top_annotation = ha, top_annotation_height = unit(4, "mm"),column_title = paste0(df$Freq[i]," genes"), column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = T)
    #     expression log2rpkm profiles, use default redBlue color
    he.scaled=Heatmap(mat_scaled, name = "scaled log2rpkm", cluster_column = FALSE, top_annotation = ha, top_annotation_height = unit(4, "mm"),column_title = "Scaled Expression", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = TRUE)
    #     expression log2rpkm profiles, use default redBlue color
    he=Heatmap(mat, name = "log2rpkm",col = hmcol, cluster_column = FALSE, top_annotation = ha, top_annotation_height = unit(4, "mm"),column_title = "Expression", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = TRUE)
    draw(he+he.scaled+hm)
}
dev.off()
## if differential MSFs presence/absence is associated with differential expression patterns, we expected to see higher expression when MSFs present in 1kb pomoter, which was not found

# another way to check: scatter plots of A vs D, color those showed different presence and absence of MSFs; if associated, red colored points should be close to x and y axis
makeTransparent = function(..., alpha=0.5) {
    if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
    alpha = floor(255*alpha)
    newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
    .makeTransparent = function(col, alpha) {
        rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
    }
    newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
    return(newColor)
}
pdf("expNuc/ogPromoterMSFs.AvsD.pdf")
plot(expr[,c("A2","D5")], pch=20, col =rgb(0,0,0,0.4))
points(expr[apply(paMSF[,c("A2","D5")],1,function(x)x[1]!=x[2]),c("A2","D5")], col=rgb(1,0,0,0.4),pch=20)
plot(expr[,c("F1.At","F1.Dt")], pch=20, col =rgb(0,0,0,0.4))
points(expr[apply(paMSF[,c("F1.At","F1.Dt")],1,function(x)x[1]!=x[2]),c("F1.At","F1.Dt")], col=rgb(1,0,0,0.4),pch=20)
plot(expr[,c("AD1.At","AD1.Dt")], pch=20, col =rgb(0,0,0,0.4))
points(expr[apply(paMSF[,c("AD1.At","AD1.Dt")],1,function(x)x[1]!=x[2]),c("AD1.At","AD1.Dt")], col=rgb(1,0,0,0.4),pch=20)
dev.off()

#################################################################################
## Visualize chromatin profiles related to duplicated gene expression patterns ##
#################################################################################
library(GenomicRanges)
#library(RColorBrewer)
#library(ggbio)
library(GenomicFeatures)
library(genomation)
options(scipen=999)

# load annotation results of DNS segment annotation
l=load("isegRes/SegAnnotation/dns/genomicAnnotation.byType.rdata")
# "sampleInfo"       "sl"               "sl.nonzero"       "Sums" "Neg_peakAnnoList" "Pos_peakAnnoList"
sampleInfo$txDB=gsub("../../","",as.character(sampleInfo$txDB))
sampleInfo$inputBG=gsub("../","isegRes/",as.character(sampleInfo$inputBG))
sampleInfo$fullBW = c("aggregation/A*n*full.bw","aggregation/combo/D*full.bw","aggregation/combo/F*full.bw","aggregation/combo/M*full.bw")
sampleInfo

## load expression datasets
# read count
load("RNAseq/Ranalysis/expression.Dref.rdata")->l
l # "all"  "info" "A2"   "D5"   "F1"   "AD1"
# expression patterns
load("RNAseq/Ranalysis/regPattern.rdata")->l
l # "res"           "dominance.AD1" "dominance.F1"
# ortholog group info
load("expNuc/ogPromoterMSFs.rdata")->l
l #  "info"   "og6"    "paMSF"  "paMSFn" "expr"
rownames(og6)=as.character(og6$D5)
res=res[as.character(og6$D5),]
res$diploid = ifelse(res$cisNtrans=="A=0","A2=D5",ifelse(res$A>0&res$A.padj<0.05,"A2>D5","A2<D5"))
res$F1 = ifelse(res$cis=="B=0","At=Dt",ifelse(res$B>0&res$B.padj<0.05,"At>Dt","At<Dt"))
res$AD1 = ifelse(is.na(res$Bp)|is.na(res$Bp.padj)|res$Bp.padj>=0.05,"At=Dt",ifelse(res$Bp>0,"At>Dt","At<Dt"))
res$diF1 = paste(res$diploid, res$F1)
dominance.AD1 = dominance.AD1[as.character(og6$D5),]
dominance.F1=dominance.F1[as.character(og6$D5),]

source("expNuc.source.r")

## expression bias in AD1, when Bp!=0
i=4
bws="aggregation/combo/McD_q20.full.bw"
txdb = loadDb(sampleInfo$txDB[i])
refgeneL=split(rownames(res), res$AD1)
lapply(refgeneL,length)
AgeneL=lapply(refgeneL,function(x)og6[x,"AD1.At"])
DgeneL=lapply(refgeneL,function(x)og6[x,"AD1.Dt"])
DtssL = lapply(DgeneL,function(genes)getTSS(txdb, include=as.character(genes)))
AtssL = lapply(AgeneL,function(genes)getTSS(txdb, include=as.character(genes)))
pdf("expNuc/plotDup.AD1bias.pdf")
for(w in names(refgeneL)){
    plotBWsWindows(bws,list(A=AtssL[[w]], D=DtssL[[w]]), compare="window", main=paste(w,length(AtssL[[w]])))
}
dev.off()

## expression bias in F1, when B!=0
i=3
bws="aggregation/combo/FcD_q20.full.bw"
txdb = loadDb(sampleInfo$txDB[i])
refgeneL=split(rownames(res), res$F1)
lapply(refgeneL,length)
AgeneL=lapply(refgeneL,function(x)og6[x,"F1.At"])
DgeneL=lapply(refgeneL,function(x)og6[x,"F1.Dt"])
DtssL = lapply(DgeneL,function(genes)getTSS(txdb, include=as.character(genes)))
AtssL = lapply(AgeneL,function(genes)getTSS(txdb, include=as.character(genes)))
pdf("expNuc/plotDup.F1bias.pdf")
for(w in names(refgeneL)){
    plotBWsWindows(bws,list(A=AtssL[[w]], D=DtssL[[w]]), compare="window", main=paste(w,length(AtssL[[w]])))
}
dev.off()
# try with qnormed F1 At vs Dt
bws=c("aggregation/parseF/FcD_q20.full.chrA_qnorm.bw","aggregation/parseF/FcD_q20.full.chrD_qnorm.bw")
names(bws)=c("F1.At_qnorm","F1.Dt_qnorm")
refgeneL=split(rownames(res), res$F1)
lapply(refgeneL,length)
i=2
txdb = loadDb(sampleInfo$txDB[i])
DgeneL=lapply(refgeneL,function(x)og6[x,"D5"])
DtssL = lapply(DgeneL,function(genes)getTSS(txdb, include=genes))
i=1
txdb = loadDb(sampleInfo$txDB[i])
AgeneL=lapply(refgeneL,function(x)og6[x,"A2"])
AtssL = lapply(AgeneL,function(genes)getTSS(txdb, include=genes))
pdf("expNuc/plotDup.F1bias_qnorm.pdf")
for(w in names(refgeneL)){
    plotBWsWindows(bws,list(A=AtssL[[w]], D=DtssL[[w]]), compare="paired", main=paste(w,length(AtssL[[w]])))
}
dev.off()

## parental expression divergence A2 vs D5, when A!=0
bws=c("aggregation/A6Dn_q20.full.bw","aggregation/combo/DcD_q20.full.bw")
names(bws)=c("A2","D5")
refgeneL=split(rownames(res), res$diploid)
lapply(refgeneL,length)
i=2
txdb = loadDb(sampleInfo$txDB[i])
DgeneL=lapply(refgeneL,function(x)og6[x,"D5"])
DtssL = lapply(DgeneL,function(genes)getTSS(txdb, include=genes))
i=1
txdb = loadDb(sampleInfo$txDB[i])
AgeneL=lapply(refgeneL,function(x)og6[x,"A2"])
AtssL = lapply(AgeneL,function(genes)getTSS(txdb, include=genes))
pdf("expNuc/plotDup.A2vsD5.pdf")
for(w in names(refgeneL)){
    plotBWsWindows(bws,list(A=AtssL[[w]], D=DtssL[[w]]), compare="paired", main=paste(w,length(AtssL[[w]])))
}
dev.off()
# try with qnormed A2 vs D5
bws=c("aggregation/parseF/AcD_q20.full_qnorm.bw","aggregation/parseF/DcD_q20.full_qnorm.bw")
names(bws)=c("A2_qnorm","D5_qnorm")
pdf("expNuc/plotDup.A2vsD5_qnorm.pdf")
for(w in names(refgeneL)){
    plotBWsWindows(bws,list(A=AtssL[[w]], D=DtssL[[w]]), compare="paired", main=paste(w,length(AtssL[[w]])))
}
dev.off()

### Hybridization effect by cis/trans, Hr, F1 dominance, inspect A2 vs F1-At and D5 vs F1-Dt
color=c("purple","brown")
bwsA=c("aggregation/parseF/AcD_q20.full_qnorm.bw","aggregation/parseF/FcD_q20.full.chrA_qnorm.bw")
names(bwsA)=c("A2","At")
bwsD=c("aggregation/parseF/DcD_q20.full_qnorm.bw","aggregation/parseF/FcD_q20.full.chrD_qnorm.bw")
names(bwsD)=c("D5","Dt")
## all OG; before qnorm both diploid>F1 probably due to normalization; now only A2>At
DgeneL=list(og=rownames(res))
AgeneL=lapply(DgeneL,function(x)og6[x,"A2"])
DtssL = lapply(DgeneL,function(genes)getTSS(txdb=loadDb(sampleInfo$txDB[2]), include=as.character(genes)))
AtssL = lapply(AgeneL,function(genes)getTSS(txdb=loadDb(sampleInfo$txDB[1]), include=as.character(genes)))
pdf("expNuc/plotDup.ratio.OG.pdf")
plotBWsWindows(bwsA, AtssL, compare="bw", line.col=color, dispersion.col=makeTransparent(color))
plotBWsWindows(bwsD, DtssL, compare="bw", line.col=color, dispersion.col=makeTransparent(color))
dev.off()
## Hr: Hr>0 indicated increased contribution of At vs Dt, may be attributed to increased sensitivity of AtvsAt, or decreased sensitivity of Dt vs D5; hard to tell
DgeneL=split(rownames(res), res$Hr.reg)
AgeneL=lapply(DgeneL,function(x)og6[x,"A2"])
DtssL = lapply(DgeneL,function(genes)getTSS(txdb=loadDb(sampleInfo$txDB[2]), include=as.character(genes)))
AtssL = lapply(AgeneL,function(genes)getTSS(txdb=loadDb(sampleInfo$txDB[1]), include=as.character(genes)))
pdf("expNuc/plotDup.ratio.Hr.pdf")
plotBWsWindows(bwsA, AtssL, compare="bw", line.col=color, dispersion.col=makeTransparent(color))
plotBWsWindows(bwsD, DtssL, compare="bw", line.col=color, dispersion.col=makeTransparent(color))
dev.off()
## Dominance: all categories mirror OG, thus not relevant
DgeneL=split(rownames(res), res$dominance.F1)
AgeneL=lapply(DgeneL,function(x)og6[x,"A2"])
DtssL = lapply(DgeneL,function(genes)getTSS(txdb=loadDb(sampleInfo$txDB[2]), include=as.character(genes)))
AtssL = lapply(AgeneL,function(genes)getTSS(txdb=loadDb(sampleInfo$txDB[1]), include=as.character(genes)))
pdf("expNuc/plotDup.ratio.Dominance.pdf")
plotBWsWindows(bwsA, AtssL, compare="bw", line.col=color, dispersion.col=makeTransparent(color))
plotBWsWindows(bwsD, DtssL, compare="bw", line.col=color, dispersion.col=makeTransparent(color))
dev.off()
## Cis and trans: if cis only, expect At resemble A2 and Dt resemble D5;
DgeneL=split(rownames(res), res$category)
AgeneL=lapply(DgeneL,function(x)og6[x,"A2"])
DtssL = lapply(DgeneL,function(genes)getTSS(txdb=loadDb(sampleInfo$txDB[2]), include=as.character(genes)))
AtssL = lapply(AgeneL,function(genes)getTSS(txdb=loadDb(sampleInfo$txDB[1]), include=as.character(genes)))
pdf("expNuc/plotDup.ratio.cistrans.pdf")
plotBWsWindows(bwsA, AtssL, compare="bw", line.col=color, dispersion.col=makeTransparent(color))
plotBWsWindows(bwsD, DtssL, compare="bw", line.col=color, dispersion.col=makeTransparent(color))
dev.off()
## A2vsD5 and AtvsDt: if cis only, expect At resemble A2 and Dt resemble D5;
DgeneL=split(rownames(res), res$diF1)
AgeneL=lapply(DgeneL,function(x)og6[x,"A2"])
DtssL = lapply(DgeneL,function(genes)getTSS(txdb=loadDb(sampleInfo$txDB[2]), include=as.character(genes)))
AtssL = lapply(AgeneL,function(genes)getTSS(txdb=loadDb(sampleInfo$txDB[1]), include=as.character(genes)))
pdf("expNuc/plotDup.ratio.diF1.pdf")
plotBWsWindows(bwsA, AtssL, compare="bw", line.col=color, dispersion.col=makeTransparent(color))
plotBWsWindows(bwsD, DtssL, compare="bw", line.col=color, dispersion.col=makeTransparent(color))
dev.off()

### Hybridization effect by cis/trans, Hr, F1 dominance, inspect log2 ratios of A2 vs F1-At and D5 vs F1-Dt
bws=c("aggregation/parseF/FAtvsA2.log2ratio.bw","aggregation/parseF/FAtvsA2.log2ratio.bw")
names(bws)=c("At/A2","Dt/D5")
pdf("expNuc/plotDup.FAA2vsFDD5.pdf")
## all OG; before qnorm both diploid>F1 probably due to normalization; now only A2>At
## Hr: Hr>0 indicated increased contribution of At vs Dt, may be attributed to increased sensitivity of AtvsAt, or decreased sensitivity of Dt vs D5; hard to tell
## Dominance: all categories mirror OG, thus not relevant
## A2vsD5 and AtvsDt: if cis only, expect At resemble A2 and Dt resemble D5;
refgeneL=c(list(og=rownames(res)), split(rownames(res), res$Hr.reg), split(rownames(res), res$dominance.F1), split(rownames(res), res$diF1), split(rownames(res), res$category))
DgeneL=lapply(refgeneL,function(x)og6[x,"D5"])
AgeneL=lapply(refgeneL,function(x)og6[x,"A2"])
DtssL = lapply(DgeneL,function(genes)getTSS(txdb= loadDb(sampleInfo$txDB[2]), include=genes))
AtssL = lapply(AgeneL,function(genes)getTSS(txdb= loadDb(sampleInfo$txDB[1]), include=genes))
for(w in names(refgeneL)){
    plotBWsWindows(bws,list(A=AtssL[[w]], D=DtssL[[w]]), compare="paired", main=paste(w,length(AtssL[[w]])), line.col=color, dispersion.col=makeTransparent(color))
}
dev.off()



