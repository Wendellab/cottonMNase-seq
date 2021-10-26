# module load r/3.5.0-py2-ufvuwmm
options(scipen=999)
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq//RanalysisAD1GNT/")
# setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase")
################################
## collect expression results ##
################################
readKallistoGNT=function(dir="", fileL=c(),gff3){
    require(rtracklayer)
    dt=readGFF(gff3)
    dt$target_id=paste0(dt$type,"::",dt$seqid,":",dt$start-1,"-",dt$end,"(",dt$strand,")")
    for ( file in fileL ) {
        print(file)
        nn<- gsub("/abundance.tsv","",file)
        x=read.table(paste0(dir,file), header=TRUE, sep="\t")
        ## check TE and gene transcripts
        x$class=ifelse(grepl("::",x$target_id),"TE","Gene")
        x$genome = ifelse(grepl("::",x$target_id),substring(gsub(".*::","",x$target_id),1,1),substring(gsub("Gohir[.]","",x$target_id),1,1))
        # a<-aggregate(x$est_counts, list(x$genome,x$class),sum)
        # a$perc=a$x/sum(a$x)
        # b<-aggregate(x$tpm, list(x$genome,x$class),sum)
        # b$perc=b$x/sum(b$x)
        # c<-merge(a,b,by=c("Group.1","Group.2"))
        # names(c)=c("genome","class","raw","rawP","tpm","tpmP")
        # print(c)
        ## aggregate TE expression by panTE
        te = x[x$class=="TE",]
        te$target_id=gsub("[.|?]","*",te$target_id) # replace (.) and (?) to (*)
        # table(te$target_id %in% dt$target_id) # should be all true
        te<-merge(te, dt[,c("target_id","Name","type","Classification")],by="target_id",all.x=T,all.y=T)
        te.count<-aggregate(te$est_counts,by=list(te$Name,te$genome,te$Classification),sum)
        names(te.count)=c("Name","genome","Classification", nn)
        te.tpm<-aggregate(te$tpm,by=list(te$Name,te$genome,te$Classification),sum)
        names(te.tpm)=c("Name","genome","Classification", nn)
        print("get TE res")
        ## gene
        gene  = x[x$class=="Gene",]
        temp.count <- gene[,c("target_id", "est_counts")]
        temp.tpm  <- gene[,c("target_id", "tpm")]
        names(temp.count)[2]<-nn
        names(temp.tpm )[2]<-nn
        print("get Gene res")
        # make big table
        if (!exists("countG")) { countG<-temp.count } else {countG<- merge(countG, temp.count, by="target_id")}
        if (!exists("tpm"))  { tpm <-temp.tpm } else {tpm <- merge(tpm, temp.tpm, by="target_id")}
        if (!exists("countT")) { countT<-te.count } else {countT <- merge(countT, te.count, by=c("Name","genome","Classification"),all.x=T,all.y=T)}
        if (!exists("tpmT"))  { tpmT <-te.tpm } else {tpmT <- merge(tpmT, te.tpm, by=c("Name","genome","Classification"),all.x=T,all.y=T)}
    }
    count=countG
    rownames(count)=gsub("[.]1$|^evm.model.","",count$target_id)
    count=count[,-1]
    if(length(grep("Gohir.1Z",rownames(count)))>0)
    {count=count[-grep("Gohir.1Z",rownames(count)),]}
    rownames(tpm)=gsub("[.]1$|^evm.model.","",tpm$target_id)
    tpm=tpm[,-1]
    if(length(grep("Gohir.1Z",rownames(tpm)))>0)
    {tpm=tpm[-grep("Gohir.1Z",rownames(tpm)),]}
    # group 1 to 4, 1 as the lowest expression, 5 the highest
    tpm$log2mean = apply(tpm,1,function(x)log2(mean(x)+1))
    breaks=quantile(tpm$log2mean[tpm$log2mean>0], probs=0:4/4)
    tpm$quantile=cut(tpm$log2mean, breaks, include.lowest=TRUE, labels=FALSE)
    tpm$quantile[is.na(tpm$quantile)]=0
    print(table(tpm$quantile))
    # group 0 as zero expression
    return(list(count=count,tpm=tpm,countTE=countT, tpmTE=tpmT))
}

## all against AD1utx ref
dir <- "/work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/MappingGNT/"
fileL=paste0(list.files(dir, patter="Maxxa.*[1-9]$"),"/abundance.tsv")
gff3="/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/AD1_UTX_v2.1.fa.mod.EDTA.TEanno.gff3"
M=readKallistoGNT(dir, fileL , gff3)
dir <- "/work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/mappingAD1gnt/"
fileL=paste0(list.files(dir, patter=".*[1-9]$"),"/abundance.tsv")
gff3="/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/AD1_UTX_v2.1.fa.mod.EDTA.TEanno.gff3"
ADF=readKallistoGNT(dir, fileL , gff3)
#  quick check up
lapply(M,dim)
lapply(ADF,dim)
table(rownames(M$count)==rownames(ADF$count))
table(rownames(M$tpm)==rownames(ADF$tpm))
table(M$countTE$Name==ADF$countTE$Name)
table(M$tpmTE$Name==ADF$tpmTE$Name)

## prep count table
## ----------------------
# gene
g<-cbind(M$count[,1:3], ADF$count[,c(1:7,9)])
g$class="Gene"
g$classM="primaryTranscript"
g$genome =substring(gsub("Gohir[.]","",rownames(g)),1,1)
# te
t<-merge(M$countTE, ADF$countTE[,-c(11)], by=c("Name","genome","Classification"))
t[is.na(t)]=0
dim(tf<-t[rowSums(t[,-c(1:3)])>0,]) # filter expressed 12191
t<-tf[,-c(1:3)]
rownames(t)=paste0(tf$Name,"_",tf$genome)
t$class="TE"
t$classM=gsub("DNA|MITE","TIR",tf$Classification)
t$genome=tf$genome
# combind
dim(count <- rbind(t,g))  # 87093    14

## examine TE TPM results
## ----------------------
# gene
g<-cbind(M$tpm[,1:3], ADF$tpm[,c(1:7,9)])
g$class="Gene"
g$classM="primaryTranscript"
g$genome =substring(gsub("Gohir[.]","",rownames(g)),1,1)
# te
t<-merge(M$tpmTE, ADF$tpmTE[,-c(11)], by=c("Name","genome","Classification"))
t[is.na(t)]=0
dim(tf<-t[rowSums(t[,-c(1:3)])>0,]) # filter expressed 12191
colSums(tf[tf$genome=="A",-c(1:3)])/colSums(tf[tf$genome=="D",-c(1:3)])
# SD5-Maxxa-S1 SD5-Maxxa-S2 SD5-Maxxa-S3    SD5-A2-S1    SD5-A2-S4    SD5-A2-S5
#   1.0547113    1.1100774    1.0760543    5.1603866    4.9879750    5.2863238
# SD5-A2xD5-S1 SD5-A2xD5-S2 SD5-A2xD5-S3    SD5-D5-S1    SD5-D5-S3
#     0.8851864    0.9508544    1.0092672    0.2020558    0.1860912
# ? more  At/Dt TE expression from F1 to AD1
t<-tf[,-c(1:3)]
rownames(t)=paste0(tf$Name,"_",tf$genome)
t$class="TE"
t$classM=gsub("DNA|MITE","TIR",tf$Classification)
t$genome=tf$genome
tg <-rbind(t,g)
nn<-rowsum(tg[,1:11],group=paste(tg$genome,tg$class));nn
#       SD5-Maxxa-S1 SD5-Maxxa-S2 SD5-Maxxa-S3 SD5-A2-S1  SD5-A2-S4  SD5-A2-S5
# A Gene    468544.10    471896.54    471929.44 917471.36 921879.097 926092.670
# A TE       21225.01     24438.80     28063.77  32955.93  30730.696  26658.747
# D Gene    490106.81    481649.21    473926.55  43186.43  41229.201  42205.590
# D TE       20124.00     22015.41     26080.26   6386.33   6160.956   5042.965
#       SD5-A2xD5-S1 SD5-A2xD5-S2 SD5-A2xD5-S3  SD5-D5-S1  SD5-D5-S3
# A Gene    500799.39     492510.8    484987.85  50120.038  56124.196
# A TE       16656.60      24298.6     31903.18   8864.564   5672.615
# D Gene    463726.95     457636.2    451498.68 897143.640 907720.205
# D TE       18817.06      25554.5     31610.24  43871.859  30482.987
nn/colSums(nn)
#      SD5-Maxxa-S1 SD5-Maxxa-S2 SD5-Maxxa-S3   SD5-A2-S1   SD5-A2-S4
#A Gene   0.46854414   0.47189656   0.47192947 0.917471405 0.921879123
#A TE     0.02122501   0.02443881   0.02806377 0.032955933 0.030730696
#D Gene   0.49010679   0.48164921   0.47392655 0.043186431 0.041229200
#D TE     0.02012400   0.02201541   0.02608026 0.006386331 0.006160957
#        SD5-A2-S5 SD5-A2xD5-S1 SD5-A2xD5-S2 SD5-A2xD5-S3   SD5-D5-S1
#A Gene 0.926092576   0.50079938   0.49251076   0.48498785 0.050120035
#A TE   0.026658747   0.01665660   0.02429860   0.03190318 0.008864565
#D Gene 0.042205594   0.46372697   0.45763619   0.45149870 0.897143665
#D TE   0.005042965   0.01881706   0.02555449   0.03161024 0.043871859
#        SD5-D5-S3
#A Gene 0.056124194
#A TE   0.005672615
#D Gene 0.907720113
#D TE   0.030482987
rowsum(tg[,1:11],group=tg$class)/colSums(tg[,1:11])
#     SD5-Maxxa-S1 SD5-Maxxa-S2 SD5-Maxxa-S3  SD5-A2-S1  SD5-A2-S4  SD5-A2-S5
#Gene   0.95865099   0.95354571   0.94585604 0.96065779 0.96310835 0.96829826
#TE     0.04134901   0.04645421   0.05414403 0.03934226 0.03689165 0.03170172
#     SD5-A2xD5-S1 SD5-A2xD5-S2 SD5-A2xD5-S3  SD5-D5-S1 SD5-D5-S3
#Gene   0.96452639    0.9501469   0.93648656 0.94726366 0.9638443
#TE     0.03547366    0.0498531   0.06351342 0.05273643 0.0361556
dim(tpm <- tg) # 87093    14

## Make some plots for TPM of TEs
## -----------------------------
library(reshape2)
# reshape
tfl<-melt(tf,id.vars=names(tf)[1:3])
tfl$classL=gsub("/.*","",tfl$Classification)
tfl$classS=gsub(".*/","",tfl$Classification)
tfl$class=gsub("DNA|MITE","TIR",tfl$Classification)
tfl$species = gsub("SD5-|-S.*","",tfl$variable)
save(tfl,file="test.rdata")
pdf("examineTEexpression.tpm.pdf")
# barplot
ggplot(tfl, aes(variable, y=value,fill=class)) +
  geom_bar(stat = "identity") +
  labs(x ="Library", y = "TPM") +
  scale_fill_brewer(palette="Paired") +
  theme_classic()+
  theme(strip.text = element_text(size = 6), axis.title.y = element_text(size=16, face="bold"),
      axis.title.x = element_text(size=14, face="bold"),
      axis.text.y = element_text(size=8),
      axis.text.x = element_text(size=8, angle=35, vjust=1, hjust=0.95),
      legend.text = element_text(size=10),
      legend.title=element_blank())
# barplot, by subgenome
ggplot(tfl, aes(variable, y=value,fill=class)) +
  geom_bar(stat = "identity") +
  facet_wrap(~genome, ncol = 5)+
labs(x ="Library", y = "TPM") +
scale_fill_brewer(palette="Paired") +
theme_classic()+
theme(strip.text = element_text(size = 6), axis.title.y = element_text(size=16, face="bold"),
      axis.title.x = element_text(size=14, face="bold"),
      axis.text.y = element_text(size=8),
      axis.text.x = element_text(size=8, angle=35, vjust=1, hjust=0.95),
      legend.text = element_text(size=10),
      legend.title=element_blank())
# grouped boxplot
ggplot(tfl, aes(x=class, y=log2(value+1), fill=class)) +
    labs(x ="Library", y = "log2(TPM+1)") +
    scale_fill_brewer(palette="Paired")+
   # geom_violin() + geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    geom_boxplot()+
    theme_classic()+
    theme(strip.text = element_text(size = 6), axis.title.y = element_text(size=16, face="bold"),
      axis.title.x = element_text(size=14, face="bold"),
      axis.text.y = element_text(size=8),
      axis.text.x = element_text(size=8, angle=35, vjust=1, hjust=0.95),
      legend.text = element_text(size=10),
      legend.title=element_blank())
# grouped boxplot
ggplot(tfl, aes(x=class, y=log2(value+1), fill=class)) +
    labs(x ="Library", y = "log2(TPM+1)") +
    scale_fill_brewer(palette="Paired")+
    facet_grid(genome~species) +
   # geom_violin() + geom_boxplot(width=0.1, color="grey", alpha=0.2) +
    geom_boxplot()+
    theme_classic()+
    theme(strip.text = element_text(size = 6), axis.title.y = element_text(size=16, face="bold"),
      axis.title.x = element_text(size=14, face="bold"),
      axis.text.y = element_text(size=8),
      axis.text.x = element_text(size=8, angle=35, vjust=1, hjust=0.95),
      legend.text = element_text(size=10),
      legend.title=element_blank())
dev.off()



# my custom function
library(ggplot2)
library(scales)
library(ape)
plotGrouping <- function(norm_log, color, shape, text, tip, save = "plotGrouping.pdf"){
    # norm<-sweep(total,2,info$lib_size,"/")*10^6
    # norm_log <- log2(norm+1)
    pca=prcomp(t(norm_log))
    dat = as.data.frame(pca$x)
    proportion<-summary(pca)$importance[2,1:2]*100
    proportion<-paste0(names(proportion)," (", proportion, "%)")
    p<-ggplot(aes(PC1, PC2, color=color,    shape=shape),data=dat) + geom_point() +xlab(proportion[1]) + ylab(proportion[2])
    pdf(save)
    print( p + geom_text(aes_string(x = "PC1", y = "PC2", label = text), color="grey", hjust = 0, nudge_x = 0.09) )
    
    hc<-hclust( dist(t(norm_log)) )
    tre <- as.phylo(hc)
    tre$tip.label <- as.character(tip)
    # try to match ggplot color: library(scale);
    # show_col(col4<-hue_pal()(4))
    tipCol <- color
    levels(tipCol) <- hue_pal()(nlevels(tipCol))
    plot(tre,tip.col =as.character(tipCol),  type="unrooted",cex=0.6, no.margin=TRUE)
    plot(tre,tip.col =as.character(tipCol), type="fan",cex=0.6, no.margin=TRUE)
    dev.off()   }

# plots
info = data.frame(sample=names(count)[1:11], genome = gsub("SD5-|-S.","",names(count[,1:11])), rep = c(rep(1:3,3),1,2))
plotGrouping(log2(tpm[,1:11]+1), color=info$genome, shape=info$genome, tip=info$sample, text=info$genome, save = "plotGrouping.AD1GNT.log2tpm.pdf")

save(all, count, tpm,file="expression.AD1GNTref.rdata")

# individal A2 ref
fileL=paste0(list.files(dir, patter="-A2-.*[1-9]$"),"/abundance.tsv")
gff3="/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/A2_WHU_v1.fa.mod.EDTA.TEanno.gff3"
A=readKallistoGNT(dir, fileL , gff3)

# individal D5 ref
fileL=paste0(list.files(dir, patter="-D5-.*[1-9]$"),"/abundance.tsv")
gff3="/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/D5_JGI_v2.fa.mod.EDTA.TEanno.gff3"
D=readKallistoGNT(dir, fileL , gff3)

# indivial F ref
fileL=paste0(list.files(dir, patter="A2xD5.*[1-9]$"),"/abundance.tsv")
gff3="/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/AD1_UTX_v2.1.fa.mod.EDTA.TEanno.gff3"
F=readKallistoGNT(dir, fileL , gff3)

#----book


####################################
## Compare AD1GNT GENE expression ##
####################################
load("expression.AD1GNTref.rdata")
# ignore TE
count<-count[count$class=="Gene",c(1:11,14)]
rownames(count)=gsub("[.].*","",gsub("Gohir.","Gohir_",rownames(count)))
tpm<-tpm[tpm$class=="Gene",c(1:11,14)]
rownames(tpm)=gsub("[.].*","",gsub("Gohir.","Gohir_",rownames(tpm)))

# OG
ogQ<-read.table("../orthohomoeolog052421.txt", sep="\t", header=FALSE,stringsAsFactor = FALSE)
og<-data.frame(gorai=c(paste0(ogQ$V5,".A"),paste0(ogQ$V5,".D")),gohir=c(ogQ$V3,ogQ$V4))
og$gohir<-gsub("[.].*","",gsub("Gohir.","Gohir_",og$gohir))
og<-og[order(og$gorai),]
head(og)
# prep
count<-cbind(count[og$gohir,],og)
rownames(count)=count$gorai
tpm<-cbind(tpm[og$gohir,],og)
rownames(tpm)=tpm$gorai

# make genome list
AD1 <-list(count = count[,1:3], tpm  =tpm[,1:3])
A2  <-list(count = count[,4:6], tpm  =tpm[,4:6])
F1  <-list(count = count[,7:9], tpm  =tpm[,7:9])
D5  <-list(count = count[,10:11], tpm  =tpm[,10:11])

getTotal=function(x){y=x[grep(".A$",rownames(x)),]+x[grep(".D$",rownames(x)),];rownames(y)=gsub(".A$","",rownames(y));return(y)}
A2$count = getTotal(A2$count)
A2$tpm = getTotal(A2$tpm)
D5$count = getTotal(D5$count)
D5$tpm = getTotal(D5$tpm)

# check A vs D
rbind(
c("Diploids colSums A2/D5",colSums(A2$tpm[,1:3])/c(colSums(D5$tpm[,1:2]),NA)),
c("AD1 colSums A/D",colSums(AD1$tpm[grep("A",rownames(AD1$tpm)),1:3])/colSums(AD1$tpm[grep("D",rownames(AD1$tpm)),1:3])),
c("F1 colSums A/D",colSums(F1$tpm[grep("A",rownames(F1$tpm)),1:3])/colSums(F1$tpm[grep("D",rownames(F1$tpm)),1:3]))
)
# F1 A/D >1 biasedly assign reads to A genome, but AD1 < or = 1; this is different from D5 ref results


#################
## DE analysis ##
#################
library(DESeq2)
library(gplots)

# count sig gene number
getSig<-function(res,fc.threshold=0,direction=NULL){
    sig<- res[res$padj<0.05 & !is.na(res$padj) & abs(res$log2FoldChange)>=fc.threshold,]
    if(is.null(direction)){
        n<-nrow(sig)
    }else if(direction=="up"){
        n<-nrow(sig[sig$log2FoldChange>0,])
    }else if(direction=="down"){
        n<-nrow(sig[sig$log2FoldChange<0,])
    }
    return(n)
}
# volcano plot
plotVolcano<-function(res, title)
{
    plot(res[,"log2FoldChange"], -log2(res[,"padj"]), main=title, xlab="log2FoldChange", ylab="-log2padj",pch=".",ylim=c(0,200))
    abline(h=-log2(0.05))
}

# A2 vs D5
count = cbind(A2$count,D5$count)
info = data.frame(sample=names(count), genome = gsub("SD5-|-S.","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res <- results(DESeq(dds), contrast=c("genome","A2","D5"))
print( summary(res,alpha=.05) ) # Higher: A 3837 D 4065; RSEM A 3320 D 3848
write.table(res, file="DE.A2vsD5.txt",  sep="\t")
# Get parental expression divergence, test for A=0
A=res

# F1: At vs Dt
count = cbind(F1$count[grep("A$",rownames(F1$count)),], F1$count[grep("D$",rownames(F1$count)),])
rownames(count) =gsub(".A$","",rownames(count))
names(count)=paste(names(count),c("A","A","A","D","D","D"),sep=".")
info = data.frame(sample=names(count), genome = c("A","A","A","D","D","D"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
colSums(count) # Kallisto assigned roughly equal reads to Dt and At, skipp follows
# libTotal = rep(colSums(F1$count),2)
# libSize = libTotal/mean(libTotal)
# sizeFactors(dds) = libSize
res = results(DESeq(dds),contrast=c("genome","A","D"))
print( summary(res,alpha=.05) ) # Higher: A 2423 D 2773; RSEM A 2084, D 2575
write.table(res, file="DE.F1.AtvsF1.Dt.txt",  sep="\t")
# Get allelic expression divergence in F1, test for B=0
B=res

# AD1: At vs Dt
count = cbind(AD1$count[grep("A$",rownames(AD1$count)),], AD1$count[grep("D$",rownames(AD1$count)),])
rownames(count) =gsub(".A$","",rownames(count))
names(count)=paste(names(count),c("A","A","A","D","D","D"),sep=".")
info = data.frame(sample=names(count), genome = c("A","A","A","D","D","D"))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
colSums(count)
#libTotal = rep(colSums(AD1$count),2)
#libSize = libTotal/mean(libTotal)
#sizeFactors(dds) = libSize
res = results(DESeq(dds),contrast=c("genome","A","D"))
print( summary(res,alpha=.05) ) # Higher: A 3048 D 3270; RSEM A 2678, D 3209
write.table(res, file="DE.AD1.AtvsAD1.Dt.txt",  sep="\t")
# homoeolog expression divergence in polyploid, test for Bp=0
Bp=res

# F1 vs AD1
count = cbind(F1$count,AD1$count)
info = data.frame(sample=names(count), genome = gsub("SD5-|-S.","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","Maxxa","A2xD5"))
print( summary(res,alpha=.05) ) # higher: AD1 5037, F1 7377; RSEM AD1 4047, F1 6389
write.table(res, file="DE.AD1vsF1.txt",  sep="\t")
# differential expression by genome doubling
W =res

# Total comparison
############################
getTotal=function(x){y=x[grep(".A$",rownames(x)),]+x[grep(".D$",rownames(x)),];rownames(y)=gsub(".A$","",rownames(y));return(y)}
F1.t = getTotal(F1$count)
AD1.t = getTotal(AD1$count)
libTotal =c(colSums(A2$count[,1:3]),colSums(D5$count[,1:2]))
libSize = libTotal/mean(libTotal)
mid = as.data.frame(cbind(A2$count[,1]/libSize[1]+D5$count[,1]/libSize[4], A2$count[,1]/libSize[1]+D5$count[,2]/libSize[5], A2$count[,2]/libSize[2]+D5$count[,1]/libSize[4], A2$count[,2]/libSize[2]+D5$count[,2]/libSize[5], A2$count[,3]/libSize[3]+D5$count[,1]/libSize[4], A2$count[,3]/libSize[3]+D5$count[,2]/libSize[5]))/2
names(mid)=paste0("SD5-mid-",c("S11","S12","S21","S22","S31","S32"))
total = cbind(A2$count,D5$count,F1.t,AD1.t,mid)
colSums(total)

# check total grouping
library(ggplot2)
library(scales)
library(ape)
plotGrouping <- function(norm_log, color, shape, text, tip, save = "plotGrouping.pdf"){
    # norm<-sweep(total,2,info$lib_size,"/")*10^6
    # norm_log <- log2(norm+1)
    pca=prcomp(t(norm_log))
    dat = as.data.frame(pca$x)
    proportion<-summary(pca)$importance[2,1:2]*100
    proportion<-paste0(names(proportion)," (", proportion, "%)")
    p<-ggplot(aes(PC1, PC2, color=color,    shape=shape),data=dat) + geom_point() +xlab(proportion[1]) + ylab(proportion[2])
    pdf(save)
    print( p + geom_text(aes_string(x = "PC1", y = "PC2", label = text), color="grey", hjust = 0, nudge_x = 0.09) )
    
    hc<-hclust( dist(t(norm_log)) )
    tre <- as.phylo(hc)
    tre$tip.label <- as.character(tip)
    # try to match ggplot color: library(scale);
    # show_col(col4<-hue_pal()(4))
    tipCol <- color
    levels(tipCol) <- hue_pal()(nlevels(tipCol))
    plot(tre,tip.col =as.character(tipCol),  type="unrooted",cex=0.6, no.margin=TRUE)
    plot(tre,tip.col =as.character(tipCol), type="fan",cex=0.6, no.margin=TRUE)
    dev.off()   }
# plots
info = data.frame(sample=names(total), genome = gsub("SD5-|-S.*","",names(total)))
libTotal=colSums(total)
libSize = libTotal/mean(libTotal)
rpm = sweep(total,2,libSize,FUN="/")
plotGrouping(log2(rpm+1), color=info$genome, shape=info$genome, tip=info$sample, text=1:17, save = "plotGrouping.total.log2rpm.pdf")


# F1.total vs Mid
count = total[,grep("A2xD5|mid",names(total))]
info = data.frame(sample=names(count), genome = gsub("SD5-|-S.*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","A2xD5","mid"))
print( summary(res,alpha=.05) ) # higher: F1 5532 Mid 4658; RSEM F1 5071, Mid 4260
write.table(res, file="DE.F1vsMid.txt",  sep="\t")
# differential expression by genome doubling
F1vsMid =res

# AD1.total vs Mid
count = total[,grep("Maxxa|mid",names(total))]
info = data.frame(sample=names(count), genome = gsub("SD5-|-S.*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","Maxxa","mid"))
print( summary(res,alpha=.05) ) # higher: AD1 5622, Mid 5795; RSEM AD1 4822, Mid 4864
write.table(res, file="DE.AD1vsMid.txt",  sep="\t")
# differential expression by genome doubling
AD1vsMid =res

# F1.total vs A2
count = total[,grep("A2xD5|-A2-",names(total))]
info = data.frame(sample=names(count), genome = gsub("SD5-|-S.*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","A2xD5","A2"))
print( summary(res,alpha=.05) ) # higher: F1 4440, A2 3050; RSEM F1 4153, A2 2783
write.table(res, file="DE.F1vsA2.txt",  sep="\t")
# differential expression by genome doubling
F1vsA2 =res

# AD1.total vs A2
count = total[,grep("Maxxa|-A2-",names(total))]
info = data.frame(sample=names(count), genome = gsub("SD5-|-S.*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","Maxxa","A2"))
print( summary(res,alpha=.05) ) # higher AD1 3903, A2 3437; RSEM AD1 3297, A2 2740
write.table(res, file="DE.AD1vsA2.txt",  sep="\t")
# differential expression by genome doubling
AD1vsA2 =res


# F1.total vs D5
count = total[,grep("A2xD5|-D5-",names(total))]
info = data.frame(sample=names(count), genome = gsub("SD5-|-S.*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","A2xD5","D5"))
print( summary(res,alpha=.05) ) # higher RSEM F1 2491, D5 1208
write.table(res, file="DE.F1vsD5.txt",  sep="\t")
# differential expression by genome doubling
F1vsD5 =res

# AD1.total vs D5
count = total[,grep("Maxxa|-D5-",names(total))]
info = data.frame(sample=names(count), genome = gsub("SD5-|-S.*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","Maxxa","D5"))
print( summary(res,alpha=.05) ) # higher RSEM AD1 2087, D5 2080
write.table(res, file="DE.AD1vsD5.txt",  sep="\t")
# differential expression by genome doubling
AD1vsD5 =res

# Mid.total vs A2
count = total[,grep("mid|-A2-",names(total))]
info = data.frame(sample=names(count), genome = gsub("SD5-|-S.*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","mid","A2"))
print( summary(res,alpha=.05) ) # higher RSEM Mid 4438, A2 2383
write.table(res, file="DE.MidvsA2.txt",  sep="\t")
# differential expression by genome doubling
MidvsA2 =res


# Mid.total vs D5
count = total[,grep("mid|-D5-",names(total))]
info = data.frame(sample=names(count), genome = gsub("SD5-|-S.*","",names(count)))
dds <- DESeqDataSetFromMatrix( countData = round(count,0), colData =info, design = ~ genome)
res = results(DESeq(dds),contrast=c("genome","mid","D5"))
print( summary(res,alpha=.05) ) # higher RSEM Mid 2755, D5 1546
write.table(res, file="DE.MidvsD5.txt",  sep="\t")
# differential expression by genome doubling
MidvsD5 =res

# save
save(list=c("A","B","Bp","W", "total", ls(pattern="vs")), file="diffExprPattern.rdata")

# summary
unique(rownames(A)==rownames(B)) # TRUE
unique(rownames(A)==rownames(Bp)) # TRUE
# make summary table
sigT<-c("sample 1","sample 2","Total","DE (q<0.05)","1>2","2>1")
for(file in list.files(pattern="^DE"))
{
    res<-read.table(file,sep="\t",header=TRUE)
    sigRes <- c(unlist(strsplit(gsub("DE.|.txt","",file),split="vs") ), nrow(res), getSig(res),getSig(res,direction="up"),getSig(res,direction="down"))
    sigT<-rbind(sigT,sigRes)
}
T<-as.data.frame(sigT[-1,], row.names=FALSE)
names(T)<-sigT[1,]
T
#pdf plot
pdf("checkDE.pdf")
# table DE
textplot(T)
mtext("DE analysis result summary")
# volcano plot
plotVolcano(A, "Diploid A2 vs D5")
plotVolcano(B, "Hybrid alleic At vs Dt")
plotVolcano(Bp, "Polyploid homoeologous At vs Dt")
plotVolcano(W, "Genome doubling F1 vs AD1")
for(i in ls(pattern="vs")){plotVolcano(get(i),i)}
# compare log2 Fold change
plot( A[,"log2FoldChange"], B[,"log2FoldChange"],pch=".",main="log2 Fold Change",xlab="A - diploids",ylab="B - F1")
lines(c(-6,6),c(-6,6),col="blue")
plot( A[,"log2FoldChange"], Bp[,"log2FoldChange"],pch=".",main="log2 Fold Change",xlab="A - diploids",ylab="Bp - Polyploids")
lines(c(-6,6),c(-6,6),col="blue")
plot( B[,"log2FoldChange"], Bp[,"log2FoldChange"],pch=".",main="log2 Fold Change",xlab="B - F1",ylab="Bp - Polyploids")
lines(c(-6,6),c(-6,6),col="blue")
plot( W[grep("A$",rownames(W)),"log2FoldChange"], W[grep("D$",rownames(W)),"log2FoldChange"],pch=".",main="log2 Fold Change",xlab="W.At",ylab="W.Dt")
lines(c(-6,6),c(-6,6),col="blue")
plot( W[grep("A$",rownames(W)),"log2FoldChange"], W[grep("D$",rownames(W)),"log2FoldChange"],pch=".",main="log2 Fold Change",xlab="W.At",ylab="W.Dt")
lines(c(-6,6),c(-6,6),col="blue")
plot( F1vsA2[,"log2FoldChange"], F1vsD5[,"log2FoldChange"],pch=".",main="log2 Fold Change")
lines(c(-6,6),c(-6,6),col="blue")
plot( AD1vsA2[,"log2FoldChange"], AD1vsD5[,"log2FoldChange"],pch=".",main="log2 Fold Change")
lines(c(-6,6),c(-6,6),col="blue")
plot( MidvsA2[,"log2FoldChange"], MidvsD5[,"log2FoldChange"],pch=".",main="log2 Fold Change")
lines(c(-6,6),c(-6,6),col="blue")

dev.off()


########################
## cis-trans analysis ##
########################
load("diffExprPattern.rdata")->l
l
### comparing A-B= 0 is tricky, both are log2FoldChange and its standard error lfcse
# maybe I can compare with t test
#### T test from means and standard errors ####
# m1, m2: the sample means
# s1, s2: the sample standard deviations
# n1, n2: the same sizes
# se1, se2: the sample standard errors
# se1 <- s1/sqrt(n)
# m0: the null value for the difference in means to be tested for. Default is 0.
# equal.variance: whether or not to assume equal variance. Default is FALSE.
t.test2 <- function(m1,m2,se1,se2,n1,n2,m0=0,equal.variance=FALSE)
{
    if( equal.variance==FALSE )
    {
        # se <- sqrt( (s1^2/n1) + (s2^2/n2) )
        se <- sqrt( (se1^2) + (se2^2) )
        # welch-satterthwaite df
        # df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
        df <- ( (se1^2 + se2^2)^2 )/( (se1^2)^2/(n1-1) + (se2^2)^2/(n2-1) )
    } else
    {
        # pooled standard deviation, scaled by the sample sizes
        # se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
        se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) )
        df <- n1+n2-2
    }
    t <- (m1-m2-m0)/se
    dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))
    names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
    return(dat)
}

x1 = rnorm(3)
x2 = rnorm(3)
# you'll find this output agrees with that of t.test when you input x1,x2
t.test(x1,x2)
t.test2( mean(x1),  mean(x2), sd(x1)/sqrt(3), sd(x2)/sqrt(3), 3,3)

# function to make a categorization table
criteria<- as.data.frame(rbind(c("A!=0;B!=0;A=B", "1.Cis only"),
c("A!=0;B=0;A!=B", "2.Trans only"),
c("A!=0;B!=0;A!=B", "Cis+Trans"),
c("A=0;B!=0;A!=B", "5.Compensatory"),
c("A=0;B=0;A=B", "6.Conserved") ))
names(criteria) <- c("class","category")
classCisTrans<-function(A.res, B.res, A.n, B.n, log2fc.threshold=0,plotTitle=NULL)
{
    # A = log2(TX2094/Maxxa), cis + trans
    A <- A.res[,c("log2FoldChange", "lfcSE", "padj")]
    names(A) <- c("A", "A.SE", "A.padj")
    # A = log2(F1_t/F1_m), cis
    B <- B.res[,c("log2FoldChange", "lfcSE", "padj")]
    names(B) <- c("B", "B.SE", "B.padj")
    
    A<-A[rownames(B),]
    table <- cbind(A,B)
    table$AminusB <- table$A - table$B
    table$AminusB.pvalue <- apply(table,1,function(x) t.test2(m1=x[1],m2=x[4],se1=x[2], se2=x[5], n1=A.n, n2=B.n)["p-value"])
    
    table$cisNtrans <- ifelse(abs(table$A)>=log2fc.threshold & table$A.padj<0.05 & !is.na(table$A.padj), "A!=0", "A=0")
    table$cis <- ifelse(abs(table$B)>=log2fc.threshold & table$B.padj<0.05 & !is.na(table$B.padj), "B!=0", "B=0")
    table$trans <- ifelse(abs(table$AminusB)>=log2fc.threshold & table$AminusB.pvalue<0.05, "A!=B", "A=B")
    table$class <- paste(table$cisNtrans,table$cis,table$trans,sep=";")
    table$category <- as.character(criteria$category[ match(table$class,criteria$class)])
    table$category[is.na(table$category)] <- "7.Ambiguous"
    table$category[ table$category=="Cis+Trans" & table$B*table$AminusB >0 ] <- "3.Cis+Trans: enhancing"
    table$category[ table$category=="Cis+Trans" & table$B*table$AminusB <0 ] <- "4.Cis+Trans: compensating"
    
    colors <- c("red","blue","purple","brown","green","black","grey")
    if(!is.null(plotTitle)){
        p<- ggplot( table, aes(x=A, y=B, color=category)) + geom_point(alpha=0.8) +  xlab("Cis + Trans") + ylab("Cis") + ggtitle(plotTitle) + scale_color_manual(values=colors) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +  geom_abline(intercept = 0, slope = 1) + theme_bw()
        # p2<-ggplot( table, aes( factor(category), fill=category)) + geom_bar(stat="count") + scale_fill_manual(values=colors)
        print(p)
    }
    return(table)
}

# function to plot cis trans results
plotCisTrans<-function(table, plotTitle="")
{
    colors <- c("red","blue","purple","brown","green","black","grey")
    p<- ggplot( table, aes(x=A, y=B, color=category)) + geom_point(alpha=0.8) +  xlab("Cis + Trans") + ylab("Cis") + ggtitle(plotTitle) + scale_color_manual(values=colors) + geom_vline(xintercept=0) + geom_hline(yintercept=0) +  geom_abline(intercept = 0, slope = 1) + theme_bw()
    # p2<-ggplot( table, aes( factor(category), fill=category)) + geom_bar(stat="count") + scale_fill_manual(values=colors)
    print(p)
}

# categorization using A and B
res <- classCisTrans(A.res = A, B.res = B, A.n = 3, B.n=3, log2fc.threshold=0)

# make plot
pdf("plotCistrans.pdf")
textplot(data.frame(table(res$category)),cex=0.6)
plotCisTrans(as.data.frame(res), "Cis/Trans")
dev.off()

# ouput
write.table(res, file ="regPattern.txt",sep="\t")
save(res, file="regPattern.rdata")

#############################
## Genome Evolution Impact ##
#############################
load("diffExprPattern.rdata")
load("regPattern.rdata")
# Categorization of additional effects according to extended analytic framework
# Hr = B - A
# Pr = Bp - A
# Wr = Bp - B
classEffects<-function(A.res, B.res, Bp.res, A.n, B.n, Bp.n, log2fc.threshold=0)
{
    # A = log2(A2/D5)
    A <- A.res[,c("log2FoldChange", "lfcSE", "padj")]
    names(A) <- c("A", "A.SE", "A.padj")
    # B = log2(F1_At/F1_Dt)
    B <- B.res[,c("log2FoldChange", "lfcSE", "padj")]
    names(B) <- c("B", "B.SE", "B.padj")
    # Bp = log2(AD1_At/AD1_Dt)
    Bp <- Bp.res[,c("log2FoldChange", "lfcSE", "padj")]
    names(Bp) <- c("Bp", "Bp.SE", "Bp.padj")
    
    # make sure all rownames identical
    table <- cbind(A,B,Bp)
    
    # Hr
    table$Hr <- table$B - table$A
    table$Hr.pvalue <- apply(table,1,function(x) t.test2(m1=x[4],m2=x[1],se1=x[5], se2=x[2], n1=B.n, n2=A.n)["p-value"])
    
    # Pr
    table$Pr <- table$Bp - table$A
    table$Pr.pvalue <- apply(table,1,function(x) t.test2(m1=x[7],m2=x[1],se1=x[8], se2=x[2], n1=Bp.n, n2=A.n)["p-value"])
    
    # Wr
    table$Wr <- table$Bp - table$B
    table$Wr.pvalue <- apply(table,1,function(x) t.test2(m1=x[7],m2=x[4],se1=x[8], se2=x[5], n1=Bp.n, n2=B.n)["p-value"])
    
    # interpretation
    table$Hr.reg <- ifelse(abs(table$Hr)>=log2fc.threshold & table$Hr.pvalue<0.05 & !is.na(table$Hr.pvalue), "Hr!=0", "Hr=0")
    table$Hr.reg[table$Hr>=log2fc.threshold & table$Hr.reg=="Hr!=0"] <- "Hr>0, stronger At"
    table$Hr.reg[table$Hr<=-log2fc.threshold & table$Hr.reg=="Hr!=0"] <- "Hr<0, stronger Dt"
    
    table$Pr.reg <- ifelse(abs(table$Pr)>=log2fc.threshold & table$Pr.pvalue<0.05 & !is.na(table$Pr.pvalue), "Pr!=0", "Pr=0")
    table$Pr.reg[table$Pr>=log2fc.threshold & table$Pr.reg=="Pr!=0"] <- "Pr>0, stronger At"
    table$Pr.reg[table$Pr<=-log2fc.threshold & table$Pr.reg=="Pr!=0"] <- "Pr<0, stronger Dt"
    
    table$Wr.reg <- ifelse(abs(table$Wr)>=log2fc.threshold & table$Wr.pvalue<0.05 & !is.na(table$Wr.pvalue), "Wr!=0", "Wr=0")
    table$Wr.reg[table$Wr>=log2fc.threshold & table$Wr.reg=="Wr!=0"] <- "Wr>0, stronger At"
    table$Wr.reg[table$Wr<=-log2fc.threshold & table$Wr.reg=="Wr!=0"] <- "Wr<0, stronger Dt"
    
    return(table)
}

# categorization using A, B and Bp
effect = classEffects(A, B, Bp, A.n=3, B.n=3, Bp.n=3, log2fc.threshold=0)
unique(effect[,1:6]==res[,1:6])
res=cbind(res, effect[,-c(1:6)])

# make plots
sumT =  rbind(c("Regulation Pattern","Measure","A","D"),
c("Diploid divergence, A2vsD5", "A", getSig(A,direction="up"), getSig(A,direction="down")),
c("Homoeolog bias, F1", "B", getSig(B,direction="up"), getSig(B,direction="down")),
c("Homoeolog bias, AD1", "Bp", getSig(Bp,direction="up"), getSig(Bp,direction="down")),
c("Hybridization effect direction", "Hr", length(grep("Hr>0",res$Hr.reg)), length(grep("Hr<0",res$Hr.reg))),
c("Allopolyploidy effect direction", "Pr", length(grep("Pr>0",res$Pr.reg)), length(grep("Pr<0",res$Pr.reg))),
c("Genome doubling effect direction", "Wr", length(grep("Wr>0",res$Wr.reg)), length(grep("Wr<0",res$Wr.reg)))
)
sumT
# "Regulation Pattern"               "Measure" "A"    "D"
# "Diploid divergence, F1"           "A"       "3320" "3848"
# "Homoeolog bias, F1"               "B"       "2084" "2575"
# "Homoeolog bias, AD1"              "Bp"      "2678" "3209"
# "Hybridization effect direction"   "Hr"      "504"  "294"
# "Allopolyploidy effect direction"  "Pr"      "1480" "1457"
# "Genome doubling effect direction" "Wr"      "1051" "1160"

pdf("plotEvoImpact.pdf")
textplot(data.frame(sumT),cex=0.6)
# compare impacts
plot( res[,"Hr"], res[,"Pr"],pch=".", main=paste("cor = ", cor(res[,"Hr"],res[,"Pr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
plot( res[,"Hr"], res[,"Wr"],pch=".", main=paste("cor = ", cor(res[,"Hr"],res[,"Wr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
plot( res[,"Wr"], res[,"Pr"],pch=".", main=paste("cor = ", cor(res[,"Wr"],res[,"Pr"],use="complete.obs")))
lines(c(-10,10),c(-10,10),col="blue")
dev.off()

# ouput
write.table(res, file ="regPattern.txt",sep="\t")
save(res, file="regPattern.rdata")

###################################
## Expression Dominance Analysis ##
###################################
load("diffExprPattern.rdata")
load("regPattern.rdata")
# TvsMid
# TvsA
# TvsD
classDominance<-function(TvsMid, TvsA, TvsD, log2fc.threshold=0)
{
    # Hybrid/polyploid vs Mid parental val
    TvsMid <- data.frame(TvsMid[,c("log2FoldChange", "lfcSE", "padj")])
    names(TvsMid) <- c("TvsMid", "TvsMid.SE", "TvsMid.padj")
    # Hybrid/polyploid vs Parent 1
    TvsA <- data.frame(TvsA[,c("log2FoldChange", "lfcSE", "padj")])
    names(TvsA) <- c("TvsA", "TvsA.SE", "TvsA.padj")
    # Hybrid/polyploid vs Parent 2
    TvsD <- data.frame(TvsD[,c("log2FoldChange", "lfcSE", "padj")])
    names(TvsD) <- c("TvsD", "TvsD.SE", "TvsD.padj")
    
    tbl = cbind(TvsMid, TvsA, TvsD)
    
    tbl$TvsMid[is.na(tbl$TvsMid)] =0
    tbl$TvsA[is.na(tbl$TvsA)] =0
    tbl$TvsD[is.na(tbl$TvsD)] =0
    
    # judge
    tbl$additivity <- ifelse(abs(tbl$TvsMid)>=log2fc.threshold & tbl$TvsMid.padj<0.05 & !is.na(tbl$TvsMid.padj), "T!=Mid", "T=Mid")
    tbl$TvsA.reg <- ifelse(abs(tbl$TvsA)>=log2fc.threshold & tbl$TvsA.padj<0.05 & !is.na(tbl$TvsA.padj), "T!=A", "T=A")
    tbl$TvsD.reg <- ifelse(abs(tbl$TvsD)>=log2fc.threshold & tbl$TvsD.padj<0.05 & !is.na(tbl$TvsD.padj), "T!=D", "T=D")
    
    # together
    tbl$class <- paste(tbl$additivity, tbl$TvsA.reg, tbl$TvsD.reg, sep=";")
    
    # assign category
    tbl$category = "Other non-additivity"
    tbl$category[grep("T=Mid",tbl$class)] = "Additivity"
    tbl$category[grep("T!=Mid;T=A;T!=D",tbl$class)] = "A-dominant"
    tbl$category[grep("T!=Mid;T!=A;T=D",tbl$class)] = "D-dominant"
    tbl$category[grepl("T!=Mid;T!=A;T!=D",tbl$class) & tbl$TvsA>0 & tbl$TvsD>0] = "Transgressive Up"
    tbl$category[grepl("T!=Mid;T!=A;T!=D",tbl$class) & tbl$TvsA<0 & tbl$TvsD<0] = "Transgressive Down"
    
    return(tbl)
}

# categorization using A, B and Bp
dominance.F1 = classDominance(F1vsMid, F1vsA2, F1vsD5, log2fc.threshold=0)
dominance.AD1 = classDominance(AD1vsMid, AD1vsA2, AD1vsD5, log2fc.threshold=0)

# make plots
sumT =  rbind(table(dominance.AD1$category), table(dominance.F1$category))
rownames(sumT)=c("AD1","F1")
pdf("plotDominance.pdf")
textplot(data.frame(sumT),cex=0.6)
dev.off()

# merge
res$dominance.AD1=dominance.AD1$category
res$dominance.F1 =dominance.F1$category

# ouput
write.table(res, file ="regPattern.txt",sep="\t")
save(res,dominance.AD1,dominance.F1, file="regPattern.rdata")

###################################
## Summarize Regulatory Patterns ##
###################################
load("regPattern.rdata")

sumTbl =as.data.frame(rbind(
c("Diploid divergence, higher expresson", "A", length(which(res$A>0 & res$A.padj<0.05 & !is.na(res$A.padj))), length(which(res$A<0 & res$A.padj<0.05 & !is.na(res$A.padj)))),
c("Homoeolog bias in F1, higher expression", "B", length(which(res$B>0 & res$B.padj<0.05 & !is.na(res$B.padj))), length(which(res$B<0 & res$B.padj<0.05 & !is.na(res$B.padj)))),
c("Homoeolog bias in AD1, higher expression", "Bp", length(which(res$Bp>0 & res$Bp.padj<0.05 & !is.na(res$Bp.padj))), length(which(res$Bp<0 & res$Bp.padj<0.05 & !is.na(res$Bp.padj)))),
c("F1 total expression vs diploids, DE" ,"F1vsA2 & F1vsD5", length(which(dominance.F1$TvsA.reg =="T!=A")), length(which(dominance.F1$TvsD.reg =="T!=D"))),
c("AD1 total expression vs diploids, DE" ,"AD1vsA2 & AD1vsD5", length(which(dominance.AD1$TvsA.reg =="T!=A")), length(which(dominance.AD1$TvsD.reg =="T!=D"))),
c("Hybridization effect, up regulate", "Hr", length(grep("Hr>0",res$Hr.reg)), length(grep("Hr<0",res$Hr.reg))),
c("Allopolyploidy effect, up regulate", "Pr", length(grep("Pr>0",res$Pr.reg)), length(grep("Pr<0",res$Pr.reg))),
c("Genome doubling effect, up regulate", "Wr", length(grep("Wr>0",res$Wr.reg)), length(grep("Wr<0",res$Wr.reg)))
))
names(sumTbl)=c("Regulation Pattern","Measure","A","D")
sumTbl[,c("A","D")]=apply(sumTbl[,c("A","D")],2,as.numeric)
sumTbl$balance=ifelse(sumTbl$A>sumTbl$D,"A>D","A<D")
sumTbl$chisqTest.pval = round(apply(sumTbl[,c("A","D")],1,function(x)chisq.test(x)$"p.value") ,6)
sumTbl
library(gplots)
pdf("plotRegSummary.pdf")
textplot(sumTbl,cex=0.6)
dev.off()

----

# restrict to Orthpolog groups
ogQ<-read.table("../../orthohomoeologQuadruplets101218.txt", sep="\t", header=TRUE)
head(ogQ)
select=as.character(ogQ$D5)
res=res[select,]
dominance.F1=dominance.F1[select,]
dominance.AD1=dominance.AD1[select,]
sumTbl =as.data.frame(rbind(
c("Diploid divergence, higher expresson", "A", length(which(res$A>0 & res$A.padj<0.05 & !is.na(res$A.padj))), length(which(res$A<0 & res$A.padj<0.05 & !is.na(res$A.padj)))),
c("Homoeolog bias in F1, higher expression", "B", length(which(res$B>0 & res$B.padj<0.05 & !is.na(res$B.padj))), length(which(res$B<0 & res$B.padj<0.05 & !is.na(res$B.padj)))),
c("Homoeolog bias in AD1, higher expression", "Bp", length(which(res$Bp>0 & res$Bp.padj<0.05 & !is.na(res$Bp.padj))), length(which(res$Bp<0 & res$Bp.padj<0.05 & !is.na(res$Bp.padj)))),
c("F1 total expression vs diploids, DE" ,"F1vsA2 & F1vsD5", length(which(dominance.F1$TvsA.reg =="T!=A")), length(which(dominance.F1$TvsD.reg =="T!=D"))),
c("AD1 total expression vs diploids, DE" ,"AD1vsA2 & AD1vsD5", length(which(dominance.AD1$TvsA.reg =="T!=A")), length(which(dominance.AD1$TvsD.reg =="T!=D"))),
c("Hybridization effect, up regulate", "Hr", length(grep("Hr>0",res$Hr.reg)), length(grep("Hr<0",res$Hr.reg))),
c("Allopolyploidy effect, up regulate", "Pr", length(grep("Pr>0",res$Pr.reg)), length(grep("Pr<0",res$Pr.reg))),
c("Genome doubling effect, up regulate", "Wr", length(grep("Wr>0",res$Wr.reg)), length(grep("Wr<0",res$Wr.reg)))
))
names(sumTbl)=c("Regulation Pattern","Measure","A","D")
sumTbl[,c("A","D")]=apply(sumTbl[,c("A","D")],2,as.numeric)
sumTbl$balance=ifelse(sumTbl$A>sumTbl$D,"A>D","A<D")
sumTbl$chisqTest.pval = round(apply(sumTbl[,c("A","D")],1,function(x)chisq.test(x)$"p.value") ,6)
sumTbl
pdf("plotRegSummary.og20362.pdf")
textplot(sumTbl,cex=0.6)
dev.off()

## draw complex heapmap check relationships between different patterns
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
load("expression.Dref.rdata")->l
l # "all"  "info" "A2"   "D5"   "F1"   "AD1"
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
for(cat in levels(factor(res$category))[1:5] ){
    print(cat)
    select = which(res$category ==cat)
    # tpm expression of A2, D5, F1.At, F1.Dt, AD1.At, AD2.Dt
    tpm = cbind(A2$tpm[,1:3],D5$tpm[,1:2],F1$tpm[grep("A$",rownames(F1$tpm)),1:3],F1$tpm[grep("D$",rownames(F1$tpm)),1:3],AD1$tpm[grep("A$",rownames(AD1$tpm)),1:3],AD1$tpm[grep("D$",rownames(AD1$tpm)),1:3])
    info = data.frame(sample=names(tpm), species = gsub("SD5-|-S.","",names(tpm)), ploidy=c(rep("diploid",5), rep("F1",6),rep("AD1",6)), genome =c(rep("A",3),rep("D",2),rep(c("A","D"),each=3,2)) )
    info$id = paste(info$sample,info$origin,sep=".")
    names(tpm) =info$id
    # head annotation
    ha = HeatmapAnnotation(df = info[,c("ploidy","genome")], col = list(ploidy=structure(brewer.pal(3, "Set2"), names = c("diploid","F1","AD1")),genome=c("A"="pink","D"="royalblue")))
    # heatmap of log2tpm, not scaled
    hm.tpm=Heatmap(as.matrix(log2(tpm[select,]+1)), name = "log2tpm",col = hmcol, cluster_column = FALSE, top_annotation = ha,  top_annotation_height = unit(4, "mm"),column_title = "log2tpm", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = FALSE)
    # hm.tpm
    
    # log2(A/D) ratios of A, B, Bp
    log2ratio = data.frame(A = res$A, B = res$B, Bp = res$Bp)
    log2ratio$A[res$cisNtrans =="A=0"] = 0
    log2ratio$B[res$cis =="B=0"] = 0
    log2ratio$Bp[res$Bp.padj>0.05 | is.na(res$Bp.padj)] = 0
    # head annotation
    ha = HeatmapAnnotation(df = data.frame(ploidy=c("diploid","F1","AD1")), col = list(ploidy=structure(brewer.pal(3, "Set2"), names = c("diploid","F1","AD1"))) )
    # heatmap of log2ratio, use default redBlue color
    hm.ratio=Heatmap(as.matrix(log2ratio[select,]), name = "log2(A/D)", cluster_column = FALSE, top_annotation = ha,  top_annotation_height = unit(4, "mm"),column_title = "log2(A/D)", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = TRUE)
    # hm.ratio
    
    #    cis trans classes, use regCol
    ct = data.frame(category=res$category[select])
    hm.ct = rowAnnotation(df = ct, col = list(category = structure( c("red","blue","purple","brown","green","black","grey"), names=levels(factor(res$category)))), width = unit(1, "cm"))
    # hm.ct
    
    # additional categories of Hr, Pr, Wr
    impact = res[,c("Hr","Pr","Wr")]
    impact$Hr[res$Hr.reg == "Hr=0"] = 0
    impact$Pr[res$Pr.reg == "Pr=0"] = 0
    impact$Wr[res$Wr.reg == "Wr=0"] = 0
    # heatmap of log2ratio, use default redBlue color
    hm.impact=Heatmap(as.matrix(impact[select,]), name = "evolution impact", cluster_column = FALSE, column_title = "evolution impact", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = TRUE)
    # hm.impact
    
    pdf(paste0("plotCH.",cat,".pdf"))
    draw(hm.tpm+hm.ratio+hm.ct+hm.impact, newpage = TRUE, column_title = "Gene regulatory divergence", column_title_gp = gpar(fontsize = 12, fontface = "bold"), heatmap_legend_side = "bottom")
    dev.off()
}
### hard to make sense
