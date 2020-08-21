# Analysis of three hierarchical layers of chromosome interactions in tetraploid cotton
# within same chromosomes - AA.sc, DD.sc  vs diploid A2.sc D5.sc
# between chromosomes of the same subgenome - AA.dc, DD.dc vs diploid A2.dc and D5.dc
# between chromosomes of different subgenomes - AD = DA
library(ggplot2)
chrs =c("01","02","03","04","05","06","07","08","09","10","11","12","13")
fl=data.frame(file=c(grep("norm.res10MB.balance.txt",list.files("HiC/A2",full=T),value=T),grep("norm.res10MB.balance.txt",list.files("HiC/D5",full=T),value=T),grep("norm.res10MB.balance.txt",list.files("HiC/AD1",full=T),value=T)),genome=rep(c("A2","D5","AD1"),each=2),rep=rep(c(1:2),3))
pdf("HiC/compareCCI.pdf")
for(i in 1:nrow(fl)){
    f=as.character(fl$file[i])
    x=read.table(f, header=T,sep="\t")
    # clean matrix
    y=as.matrix(x[,-c(1:2)])
    # chr
    indx = gsub("-.*","",x[,1])
    
    if(fl$genome[i]=="AD1"){
        # inter-subgenome
        AD = as.numeric(as.matrix( y[grep("A",indx), grep("D",indx)] ))
        DA = as.numeric(as.matrix( y[grep("D",indx), grep("A",indx)] ))
        # intra-subgenome
        mat = as.matrix( y[grep("A",indx), grep("A",indx)] )
        AA = as.numeric(mat[upper.tri(mat,diag=T)])
        mat = as.matrix( y[grep("D",indx), grep("D",indx)] )
        DD = as.numeric(mat[upper.tri(mat,diag=T)])
        
        # within same chromosomes
        AA.sc=c()
        DD.sc=c()
        AA.dc=c()
        DD.dc=c()
        for(chr in chrs){
            mat = as.matrix( y[indx==paste0("A",chr), indx==paste0("A",chr)] )
            As = as.numeric(mat[upper.tri(mat,diag=T)])
            mat = as.matrix( y[indx==paste0("D",chr), indx==paste0("D",chr)] )
            Ds = as.numeric(mat[upper.tri(mat,diag=T)])
            Ad = as.numeric(as.matrix( y[indx==paste0("A",chr), indx!=paste0("A",chr)&grepl("A",indx)] ))
            Dd = as.numeric(as.matrix( y[indx==paste0("D",chr), indx!=paste0("D",chr)&grepl("D",indx)] ))
            AA.sc=c(AA.sc,As)
            DD.sc=c(DD.sc,Ds)
            AA.dc=c(AA.dc,Ad)
            DD.dc=c(DD.dc,Dd)
        }
        cciL = list (AA.sc=AA.sc,DD.sc=DD.sc,AA.dc=AA.dc,DD.dc=DD.dc,AA=AA,DD=DD,AD=AD)
        cciD=data.frame( cci=c(AA.sc,DD.sc,AA.dc,DD.dc,AA,DD,AD), genome=c(rep("A",length(AA.sc)),rep("D",length(DD.sc)),rep("A",length(AA.dc)),rep("D",length(DD.dc)),rep("A",length(AA)),rep("D",length(DD)),rep("A",length(AD))),  type=c(rep("withChr",length(AA.sc)+length(DD.sc)),rep("betweenChr",length(AA.dc)+length(DD.dc)),rep("withinSubgenome",length(AA)+length(DD)),rep("betweenSubgenome",length(AD))),file=f)
        #ggplot(cciD,aes(x=genome, y=cci,fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Dark2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6)) +labs(title=f)  + facet_wrap(~ type,nrow=1)+ ylab("Chromatin conformatin interaction")+scale_y_continuous(trans='log2')

    }else{
        # within same chromosomes
        sc=c()
        dc=c()
        for(chr in chrs){
            mat = as.matrix( y[indx==paste0("Chr",chr), indx==paste0("Chr",chr)] )
            sc1 = as.numeric(mat[upper.tri(mat,diag=T)])
            dc1 = as.numeric(as.matrix( y[indx==paste0("Chr",chr), indx!=paste0("Chr",chr)] ))
            sc=c(sc,sc1)
            dc=c(dc,dc1)
        }
        cciL = list (sc=sc,dc=dc)
        g=substring(fl$genome[i],1,1)
        cciD=data.frame(cci=c(sc,dc), genome=g, type=c(rep("withChr",length(sc)),rep("betweenChr",length(dc))),file=f)
    }
   boxplot(cciL,main = f)
   boxplot(lapply(cciL,log2),main = f)
   assign(paste0(fl$genome[i],"Rep",fl$rep[i]),cciL)
   assign(paste0("df",fl$genome[i],"Rep",fl$rep[i]),cciD)
}
dev.off()

df1 = rbind(dfA2Rep1,dfA2Rep2,dfD5Rep1,dfD5Rep2); df1$ploidy="Diploid"
df2 = rbind(dfAD1Rep1,dfAD1Rep2); df2$ploidy="AD1"
df=rbind(df1, df2)
df$type=factor(df$type,levels=c("withChr","betweenChr","withinSubgenome","betweenSubgenome"))
df$ploidy=factor(df$ploidy,levels=c("Diploid","AD1"))
p=ggplot(df,aes(x=ploidy, y=log2(cci),fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Dark2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6)) +labs(title="3D interaction")  + facet_wrap(~ type,nrow=1)+ ylab("Observed/expected ratio, log2")

library(gridExtra)
pdf("HiC/compareCCI.1.pdf")
grid.arrange(p, nrow = 2)
dev.off()

#####################
## A/B compartment ##
#####################
library(GenomicFeatures)
options(stringsAsFactors=F)

setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/HiC/")
fl=data.frame(file=grep(".pca.PC1.txt",c(list.files("A2",full=T),list.files("D5",full=T),list.files("AD1",full=T)),value=T),genome=rep(c("A2","D5","AD1"),each=2),rep=rep(c(1:2),3),txDB = rep(c("../refGenomes/txdb.A2WHU.sqlite","../refGenomes/txdb.D5.sqlite","../refGenomes/txdb.AD1utx.sqlite"),each=2))
# load young leaf expression
l=load("../DNase/Ranalysis/expression.Indiv.rdata");l # A2 D5, AD1
# l=load("../RNAseq/Ranalysis/expression.Indiv.rdata")
res=c("genome","rep","subgenome","B%","A%","NAgene%","Agene%","Bgene%")
pdf("compareAB.pdf")
for(i in 1:nrow(fl)){
    f=as.character(fl$file[i])
    print(fl[i,])
    x=read.table(f, header=F,sep="\t")
    names(x)=c("idx","seqnames","start","end","strand","PC1")
    c=makeGRangesFromDataFrame(x,keep.extra.columns=T)
    start(c)=start(c)+1
    #c.cov <- GenomicRanges::coverage(c, weight="PC1")
    
    ## check gene expression in A vs B compartments
    # get genesstart(c)
    txdb = loadDb(as.character(fl$txDB[i]))
    genes.gr =  genes(txdb)
    genes.gr=genes.gr[grep("^Chr|^A|^D",seqnames(genes.gr)),]
    genes.gr$gene_id =gsub("[.]gene","",genes.gr$gene_id)
    # print(genes.gr)
    # assign expression
    expr=get(as.character(fl$genome[i]))$tpm
    genes.gr$tpm = expr[genes.gr$gene_id,"log2mean"]
    # assign A/B
    #genes.gr$PC1 = apply(data.frame(genes.gr),1,function(x){mean(c.cov[[x[1]]][x[2]:x[3]])})
    hits = GenomicRanges::findOverlaps(genes.gr,c,ignore.strand=T)
    pc1 = aggregate(c$PC1[subjectHits(hits)], list(queryHits(hits)),mean)
    genes.gr$PC1 = NA
    genes.gr$PC1[pc1$Group.1]=pc1$x
    genes.gr$compartment=ifelse(is.na(genes.gr$PC1),"-",ifelse(genes.gr$PC1>0,"A","B"))
    # compare
    a1=aov(genes.gr$tpm~genes.gr$compartment)
    print(TukeyHSD(a1))
    boxplot(genes.gr$tpm~genes.gr$compartment,main=fl$file[i])
    
    
    ## proportion of A compartment
    ss=c(fl$genome[i],fl$rep[i],"",table(x$PC1>0)/nrow(x), table(genes.gr$compartment)/length(genes.gr))
    res=rbind(res,ss)
    if(fl$genome[i]=="AD1"){
        select=grep("A",x$seqnames)
        gr= genes.gr[grep("A",seqnames(genes.gr)),]
        ss=c(fl$genome[i],fl$rep[i],"A",table(x$PC1[select]>0)/nrow(x[select,]), table(gr$compartment)/length(gr))
        res=rbind(res,ss)
        select=grep("D",x$seqnames)
        gr= genes.gr[grep("D",seqnames(genes.gr)),]
        ss=c(fl$genome[i],fl$rep[i],"D",table(x$PC1[select]>0)/nrow(x[select,]), table(gr$compartment)/length(gr))
        res=rbind(res,ss)
        
    }
}
dev.off()
# post hoc tests show that genes in A are more highly expressed than genes in B compartment
res
colnames(res)=res[1,]
res=data.frame(res[-1,])
write.table(res,file="compareAB.txt",sep="\t",row.names=FALSE)
