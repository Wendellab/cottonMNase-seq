# R
# Generating gene annotations for peak calls and evaluating overall distributions across dataset
# Functional enrichment of gene annotations obtained from peak calls

# Load libraries
library(GenomicFeatures)
library(genomation)

#####################
## Gene Annotation ##
#####################

## TM1 saski
txdb<- makeTxDbFromGFF("~/jfw-lab/GenomicResources/archived_resources/AD1Saski/annotation/Ghirsutum_458_v1.1.gene_exons.gff3", format="gff3")
txdb
## restrict features to chromomes, excluding scaffolds
seqlevels(txdb)
seqlevels(txdb)<-seqlevels(txdb.TM1.saski)[1:26]
# original levels are seqlevels0(txdb)
saveDb(txdb, file="refGenomes/txdb.TM1saski.sqlite")
txdb.TM1saski <- loadDb("refGenomes/txdb.TM1saski.sqlite")

## TM1 UTX 2020
txdb<- makeTxDbFromGFF("/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/Ghirsutum_527_v2.1.gene_exons.gff3", format="gff3")
txdb
## restrict features to chromomes, excluding scaffolds
seqlevels(txdb)
seqlevels(txdb)<-seqlevels(txdb)[1:26]
# original levels are seqlevels0(txdb)
saveDb(txdb, file="refGenomes/txdb.AD1utx.sqlite")
txdb.AD1utx <- loadDb("refGenomes/txdb.AD1utx.sqlite")

## A2 Du et al. 2018
txdb <- makeTxDbFromGFF("~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Du2018/1.PacBio-Gar-Assembly-v1.0/G.arboreum.Chr.v1.0.gff", format="gff3")
txdb
## restrict features to chromomes, excluding scaffolds
seqlevels(txdb)
seqlevels(txdb)<-seqlevels(txdb)[1:13]
# original levels are seqlevels0(txdb)
saveDb(txdb, file="refGenomes/txdb.A2du.sqlite")
txdb.A2 <- loadDb("refGenomes/txdb.A2du.sqlite")

## A2 WHU Huang et al. 2020
txdb <- makeTxDbFromGFF("/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/Garboreum_Shixiya1_WHUv3.0rc.gene.standard.gff3", format="gff3")
txdb
## restrict features to chromomes, excluding scaffolds
seqlevels(txdb)
seqlevels(txdb)<-seqlevels(txdb)[1:13]
# original levels are seqlevels0(txdb)
saveDb(txdb, file="refGenomes/txdb.A2WHU.sqlite")
txdb.A2WHU <- loadDb("refGenomes/txdb.A2WHU.sqlite")

## D5
txdb <- makeTxDbFromGFF("~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/Graimondii_221_v2.1.gene_exons.gff3", format="gff3")
txdb
## restrict features to chromomes, excluding scaffolds
seqlevels(txdb)
seqlevels(txdb)<-seqlevels(txdb)[1:13]
# original levels are seqlevels0(txdb)
saveDb(txdb, file="refGenomes/txdb.D5.sqlite")
txdb.D5 <- loadDb("refGenomes/txdb.D5.sqlite")

## F1 2020 = A2 WHU x D5; need to remove D5 "Name=.." to keep compartible
system("sed 's/^Chr/A/g' /work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/Garboreum_Shixiya1_WHUv3.0rc.gene.standard.gff3 >refGenomes/At.gff")
#|sed 's/[.]gene//g'
system("sed 's/^Chr/D/g' /work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/G.raimondii_JGI_221_v2.1.transcripts_exons.gff3 |sed 's/-JGI_221_v2.1//g' |sed 's/Name.*;Parent/Parent/g' | sed 's/Name.*//g'>refGenomes/Dt.gff")
system("cat refGenomes/At.gff refGenomes/Dt.gff >refGenomes/F2020.gff3")
txdb <- makeTxDbFromGFF("refGenomes/F2020.gff3", format="gff3")
txdb
## restrict features to chromomes, excluding scaffolds
seqlevels(txdb)
txdb1<-keepSeqLevels(txdb, grep("^(A|D)",seqlevels(txdb), value=T))
# original levels are seqlevels0(txdb)
saveDb(txdb, file="refGenomes/txdb.F2020.sqlite")
txdb.F2020 <- loadDb("refGenomes/txdb.F2020.sqlite")

## F1 = A2 x D5
system("sed 's/^Chr/A/g' ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Du2018/1.PacBio-Gar-Assembly-v1.0/G.arboreum.Chr.v1.0.rn.gff >refGenomes/At.gff")
system("sed 's/^Chr/D/g' ~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/Graimondii_221_v2.1.gene_exons.gff3 >refGenomes/Dt.gff")
system("cat refGenomes/At.gff refGenomes/Dt.gff >refGenomes/F1.gff ")
txdb <- makeTxDbFromGFF("refGenomes/F1.gff", format="gff3")
txdb
## restrict features to chromomes, excluding scaffolds
seqlevels(txdb)
seqlevels(txdb)<-seqlevels(txdb)[1:26]
# original levels are seqlevels0(txdb)
saveDb(txdb, file="refGenomes/txdb.F1.sqlite")
txdb.F1 <- loadDb("refGenomes/txdb.F1.sqlite")

## Check Features
promoters(txdb, upstream=2000, downstream=400)
genes(txdb)
transcripts(txdb)
exons(txdb)
cds(txdb)


################################
## Prepare genomic annotation ##
################################
library(GenomicRanges)
#library(RColorBrewer)
#library(ggbio)
library(GenomicFeatures)
library(genomation)
options(scipen=999)

## FUN
# prepare quantile BED files from GR
writeGRtoBED = function(gr, group="quantile",out="")
{
    df <- data.frame(seqnames=seqnames(gr), starts=start(gr)-1, ends=end(gr), names=c(rep(".", length(gr))),scores=c(rep(".", length(gr))),strands=strand(gr))
    if(!is.null(group)){
        id = gsub("[.]gene","",gr$gene_id)
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

txDB=c("refGenomes/txdb.A2WHU.sqlite","refGenomes/txdb.D5.sqlite", "refGenomes/txdb.F2020.sqlite","refGenomes/txdb.AD1utx.sqlite")
ref = data.frame(txDB, genome = gsub("refGenomes/txdb.|.sqlite","",txDB), subgenome =c(1,1,2,2))
# write gene annotation to bed
for(i in 1:4)
{
    # genomic features
    print(ref[i,])
    txdb = loadDb(as.character(ref$txDB[i]))
    genes.gr =  genes(txdb)
    print(genes.gr)
    genes.gr=genes.gr[grep("^Chr|^A|^D",seqnames(genes.gr)),]
    print(genes.gr)
    
    writeGRtoBED(genes.gr, group=NULL, out=paste0("refGenomes/", ref$genome[i], ".gene"))
    
    if(ref$subgenome[i]==2){
        writeGRtoBED(genes.gr[grep("^A",seqnames(genes.gr)),], group=NULL,out=paste0("refGenomes/", ref$genome[i],".gene.A"))
        writeGRtoBED(genes.gr[grep("^D",seqnames(genes.gr)),], group=NULL,out=paste0("refGenomes/", ref$genome[i],".gene.D"))
    }
}
list.files("refGenomes",pattern="gene.*bed")
# "A2WHU.gene.bed"    "AD1utx.gene.A.bed" "AD1utx.gene.bed" "AD1utx.gene.D.bed" "D5.gene.bed"       "F2020.gene.A.bed" "F2020.gene.bed"    "F2020.gene.D.bed"

# write quantile gene annotation to bed
l=load("RNAseq/Ranalysis/expression.Indiv.rdata");l
Exp_genome=c("A2",  "D5",  "F1",  "AD1")
for(i in 1:4)
{
    # genomic features
    print(ref[i,])
    txdb = loadDb(as.character(ref$txDB[i]))
    genes.gr =  genes(txdb)
    print(genes.gr)
    genes.gr=genes.gr[grep("^Chr|^A|^D",seqnames(genes.gr)),]
    print(genes.gr)
    length(genes.gr)
    expr=get(Exp_genome[i])$tpm
    
    writeGRtoBED(genes.gr, group="quantile",out=paste0("refGenomes/bedQ/", ref$genome[i]))
    
    if(ref$subgenome[i]==2){
        writeGRtoBED(genes.gr[grep("^A",seqnames(genes.gr)),], group="quantile",out=paste0("refGenomes/bedQ/", ref$genome[i],".A"))
        writeGRtoBED(genes.gr[grep("^D",seqnames(genes.gr)),], group="quantile",out=paste0("refGenomes/bedQ/", ref$genome[i],".D"))
    }
}
list.files("refGenomes/bedQ/",pattern="Q*bed")

# write young leaf based quantile gene annotation to bed
l=load("DNase/Ranalysis/expression.Indiv.rdata");l
Exp_genome=c("A2",  "D5",  "F1",  "AD1")
for(i in c(1,2,4))
{
    # genomic features
    print(ref[i,])
    txdb = loadDb(as.character(ref$txDB[i]))
    genes.gr =  genes(txdb)
    print(genes.gr)
    genes.gr=genes.gr[grep("^Chr|^A|^D",seqnames(genes.gr)),]
    print(genes.gr)
    length(genes.gr)
    expr=get(Exp_genome[i])$tpm
    
    writeGRtoBED(genes.gr, group="quantile",out=paste0("refGenomes/bedQ_YL/", ref$genome[i]))
    
    if(ref$subgenome[i]==2){
        writeGRtoBED(genes.gr[grep("^A",seqnames(genes.gr)),], group="quantile",out=paste0("refGenomes/bedQ_YL/", ref$genome[i],".A"))
        writeGRtoBED(genes.gr[grep("^D",seqnames(genes.gr)),], group="quantile",out=paste0("refGenomes/bedQ_YL/", ref$genome[i],".D"))
    }
}
list.files("refGenomes/bedQ/",pattern="Q*bed")

---book 060520
## Write gene bed files for ortho-homoelogs ##
# /lss/research/jfw-lab/Projects/MNase-seq/orthohomoeologQuadruplets101218.txt
# This list was made for A2du, D5, TM1saski
ogQ<-read.table("orthohomoeologQuadruplets101218.txt", sep="\t", header=TRUE)
head(ogQ)
# A
txdb=loadDb("refGenomes/txdb.A2du.sqlite")
genes.gr = genes(txdb)
genes.gr = genes.gr[genes.gr$gene_id %in% ogQ$A2,]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/A2du.og")
# D
txdb=loadDb("refGenomes/txdb.D5.sqlite")
genes.gr = genes(txdb)
genes.gr = genes.gr[genes.gr$gene_id %in% ogQ$D5,]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/D5.og")
# M
txdb=loadDb("refGenomes/txdb.TM1saski.sqlite")
genes.gr = genes(txdb)
genes.gr = genes.gr[genes.gr$gene_id %in% c(as.character(ogQ$At),as.character(ogQ$Dt)),]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/TM1saski.og")
# F
txdb=loadDb("refGenomes/txdb.F1.sqlite")
genes.gr = genes(txdb)
genes.gr = genes.gr[genes.gr$gene_id %in% c(as.character(ogQ$A2),as.character(ogQ$D5)),]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/F1.og")
# check
list.files("refGenomes",pattern="*og.bed")
# "A2du.og.bed"     "D5.og.bed"       "F1.og.bed"       "TM1saski.og.bed"

###################
## TE Annotation ##
###################
## FUN
# extract genomic regions from txDB
getPoint=function(gr,type=c("start","end"),be=1000, af=1000)
{
    grP =gr
    if(type=="start")
    {
        start(grP[strand(grP)=="+",]) = start(gr[strand(gr)=="+",])-be
        # end of the + strand genes must be equalized to start pos
        end(grP[strand(grP)=="+",])  = start(gr[strand(gr)=="+",])+af
        # start of the - strand genes must be equalized to end pos
        start(grP[strand(grP)=="-",]) = end(gr[strand(gr)=="-",])-af
        end(grP[strand(grP)=="-",]) = end(gr[strand(gr)=="-",])+be
    }
    if(type=="end")
    {
        # start of the + strand genes must be equalized to end pos
        start(grP[strand(grP)=="+",])  =end(gr[strand(gr)=="+",]) -be
        end(grP[strand(grP)=="+",]) = end(gr[strand(gr)=="+",]) + af
        # end of the - strand genes must be equalized to start pos
        start(grP[strand(grP)=="-",]) = start(gr[strand(gr)=="-",])-af
        end(grP[strand(grP)=="-",]) = start(gr[strand(gr)=="-",])+be
     }
    # remove duplicated
    grP=grP[!duplicated(grP),]
    seqlevels(grP)=unique(as.character(seqnames(grP)))
    return(grP)
}
# module load r-udunits2
library(ChIPseeker)
rns =c("seqnames","source","type","start","end","diversity","strand","sw_score","family")

####### TM1 saski
df=read.table("/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/TEannotation/Tx-JGI_G.hirsutum_v1.1.fa.out.gff",sep="\t", comment.char = "#",col.names=rns,stringsAsFactors =FALSE)
df$seqnames = gsub("AD1","",df$seqnames)
df$type2=gsub("[/].*","",df$type)
df$type3 =as.character(df$type)
df$type3[grep("LTR",df$type)]<-df$type[grep("LTR",df$type)]
df$type3[grep("DNA|MITE",df$type)]="TIR"
df$type3[grep("Helitron",df$type)]="Helitron"
gr=makeGRangesFromDataFrame(df,keep.extra.columns=TRUE)
gr.AD1= gr

####### D5
df=read.table("/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/TEannotation/G.raimondii_JGI_221_v2.0.assembly.fasta.out.gff",sep="\t", comment.char = "#",col.names=rns,stringsAsFactors =FALSE)
df$seqnames = gsub(".*_","",df$seqnames)
df$type2=gsub("[/].*","",df$type)
df$type3 =as.character(df$type)
df$type3[grep("LTR",df$type)]<-df$type[grep("LTR",df$type)]
df$type3[grep("DNA|MITE",df$type)]="TIR"
df$type3[grep("Helitron",df$type)]="Helitron"
gr=makeGRangesFromDataFrame(df,keep.extra.columns=TRUE)
gr.D5= gr

####### A2
df=read.table("/work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/TEannotation/G.raimondii_JGI_221_v2.0.assembly.fasta.out.gff",sep="\t", comment.char = "#",col.names=rns,stringsAsFactors =FALSE)
df$seqnames = gsub(".*_","",df$seqnames)
df$type2=gsub("[/].*","",df$type)
df$type3 =as.character(df$type)
df$type3[grep("LTR",df$type)]<-df$type[grep("LTR",df$type)]
df$type3[grep("DNA|MITE",df$type)]="TIR"
df$type3[grep("Helitron",df$type)]="Helitron"
gr=makeGRangesFromDataFrame(df,keep.extra.columns=TRUE)
gr.A2= gr

####### F1
gr.At=gr.A2
gr.Dt=gr.D5
seqlevels(gr.At) = gsub("Chr","A",seqlevels(gr.At))
seqlevels(gr.Dt) = gsub("Chr","D",seqlevels(gr.Dt))
gr.F1 = c(gr.At, gr.Dt)

save(gr.F1, gr.A2, gr.D5,gr.AD1,file="refGenomes/TEannotation.rdata")

