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
#-------
system("grep '^A' /work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/Ghirsutum_527_v2.1.gene_exons.gff3 |sed s/^A/Chr/g >temp.gff3")
txdb<- makeTxDbFromGFF("temp.gff3", format="gff3")
saveDb(txdb, file="refGenomes/txdb.AD1utxA.sqlite")
#-------
system("grep '^D' /work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/Ghirsutum_527_v2.1.gene_exons.gff3 |sed s/^D/Chr/g >temp.gff3")
txdb<- makeTxDbFromGFF("temp.gff3", format="gff3")
saveDb(txdb, file="refGenomes/txdb.AD1utxD.sqlite")

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

## Write gene bed files for ortho-homoelogs ##
# v101218, made for A2du, D5, TM1saski: /lss/research/jfw-lab/Projects/MNase-seq/orthohomoeologQuadruplets101218.txt
# v063020, made for A2WHU, D5NSF, AD1UTX: /work/LAS/jfw-lab/jconover/Joel_cytonuclear/Synteny_fastas/Cotton_Jing/Shit-a-logs_out.txt
# v020321, made for A2WHU, D5NSF, AD1UTX: /work/LAS/jfw-lab/jconover/Joel_cytonuclear/Synteny_fastas/Cotton_Jing/Shit-a-logs_singletons.txt
# v052421, made for A2WHU, D5NSF, AD1UTX: /work/LAS/jfw-lab/jconover/Joel_cytonuclear/Synteny_fastas/Cotton_Jing/pSONIC.singeltons.txt
system("cp /work/LAS/jfw-lab/jconover/Joel_cytonuclear/Synteny_fastas/Cotton_Jing/pSONIC.singletons.txt orthohomoeolog052421.txt")
ogQ<-read.table("orthohomoeolog052421.txt", sep="\t", header=FALSE)
dim(ogQ) # 22889
head(ogQ)
# V1          V2                      V3                      V4
# 1 OG0002592 Gar01G00020 Gohir.A01G001000.1.v2.1 Gohir.D01G000700.1.v2.1
# 2 OG0002593 Gar01G00030 Gohir.A01G001100.1.v2.1 Gohir.D01G000800.1.v2.1
# 3 OG0002594 Gar01G00040 Gohir.A01G001200.1.v2.1 Gohir.D01G000900.1.v2.1
names(ogQ)<-c("og","A2","At","Dt","D5")
# A2WHU
txdb=loadDb("refGenomes/txdb.A2WHU.sqlite")
genes.gr = genes(txdb)
table(genes.gr$gene_id %in% paste0(ogQ$A2,".gene")) # F 20389, T 22889
genes.gr = genes.gr[genes.gr$gene_id %in% paste0(ogQ$A2,".gene"),]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/A2WHU.og")
# D
txdb=loadDb("refGenomes/txdb.D5.sqlite")
genes.gr = genes(txdb)
table(genes.gr$gene_id %in% gsub("[.]1$","",ogQ$D5)) # F 14616 T 22889
genes.gr = genes.gr[genes.gr$gene_id %in% gsub("[.]1$","",ogQ$D5),]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/D5.og")
# AD1utcx
txdb=loadDb("refGenomes/txdb.AD1utx.sqlite")
genes.gr = genes(txdb)
table(genes.gr$gene_id %in% substr(ogQ$At,1,16)) # F 52487 T 22889
table(genes.gr$gene_id %in% substr(ogQ$Dt,1,16)) # F 52487 T 22889
table(genes.gr$gene_id %in% c(substr(ogQ$At,1,16),substr(ogQ$Dt,1,16))) # F 29598 T 45778
genes.gr = genes.gr[genes.gr$gene_id %in% c(substr(ogQ$At,1,16),substr(ogQ$Dt,1,16)),]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/AD1utx.og")
# F12020
txdb=loadDb("refGenomes/txdb.F2020.sqlite")
genes.gr = genes(txdb)
table(genes.gr$gene_id %in% c(paste0(ogQ$A2,".gene"),gsub("[.]1$","",ogQ$D5))) # F35005 T 45778
genes.gr = genes.gr[genes.gr$gene_id %in% c(paste0(ogQ$A2,".gene"),gsub("[.]1$","",ogQ$D5)),]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/F2020.og")
# check
list.files("refGenomes",pattern="*og.bed")
# "A2WHU.og.bed"  "AD1utx.og.bed" "D5.og.bed"     "F2020.og.bed"

#####################
## Func Annotation ##
#####################





###################
## TE Annotation ##
###################
## FUN
# extract genomic regions from txDB
library(rtracklayer)
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


# AD1
zz=gzfile('refGenomes/TEannotation/v010721/AD1_UTX_v2.1.fa.mod.EDTA.TEanno.gff3.gz','rt')
df=read.table(zz,sep="\t", comment.char = "#",col.names=rns,stringsAsFactors =FALSE)
df$seqnames = gsub("AD1UTX_","",df$seqnames)
# Method=structural indicates intact TE
table(df$type,grepl("structural",df$family))
#                                FALSE   TRUE
#  CACTA_TIR_transposon         116242    537
#  Copia_LTR_retrotransposon    177190   4999
#  Gypsy_LTR_retrotransposon    962485   9514
#  hAT_TIR_transposon           100369   2481
#  helitron                       8893    137
#  long_terminal_repeat              0  35394
#  LTR_retrotransposon          980116   3184
#  Mutator_TIR_transposon       449918   2030
#  PIF_Harbinger_TIR_transposon  56514    219
#  repeat_region                     0  17697
#  target_site_duplication           0  35362
#  Tc1_Mariner_TIR_transposon    17487     94
table(df$strand)
#      -       ?       .       +
#1475551   19086    5361 1480864
df$strand[df$strand%in%c("?",".")]="*"
gr=makeGRangesFromDataFrame(df,keep.extra.columns=TRUE)
grG=gr[grep("scf|tig|ctg|scaffold",seqnames(gr),invert=TRUE),]
grG=grG[grG$type=="Gypsy_LTR_retrotransposon",]
export(grG,'refGenomes/AD1utx.gypsy.bed')
export(grG[grep("structural",grG$family),],'refGenomes/AD1utx.gypsy.intact.bed')
gr.AD1= gr

# A2
zz=gzfile('refGenomes/TEannotation/v010721/A2_WHU_v1.fa.mod.EDTA.TEanno.gff3.gz','rt')
df=read.table(zz,sep="\t", comment.char = "#",col.names=rns,stringsAsFactors =FALSE)
df$seqnames = gsub("A2WHU_","",df$seqnames)
table(df$type,grepl("structural",df$family))
#                               FALSE   TRUE
# CACTA_TIR_transposon          57915    274
# Copia_LTR_retrotransposon     93471   1445
# Gypsy_LTR_retrotransposon    778308  14193
# hAT_TIR_transposon            49277   1156
# helitron                       4218     69
# long_terminal_repeat              0  35040
# LTR_retrotransposon          804683   1882
# Mutator_TIR_transposon       203211    867
# PIF_Harbinger_TIR_transposon  27086    111
# repeat_region                     0  17520
# target_site_duplication           0  35006
# Tc1_Mariner_TIR_transposon     8901     48
table(df$strand)
#      -       ?       .       +
#1060994   11466    2456 1059765
df$strand[df$strand%in%c("?",".")]="*"
gr=makeGRangesFromDataFrame(df,keep.extra.columns=TRUE)
grG=gr[grep("scf|tig|ctg|scaffold",seqnames(gr),invert=TRUE),]
grG=grG[grG$type=="Gypsy_LTR_retrotransposon",]
export(grG,'refGenomes/A2WHU.gypsy.bed')
export(grG[grep("structural",grG$family),],'refGenomes/A2WHU.gypsy.intact.bed')
gr.A2= gr

# D5
zz=gzfile('refGenomes/TEannotation/v010721/D5_JGI_v2.fa.mod.EDTA.TEanno.gff3.gz','rt')
df=read.table(zz,sep="\t", comment.char = "#",col.names=rns,stringsAsFactors =FALSE)
df$seqnames = gsub("D5JGI_","",df$seqnames)
table(df$type,grepl("structural",df$family))
#                             FALSE   TRUE
# CACTA_TIR_transposon          59970    359
# Copia_LTR_retrotransposon     84493    992
# Gypsy_LTR_retrotransposon    232925    592
# hAT_TIR_transposon            51138   1630
# helitron                       4708     90
# long_terminal_repeat              0   3612
# LTR_retrotransposon          212563    222
# Mutator_TIR_transposon       222762   1080
# PIF_Harbinger_TIR_transposon  29283    113
# repeat_region                     0   1806
# target_site_duplication           0   3600
# Tc1_Mariner_TIR_transposon     8887     64
table(df$strand)
#      -       ?       .       +
#    458544   1266   3246 457833
df$strand[df$strand%in%c("?",".")]="*"
gr=makeGRangesFromDataFrame(df,keep.extra.columns=TRUE)
grG=gr[grep("scf|tig|ctg|scaffold",seqnames(gr),invert=TRUE),]
grG=grG[grG$type=="Gypsy_LTR_retrotransposon",]
export(grG,'refGenomes/D5.gypsy.bed')
export(grG[grep("structural",grG$family),],'refGenomes/D5.gypsy.intact.bed')
gr.D5= gr

####### F1
gr.At=gr.A2
gr.Dt=gr.D5
gr.F1 = c(gr.At, gr.Dt)
system("cat refGenomes/A2WHU.gypsy.bed refGenomes/D5.gypsy.bed > refGenomes/F12020.gypsy.bed")
system("cat refGenomes/A2WHU.gypsy.intact.bed refGenomes/D5.gypsy.intact.bed > refGenomes/F12020.gypsy.intact.bed")

save(gr.F1, gr.A2, gr.D5,gr.AD1,file="refGenomes/TEannotation.rdata")


############################################################
## Below for previous versions of TE annotation, archived ##
############################################################

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

