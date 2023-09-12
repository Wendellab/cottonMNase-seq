# Motif Discovery and Enrichment

Guanjing Hu

10-20-2021

## Objectives

**Questions**

* How did the promoters differentiate between cotton (sub)genomes?
* What motifs are differentially enriched between cotton (sub)genomes?

## Input Data

### 1. Promoter accessible chromatin regions (ACRs) identified 

* `A2_ACRs_promoter1k.fasta`
* `D5_ACRs_promoter1k.fasta`
* `F1a_ACRs_promoter1k.fasta`
* `F1d_ACRs_promoter1k.fasta`
* `AD1a_ACRs_promoter1k.fasta`
* `AD1d_ACRs_promoter1k.fasta`

```R
df[df$Feature=="Promoter (<=1kb)",]
#            Feature  Number RegionSize NumberPerc RegionSizePerc Sample  ploidy
# 1  Promoter (<=1kb) 1046436      17894  0.3782327      0.3597869     A2 Diploid
# 12 Promoter (<=1kb)  849325      17826  0.3649802      0.3628998     D5 Diploid
# 23 Promoter (<=1kb)  863736      15950  0.3410552      0.3335355    F1a      F1
# 34 Promoter (<=1kb)  907249      17345  0.3687499      0.3653656    F1d      F1
# 45 Promoter (<=1kb) 1830576      24044  0.4619651      0.4185423   AD1a     AD1
# 56 Promoter (<=1kb) 1529891      22074  0.4615450      0.4205853   AD1d     AD1
```

### 2. Promter sequences as control

* `A2_promoter1k.fasta` - A2 WHU ref
* `D5_promoter1k.fasta` - D5 JGI ref
* `AD1a_promoter1k.fasta` - AD1 TM-1 UTX ref
* `AD1d_promoter1k.fasta` - AD1 TM-1 UTX ref

### 3. Motif database

* `Ath_TF_binding_motifs.meme` downloaded from http://planttfdb.gao-lab.org/download/motif/Ath_TF_binding_motifs.meme.gz , 619 motifs; `Ath_TF_binding_motifs_information.txt` downloaded from http://planttfdb.gao-lab.org/download/motif/Ath_TF_binding_motifs_information.txt

```R
## load motifi info
setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/Qregulation3")
x<-read.table("0.input/Ath_TF_binding_motifs_information.txt",header=T,sep="\t")
head(x)
#     Gene_id      Family Matrix_id              Species Method    Datasource Datasource_ID
# 1 AT1G01060 MYB_related   MP00119 Arabidopsis thaliana    DAP PMID:27203113     AT1G01060
# 2 AT1G01250         ERF   MP00120 Arabidopsis thaliana    DAP PMID:27203113     AT1G01250
# 3 AT1G01260        bHLH   MP00100 Arabidopsis thaliana    PBM PMID:26531826      MA0958.1
table(x$Family)
#         AP2         ARF       ARR-B          B3     BBR-BPC 
#          8           6           5          10           3 
#       BES1        bHLH        bZIP        C2H2         C3H 
#          6          32          35          27           4 
#      CAMTA         CPP         Dof      E2F/DP         EIL 
#          4           5          22           5           2 
#        ERF        FAR1     G2-like        GATA        GeBP 
#         80           2          21          15           3 
#       GRAS         GRF      HD-ZIP         HSF         LBD 
#          1           1          22          10           8 
#        LFY M-type_MADS   MIKC_MADS         MYB MYB_related 
#          1           2          18          67          23 
#        NAC       NF-YB    Nin-like         RAV   S1Fa-like 
#         59           2           3           4           1 
#        SBP         SRS         TCP    Trihelix         VOZ 
#         12           2          18          14           1 
#        WOX        WRKY       YABBY       ZF-HD 
#          3          43           2           7 
```

* `JASPAR2018_CORE_plants_non-redundant.meme` downloaded from https://meme-suite.org/meme/db/motifs

* `JASPAR2022_CORE_plants_non-redundant_pfms_transfac.txt` and `JASPAR2022_CORE_non-redundant_pfms_meme.txt` downloaded from https://jaspar.genereg.net/downloads/, contains 656 motifs with TF family information

  ```R
  ## load motifi info
   y<-as.data.frame(matrix(x,ncol=3,byrow=T))
  names(y)<-c("id","family","class")
  head(y)
  #         id family                             class
  # 1 MA0020.1    DOF Other C4 zinc finger-type factors
  # 2 MA0021.1    DOF Other C4 zinc finger-type factors
  # 3 MA0034.1    Myb        Tryptophan cluster factors
  # 4 MA0053.1    DOF Other C4 zinc finger-type factors
  motifs<-y
  ```
  
  

**MEME Installation**

Download `meme-5.4.1.tar.gz` from  https://meme-suite.org/meme/doc/download.html and follow the [online instruction](https://meme-suite.org/meme/doc/install.html?man_type=web#quick) 

```bash
cd /Users/guanjinghu/tools/meme-5.4.1
./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
make
make test
make install
export PATH=/Users/guanjinghu/meme/bin:/Users/guanjinghu/meme/libexec/meme-5.4.1:$PATH
```



## Analysis

### 1. Compare promoters

FIMO scanning motifs in 1kb promoters

```bash
## scan JASPAR 2018 nr
# promter
fimo --o ~/Desktop/jA2.p1k 0.input/JASPAR2018_CORE_plants_non-redundant.meme 0.input/A2_promoter1k.fasta
fimo --o ~/Desktop/jD5.p1k 0.input/JASPAR2018_CORE_plants_non-redundant.meme 0.input/D5_promoter1k.fasta
fimo --o ~/Desktop/jAD1a.p1k 0.input/JASPAR2018_CORE_plants_non-redundant.meme 0.input/AD1a_promoter1k.fasta
fimo --o ~/Desktop/jAD1d.p1k 0.input/JASPAR2018_CORE_plants_non-redundant.meme 0.input/AD1d_promoter1k.fasta

## scan JASPAR 2022 nr
# promter
fimo --o ~/Desktop/jA2.p1k 0.input/JASPAR2022_CORE_non-redundant_pfms_meme.txt 0.input/A2_promoter1k.fasta
fimo --o ~/Desktop/jD5.p1k 0.input/JASPAR2022_CORE_non-redundant_pfms_meme.txt 0.input/D5_promoter1k.fasta
fimo --o ~/Desktop/jAD1a.p1k 0.input/JASPAR2022_CORE_non-redundant_pfms_meme.txt 0.input/AD1a_promoter1k.fasta
fimo --o ~/Desktop/jAD1d.p1k 0.input/JASPAR2022_CORE_non-redundant_pfms_meme.txt 0.input/AD1d_promoter1k.fasta
# failed, error in meme
#The PSPM of motif MA0207.1 has probabilities which don't sum to 1 on row 6.
#The PSPM of motif MA0207.1 has unparsable characters "URL #http://jaspar.genereg.net/matrix/MA0207.1" on the end of row 6. Only spaces or tabs are allowed on the same line as the numbers.
#FATAL: An error occurred reading the motif file.

## scan At TFDB
# promoter
fimo --o ~/Desktop/jA2.p1k 0.input/Ath_TF_binding_motifs.meme 0.input/A2_promoter1k.fasta
fimo --o ~/Desktop/jD5.p1k 0.input/Ath_TF_binding_motifs.meme 0.input/D5_promoter1k.fasta
fimo --o ~/Desktop/jAD1a.p1k 0.input/Ath_TF_binding_motifs.meme 0.input/AD1a_promoter1k.fasta
fimo --o ~/Desktop/jAD1d.p1k 0.input/Ath_TF_binding_motifs.meme 0.input/AD1d_promoter1k.fasta
```

Compare motif occurance between genomes

```R
## load scanning results
setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/Qregulation3/1.PromoterComparison")
# loop input
fl<-paste0(list.files(path="~/Desktop",pattern="p1k",full.names=T),"/fimo.tsv"); fl
# "jA2.p1k/fimo.tsv"   "jAD1a.p1k/fimo.tsv" "jAD1d.p1k/fimo.tsv" "jD5.p1k/fimo.tsv" 
for(i in fl)
  {
  x<-fread(i,select=1)
  y<-as.data.frame(table(x$motif_id))
  names(y)=c("MOTIF",gsub(".acr/.*","",i))
  if(exists("res")){res=merge(res,y,all.x=T,all.y=T)}else{res=y}
}
head(res)
#      MOTIF jA2.p1k/fimo.tsv jAD1a.p1k/fimo.tsv jAD1d.p1k/fimo.tsv jD5.p1k/fimo.tsv
# 1 MA0001.2            16660              14521              15870            15278
# 2 MA0005.2            15275              13314              14691            14240
# 3 MA0008.2            24341              25420              25707            24890

## compare promoters
x<-res[,-1]
names(x)<-gsub(".p1k/fimo.tsv|j","",names(x))
rownames(x)<-res$MOTIF
write.table(x,"comparePromoter.jaspar2018.txt",sep="\t") 

# repeat for "comparePromoter.athMotif619.txt"
```

Examine promoter results

```R
library("pheatmap")
for(fl in c("comparePromoter.jaspar2018.txt","comparePromoter.athMotif619.txt"))
  {
  print(fl)
  x<-read.table(fl,header=T, sep="\t")
  xx=log2(x)
  # colSums
  print(colSums(x))
  print(colSums(xx))
  # correlations
  print(cor(x))
  print(cor(xx))
  # plot
  pdf(gsub(".txt",".pdf",fl))
  boxplot(x,main="Freq")
  boxplot(xx, main="log2 Freq")
  plot(x,main="Freq")
  pheatmap(cor(x),main="Correlation matrix")
  plot(xx, main="log2 Freq")
  pheatmap(cor(xx),main="Correlation matrix, log2")
  pheatmap(x, scale = "none", fontsize_row=3, main="Euclidean, unscaled")  # default euclidean distance
  pheatmap(xx, scale = "none", fontsize_row=3, main="log2, Euclidean, unscaled")  # default euclidean distance
  pheatmap(xx, scale = "row", fontsize_row=3, main="log2, Euclidean, row scaled")  # default euclidean distance
  pheatmap(xx, scale = "column", fontsize_row=3, main="log2, Euclidean, column scaled")  # default euclidean distance
  dev.off()
}
```

Print Output

```
[1] "comparePromoter.jaspar2018.txt"
     A2    AD1a    AD1d      D5 
6661743 5984436 6226172 6205415 
      A2     AD1a     AD1d       D5 
6215.696 6089.921 6169.104 6135.315 
            A2      AD1a      AD1d        D5
A2   1.0000000 0.9399910 0.9433270 0.9609195
AD1a 0.9399910 1.0000000 0.9859531 0.9740298
AD1d 0.9433270 0.9859531 1.0000000 0.9594149
D5   0.9609195 0.9740298 0.9594149 1.0000000
            A2      AD1a      AD1d        D5
A2   1.0000000 0.9853390 0.9839880 0.9853099
AD1a 0.9853390 1.0000000 0.9970224 0.9955319
AD1d 0.9839880 0.9970224 1.0000000 0.9968052
D5   0.9853099 0.9955319 0.9968052 1.0000000
[1] "comparePromoter.athMotif619.txt"
     A2    AD1a    AD1d      D5 
8793873 8027567 8324093 8131779 
      A2     AD1a     AD1d       D5 
8029.477 7860.305 7949.761 7902.994 
            A2      AD1a      AD1d        D5
A2   1.0000000 0.9579867 0.9653300 0.9600881
AD1a 0.9579867 1.0000000 0.9828834 0.9993615
AD1d 0.9653300 0.9828834 1.0000000 0.9820444
D5   0.9600881 0.9993615 0.9820444 1.0000000
            A2      AD1a      AD1d        D5
A2   1.0000000 0.9872215 0.9867716 0.9854858
AD1a 0.9872215 1.0000000 0.9975637 0.9983553
AD1d 0.9867716 0.9975637 1.0000000 0.9986782
D5   0.9854858 0.9983553 0.9986782 1.0000000
```

Above are based on motif occurance, let us check TF types

```R
info<-read.table("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/Qregulation3/0.input/Ath_TF_binding_motifs_information.txt",header=T, sep="\t")
#
fl="comparePromoter.athMotif619.txt"
x<-read.table(fl,header=T, sep="\t")
f<-info$Family[match(rownames(x),info$Gene_id)]
tf<-aggregate(x,by=list(f),sum)
x<-tf[,-1]
rownames(x)<-tf$Group.1
xx<-log2(x)
#
pdf(gsub(".txt",".TF.pdf",fl))
boxplot(x,main="Freq")
boxplot(xx, main="log2 Freq")
plot(x,main="Freq")
pheatmap(cor(x),main="Correlation matrix")
plot(xx, main="log2 Freq")
pheatmap(cor(xx),main="Correlation matrix, log2")
pheatmap(x, scale = "none", fontsize_row=8, main="Euclidean, unscaled")  # default euclidean distance
pheatmap(xx, scale = "none", fontsize_row=8, main="log2, Euclidean, unscaled")  # default euclidean distance
pheatmap(xx, scale = "row", fontsize_row=8, main="log2, Euclidean, row scaled")  # default euclidean distance
pheatmap(xx, scale = "column", fontsize_row=8, main="log2, Euclidean, column scaled")  # default euclidean distance
dev.off()
```

**Conclusions**: A2 promoter is different from the other three, so the background should be considered for comparing motif enrichment.  

### 2. Obtain enriched motifs for each ACR list

2.1 Run [AME](https://meme-suite.org/meme/tools/ame) for each ACR list, using the corresponding (sub)genome promoter sequences as control to obtain enriched motifs from `Ath_TF_binding_motifs.meme`

```R
info<-read.table("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/Qregulation3/0.input/Ath_TF_binding_motifs_information.txt",header=T, sep="\t")
#
setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/Qregulation3/2.MotifEnrich")
# loop input
fl<-list.files(pattern="tsv"); fl
# "ame_A2.tsv"  "ame_D5.tsv"  "ame_F1a.tsv" "ame_F1d.tsv" 
for(i in fl)
  {
  x<-read.table(i,header=T,sep="\t")
  # x<-x[,c("motif_ID","adj_p.value")]
  x<-x[,c("motif_ID","rank")]
  names(x)[2]<-gsub("ame_|.tsv","",i)
  if(exists("res")){res=merge(res,x,all.x=T,all.y=T)}else{res=x}
}
x<-res[,-1]
rownames(x)<-res$motif_ID
head(x)
# rank
z<-x
z[is.na(z)]<-600 # NA as insignaificant
# binary
y<-x[,c("A2","D5","F1a","F1d","AD1a","AD1d")]
y[!is.na(y)]<-1
y[is.na(y)]<-0
colSums(y)
#  A2   D5  F1a  F1d AD1a AD1d 
# 329  307  274  313  380  388 
dim(y)
# 423   6
family= info$Family[match(rownames(y),info$Gene_id)]
pattern=apply(y[,1:6],1, function(x){paste0(x,collapse="")})
table(pattern)
# 000001 000010 000011 000100 000101 000110 000111 001000 
#     15     13     33      2      2      1      2      1 
# 001011 001111 010000 010010 010011 010101 010111 011000 
#      1      2      1      2      4      3      8      1 
# 011111 100000 100001 100010 100011 100100 100101 100111 
#      3      7      3      2      9      1      2      6 
# 101001 101010 101011 101111 110000 110011 110101 110111 
#      1      1      5      7      1      7      2     23 
# 111011 111100 111110 111111 
#      3      1      1    247 

# cross tabulation of motif family and grouping
mp<-as.matrix(table(family,pattern))
head(mp)

## plot
library(pheatmap)
library(RColorBrewer)
library(UpSetR)
pdf("compareMotifs.pdf")
# color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
pheatmap(z, scale = "none", fontsize_row=3, main="Rank: Euclidean, unscaled",colorRampPalette(rev(brewer.pal(n = 11, name =
  "Spectral")))(100))
pheatmap(cor(z), scale = "none", fontsize_row=8, main="Correlation",colorRampPalette(brewer.pal(n = 8, name =   "Oranges"))(50))
upset(y,nset=ncol(y),order.by="freq",decreasing=TRUE)
upset(y,nset=ncol(y))
pie(table(family),main=paste("Enriched motifs"))
for(i in unique(pattern))
{
  use=family[pattern==i]
  pie(table(use),main=paste(i,length(use)))
}
pheatmap(mp[,colnames(mp)!="111111"],main="Cross tabulation of motif TF family and grouping")
pheatmap(mp,scale="column",main="Cross tabulation of motif TF family and grouping, column scaled")
dev.off()

## extract motifs not common to all
mm<-as.data.frame(pattern[pattern!="111111"])
write.table(mm,"diffEnrichMotifs.txt",sep="\t")

write.table(cbind(x,pattern,family),"EnrichMotifs.txt",sep="\t")
```



2.2 Run [XSTREME](https://meme-suite.org/meme/tools/xstreme) for each ACR list, using the corresponding (sub)genome promoter sequences as control to obtain enriched motifs 

```bash
cd /Users/guanjinghu/Nutstore\ Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/Qregulation3/2.MotifEnrich

# examine output file
ls XSTREME/*N/xstreme.txt 
# 1.XSTREME/A2/xstreme.txt	1.XSTREME/D5/xstreme.txt
# 1.XSTREME/AD1a/xstreme.txt	1.XSTREME/F1a/xstreme.txt
# 1.XSTREME/AD1d/xstreme.txt	1.XSTREME/F1d/xstreme.txt

# check command to make sure input and parameters correct
tail -2 XSTREME/*N/xstreme.tsv

# check number of motif find
for j in $(ls XSTREME/*N/xstreme.tsv);do;
echo $j
tail -5 $j |head -1|cut -f1
done
# 1.XSTREME/A2N/xstreme.tsv		117
# 1.XSTREME/AD1aN/xstreme.tsv	367
# 1.XSTREME/AD1dN/xstreme.tsv	352
# 1.XSTREME/D5N/xstreme.tsv		133
# 1.XSTREME/F1a/xstreme.tsv		3    ??????? what is wrong
# 1.XSTREME/F1d/xstreme.tsv		143

# check number of seed motif find
for j in $(ls XSTREME/*N/xstreme.txt);do;
echo $j
grep -c "^MOTIF" $j
done
# XSTREME/A2N/xstreme.txt 25
# 1.XSTREME/AD1aN/xstreme.txt 47
# 1.XSTREME/AD1dN/xstreme.txt 45
# 1.XSTREME/D5N/xstreme.txt 36
# 1.XSTREME/F1aN/xstreme.txt 3
# 1.XSTREME/F1dN/xstreme.txt 28

# check Background frequency
for j in $(ls XSTREME/*N/xstreme.txt);do;
echo $j
sed -n 27p $j
done
# Background letter frequencies (from file `./background'):
#1.XSTREME/A2N/xstreme.txt   A 0.34720 C 0.15280 G 0.15280 T 0.34720 #1.XSTREME/AD1aN/xstreme.txt A 0.35960 C 0.14040 G 0.14040 T 0.35960 
#1.XSTREME/D5N/xstreme.txt   A 0.35640 C 0.14360 G 0.14360 T 0.35640 
#1.XSTREME/AD1dN/xstreme.txt A 0.35510 C 0.14490 G 0.14490 T 0.35510 
#. see the drop of GC. content in AD1a promoters?? What does it mean by openness?

# merge meme
ls XSTREME/*N/xstreme.txt
cat XSTREME/*N/xstreme.txt |grep -c "^MOTIF" 
# 184
cat XSTREME/*N/xstreme.txt |grep "^MOTIF" |sort |uniq |wc -l 
# 121
meme2meme XSTREME/*N/xstreme.txt >1.XSTREME/merged.txt
grep -c "^MOTIF" XSTREME/merged.txt #184

# manually remove reduancy
grep -c "^MOTIF" XSTREME/merged_nr.txt #121

```

### 3. Cluster enriched motif

3.1 Try STAMP http://www.benoslab.pitt.edu/stamp/ with default

3.2 Run [RSAT Matrix-Clustering](http://rsat.eead.csic.es/plants/matrix-clustering_form.cgi) for XTREME results to cluster all seed motifs obtained

```bash
pwd
# /Users/guanjinghu/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/Qregulation3/2.MotifEnrich
# prep xtreme results for RSAT, must contain the Background line
cat XSTREME/A2N/xstreme.txt | sed s/'^MOTIF '/'MOTIF A2-'/g >../3.Clustering/XTREMEseedFull/xstremeA2.txt
cat XSTREME/D5N/xstreme.txt | sed s/'^MOTIF '/'MOTIF D5-'/g  >../3.Clustering/XTREMEseedFull/xstremeD5.txt
cat XSTREME/F1aN/xstreme.txt | sed s/'^MOTIF '/'MOTIF F1a-'/g  >../3.Clustering/XTREMEseedFull/xstremeF1a.txt
cat XSTREME/F1dN/xstreme.txt | sed s/'^MOTIF '/'MOTIF F1d-'/g  >../3.Clustering/XTREMEseedFull/xstremeF1d.txt
cat XSTREME/AD1aN/xstreme.txt | sed s/'^MOTIF '/'MOTIF AD1a-'/g >../3.Clustering/XTREMEseedFull/xstremeAD1a.txt
cat XSTREME/AD1dN/xstreme.txt | sed s/'^MOTIF '/'MOTIF AD1d-'/g >../3.Clustering/XTREMEseedFull/xstremeAD1d.txt
# meme combined
meme2meme ../3.Clustering/XTREMEseedFull/xstreme*txt >../3.Clustering/XTREMEseedFull/combined.txt
grep -c "MOTIF" ../3.Clustering/XTREMEseedFull/combined.txt
# 184 motifs
```

3.3 Run [RSAT Matrix-Clustering](http://rsat.eead.csic.es/plants/matrix-clustering_form.cgi) for XTREME seed motifs and AME motifs

```bash
cd /Users/guanjinghu/Nutstore\ Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/Qregulation3/3.Clustering
cp ../0.input/Ath_TF_binding_motifs.meme .
```

```R
m<-read.table("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/Qregulation3/2.MotifEnrich/diffEnrichMotifs.txt", header=T, sep="\t",colClasses = 'character')
x<-readLines("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/Qregulation3/3.Clustering/Ath_TF_binding_motifs.meme")
y<-x
ss<-grep("MOTIF",x)
for(i in ss){
  #print(x[i])
  tt<-unlist(strsplit(x[i]," "))
  tt[2]<- paste0(tt[2],"-",m[tt[2],])
  y[i]<-paste(tt,collapse=" ")
}
writeLines(y,"~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/Qregulation3/3.Clustering/Ath_TF_binding_motifs.rename.meme")
```

Result:http://rsat.eead.csic.es/plants/tmp/apache/2021/11/23/matrix-clustering_2021-11-23.141735_M81GmO/



### 4. Scan enriched motifs for comparison

```bash
## Scan 121 enriched motifs
#acr
fimo --o 3.FIMO/A2.acr 1.enrich\ with\ XSTREME/merged_nr.txt 0.input/A2_ACRs_promoter1k.fasta
fimo --o 3.FIMO/D5.acr 1.enrich\ with\ XSTREME/merged_nr.txt 0.input/D5_ACRs_promoter1k.fasta
fimo --o 3.FIMO/F1a.acr 1.enrich\ with\ XSTREME/merged_nr.txt 0.input/F1a_ACRs_promoter1k.fasta
fimo --o 3.FIMO/F1d.acr 1.enrich\ with\ XSTREME/merged_nr.txt 0.input/F1d_ACRs_promoter1k.fasta
fimo --o 3.FIMO/AD1a.acr 1.enrich\ with\ XSTREME/merged_nr.txt 0.input/AD1a_ACRs_promoter1k.fasta
fimo --o 3.FIMO/AD1d.acr 1.enrich\ with\ XSTREME/merged_nr.txt 0.input/AD1d_ACRs_promoter1k.fasta
#promoter
fimo --o 3.FIMO/A2.p1k 1.enrich\ with\ XSTREME/merged_nr.txt 0.input/A2_promoter1k.fasta
fimo --o 3.FIMO/D5.p1k 1.enrich\ with\ XSTREME/merged_nr.txt 0.input/D5_promoter1k.fasta
fimo --o 3.FIMO/AD1a.p1k 1.enrich\ with\ XSTREME/merged_nr.txt 0.input/AD1a_promoter1k.fasta
fimo --o 3.FIMO/AD1d.p1k 1.enrich\ with\ XSTREME/merged_nr.txt 0.input/AD1d_promoter1k.fasta
```

The motif occurance in ACRs were normalized by by total ACR sizes, and then divided by occurance in 1k promoter regions. 

```R
setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/Qregulation3/4.FIMO")

library(data.table)
## loop input FIMO scanning results and sum into table

fl<-paste0(list.files(pattern="acr"),"/fimo.tsv"); fl
# "A2.acr/fimo.tsv"   "AD1a.acr/fimo.tsv" "AD1d.acr/fimo.tsv" "D5.acr/fimo.tsv"   "F1a.acr/fimo.tsv"  "F1d.acr/fimo.tsv" 
for(i in fl)
  {
  x<-fread(i,select=1)
  y<-as.data.frame(table(x$motif_id))
  names(y)=c("MOTIF",gsub(".acr/.*","",i))
  if(exists("res")){res=merge(res,y,all.x=T,all.y=T)}else{res=y}
}

fl<-paste0(list.files(pattern="p1k"),"/fimo.tsv"); fl
#"A2.p1k/fimo.tsv"   "AD1a.p1k/fimo.tsv" "AD1d.p1k/fimo.tsv" "D5.p1k/fimo.tsv" 
for(i in fl)
  {
  x<-fread(i,select=1)
  y<-as.data.frame(table(x$motif_id))
  names(y)=c("MOTIF",gsub("/.*","",i))
  res=merge(res,y,all.x=T,all.y=T)
}
write.table(res,file="rawFreq.txt",sep="\t",row.names=FALSE)

## compare freq
x<-res[,-1]
rownames(x)<-res$MOTIF
colSums(x)
#     A2     AD1a     AD1d       D5      F1a      F1d   A2.p1k AD1a.p1k AD1d.p1k   D5.p1k 
#   31525    67609    56027    21045    22621    24704  1213512   949836  1037081   986457 
boxplot(log2(x[,1:6]),main="Freq in ACR, log2")
boxplot(log2(x[,7:10]),main="Freq in Promoter, log2")

# normalize by total acr length
freq_acr<-x[,1:6]
# acr size 
lz<-c(17894,17826,15950,17345,24044,22074)
names(lz)<-c("A2","D5","F1a","F1d","AD1a","AD1d")
l<-lz[names(freq_acr)]
sizeFactor<-l/mean(l)
freq_acr_norm<-sweep(freq_acr, 2 ,sizeFactor, "/")

## then check enrichment in1kb promoter
y=freq_acr_norm
y$A2 = y$A2/x$A2.p1k
y$F1a = y$F1a/x$A2.p1k
y$D5 = y$D5/x$D5.p1k
y$F1d = y$F1d/x$D5.p1k
y$AD1a = y$AD1a/x$AD1a.p1k
y$AD1d = y$AD1d/x$AD1d.p1k
freqRatio_acr_norm= y

## then check enrichment in1kb promoter
y=freq_acr
y$A2 = y$A2/x$A2.p1k
y$F1a = y$F1a/x$A2.p1k
y$D5 = y$D5/x$D5.p1k
y$F1d = y$F1d/x$D5.p1k
y$AD1a = y$AD1a/x$AD1a.p1k
y$AD1d = y$AD1d/x$AD1d.p1k
freqRatio_acr= y


pdf("plotMotifEnich.pdf")
boxplot(log2(x[,1:6]),main="Freq in ACR, log2")
boxplot(log2(x[,7:10]),main="Freq in Promoter")
#
boxplot(log2(freq_acr),main="Freq in ACR, log2")
heatmap(as.matrix(log2(freq_acr)),scale="none", main="Freq in ACR, log2",cexRow=0.2)
#
boxplot(log2(freq_acr_norm),main="Freq in ACR normalized by size, log2")
heatmap(as.matrix(log2(freq_acr_norm)),scale="none", main="Freq in ACR normalized by size, log2", cexRow=0.2)
#
boxplot(freqRatio_acr,main="Freq ratio of acr/promoter")
heatmap(as.matrix(freqRatio_acr),scale="none", main="Freq ratio of acr/promoter", cexRow=0.2)
#
boxplot(freqRatio_acr_norm,main="Freq ratio of acr_norm/promoter")
heatmap(as.matrix(freqRatio_acr_norm),scale="none", main="Freq ratio of acr_norm/promoter", cexRow=0.2)
dev.off()
```



## Tools

Main tasks include:

* De novo discovery
* Enrichment of known motifs
* Sanning to derive frequency
* Comparison and visualization

**Popular suites**

MEME http://meme-suite.org/index.html

RSAT http://rsat.eead.csic.es/plants/

**Other tools**

* DiNAMO (2018) highly sensitive DNA motif discovery in high-throughput sequencing data https://github.com/bonsai-team/DiNAMO
* DMINDA 2.0 (2017) integrated and systematic views of regulatory DNA motif identification and analyses http://bmbl.sdstate.edu/DMINDA2
* ....



### RSAT Matrix-clustering Usage

**Input**

choose meme format

**Output files**

`matrix-clustering_radial_tree.html`

`matrix-clustering_portable_logo_tree.html`

`matrix-clustering_aligned_logos/`pairwise alignment

`matrix-clustering_html/pairwise_compa.html`summary of pairwise alignment



### XSTREME Usage

**Output files**

`xstreme.html` report page

`xstreme.txt` contains the 'seed' motif from each cluster of similar (redundant) motifs. This column is 1 if this motif is the seed motif of its cluster; otherwise the column is 0. The seed motif in a cluster is either the most enriched discovered motif, or, if the cluster contains no discovered motifs, the most enriched known motif in the cluster. Enrichment is determined by the SEA motif enrichment analysis program. The motifs are in Minimal MEME Motif format.

`xstreme.tsv`  contains one line for each motif found to be significantly enriched, sorted in order of decreasing statistical significance. The first line in the file contains the (tab-separated) names of the fields. Your command line is given at the end of the TSV file in a comment line starting with the character '#'. The names and meanings of each of the fields are described in the table below.
field	name	contents

1. RANK	The rank of the enrichment of the motif in the primary sequences relative to the control sequences according to the SEA motif enrichment analysis program. The rank is relative to all of the ab initio motifs discovered by STREME and MEME, as well as all the known motifs you provided to XSTREME. The rank is based on the SEA enrichment p-values of the motifs. Motifs with tied SEA enrichment p-values are given the same rank. Click on the rank to see the motif in the SEA output. (There will be no link to the SEA output if the SEA E-value of the motif is larger than the E-value threshold.) Caution:The SEA statistical significance statistics (p-, E- and q-values) are NOT accurate for discovered motifs because they are based on the same sequences as input to STREME and MEME.
2. SEED_MOTIF	This column is 1 if this motif is the seed mo tif of its cluster; otherwise the column is 0. The seed motif in a cluster is either the most enriched discovered motif, or, if the cluster contains no discovered motifs, the most enriched known motif in the cluster. Enrichment is determined by the SEA motif enrichment analysis program.
3. CLUSTER	This column contains the cluster number of the cluster of similar motifs that the motif belongs to. Clusters are numbered in order of the enrichment of their seed motif.
4. SOURCE	The name of of the program that found the de novo motif, or the name of a file of motifs ("motif database file") that contains the motif.
5. ID	The name of the motif, which is unique in the motif discovery program output or in the motif database file.
6. ALT_ID	An alternate name for the motif that may be provided in the motif discovery program output or in the motif database file.
7. CONSENSUS	A consensus sequence computed from the motif (as described below).
8. WIDTH	The width of the motif.
9. SITES	The number of primary sequences matching the motif according to the motif discovery program (discovered motifs), or according to the SEA enrichment analysis program (database motifs).
10. SEA_PVALUE	This is the p-value of the motif according to the SEA motif enrichment analysis program. The RANK of the motif is based on this value. This is NOT an accurate measure of the statistical significance of the motif because it is not adjusted for multiple tests. SEA estimates the p-value of motifs using either the Fisher exact test or the Binomial test, as is described here.
11. EVALUE	The E-value of the motif according to the motif discovery program (ab initio motifs), or according to the SEA motif enrichment analysis program (known motifs). It is an accurate measure of the statistical significance of the motif unless the value of EVALUE_ACC is "0". The different programs estimate E-values as follows:
    * SEA estimates the E-value of the motif using either the Fisher exact test or the Binomial test, as is described here.
    * STREME also estimates the E-value of the motif using either the Fisher exact test or the Binomial test, as is described here.
    * MEME estimates the E-value of the motif based on its log likelihood ratio, as is described here. Note: MEME E-values are based on a different null model than STREME and SEA E-values and are therefore not directly comparable to them.
12. EVALUE_ACC	The E-value is an accurate measure of motif significance if this is "1"; otherwise the EVALUE is not an accurate measure of motif significance.
13. SIM_SOURCE	Motifs discovered by a motif discovery program (STREME or MEME) are compared with known motifs in the motif database(s) you specified. This column and the next two columns identify the most similar known motif to the discovered motif. This column gives the similar motif's database name, and the next two columns give its ID and its URL, respectively. Only known motifs with Tomtom similarity E-values of less than 1.0 to the discovered motif will be shown here. For known motifs, this column and the next two give the database name, ID and URL of the known motif.
14. SIM_MOTIF	For discovered motifs this column gives the ID (with any alternate ID in parentheses) of the most similar known motif in the motif database(s) that you provided to XSTREME. For known motifs it gives the known motif's ID.
15. MOTIF_URL	For discovered motifs this columns gives the URL of the most similar known motif in the motif database(s) that you provided to XSTREME. For known motifs it gives the known motif's URL in its online motif database.

A consensus sequence is constructed from each column in a motif's frequency matrix using the "50% rule" as follows:

1. The letter frequencies in the column are sorted in decreasing order.
2. Letters with frequency less 50% of the maximum are discarded.
3. The letter used in this position in the consensus sequence is determined by the first rule below that applies:
   * If there is only one letter left, or if the remaining letters exactly match an ambiguous symbol in the alphabet, the letter or ambiguous symbol, respectively, is used.
   * Otherwise, if the remaining set contains at least 50% of the core symbols in the alphabet, the alphabet's wildcard (e.g., "N" for DNA or RNA, and "X" for protein) is used.
   * Otherwise, the letter with the maximum frequency is used.