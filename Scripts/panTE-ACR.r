####################
## ACR panTE annotation ##
####################
# module load r-udunits2
# module load r/3.5.0-py2-ufvuwmm
library(ChIPseeker) # require r-udunits2
library(GenomicFeatures)
library(ggthemes)
library(gridExtra)
library(reshape2)
options(stringsAsFactors=FALSE)


### ACR, combined MSF and SPO, genomic and TE annotation
load("Qregulation3/ACRsC.rdata")
#
load("refGenomes/TEannotation.rdata")
# "gr.F1"  "gr.A2"  "gr.D5"  "gr.AD1"
gr.F1a=gr.F1[grep("^A",seqnames(gr.F1))]
gr.F1d=gr.F1[grep("^D",seqnames(gr.F1))]
gr.AD1a=gr.AD1[grep("^A",seqnames(gr.AD1))]
gr.AD1d=gr.AD1[grep("^D",seqnames(gr.AD1))]

# loop
for(g in c("A2","D5","F1a","F1d","AD1a","AD1d")){
    print(g)
    peaks=get(paste0("ACR.",g))
    if(g=="A2"){seqlevels(peaks)<-gsub("Chr","A",seqlevels(peaks))}
    if(g=="D5"){seqlevels(peaks)<-gsub("Chr","D",seqlevels(peaks))}
    tes=get(paste0("gr.",g))
    seqlevels(peaks) <- seqlevelsInUse(peaks)
    seqlevels(tes) <- seqlevelsInUse(tes)
    # TE categorization
    tes$class <- gsub(".*Classification=|;Sequence_ontology.*","",tes$family)
    tes$TEid <- gsub(".*Name=|;Classification.*","",tes$family)
    # annotate ACRs with TE overlap
    hit=findOverlaps(peaks,tes)
    d=data.frame(hit)
    # TE class
    d$class =tes$class[subjectHits(hit)]
    type =  aggregate(d$class,list(d$queryHits),function(x)paste0(unique(sort(x)),collapse="@@"))
    peaks$TEclass =NA
    peaks$TEclass[type$Group.1] =type$x
    # TE panID
    d$TEid =tes$TEid[subjectHits(hit)]
    type =  aggregate(d$TEid,list(d$queryHits),function(x)paste0(unique(sort(x)),collapse="@@"))
    peaks$TEid =NA
    peaks$TEid[type$Group.1] =type$x
    assign(paste0("ACRte.",g),peaks)
    
    # annotation by family, how each TE annotated peak correspond to multiple TE families
    number = dcast(d,queryHits~TEid,length)
    rownames(number) =number$queryHits
    number=number[,-1]
    number[number>1]=1
    # ACR proportion co in TE families, whether certain TE families enriched by peaks
    en = aggregate(width(tes),list(tes$TEid),sum)
    ACRinTEsize  = apply(number,2,function(x)sum(x*width(peaks)[as.numeric(rownames(number))]))
    names(en) =c("TE_family","Bp")
    en$TE_family=as.character(en$TE_family)
    en$ACR_in_TE.Bp = 0
    en$ACR_in_TE.Bp = ACRinTEsize[en$TE_family]
    en$ACR_in_TE.Bp[is.na(en$ACR_in_TE.Bp)]=0
    print(cor.test(en$Bp,en$ACR_in_TE.Bp ))
    en1<-en[grep("^TE",en$TE_family,invert=FALSE),]
    nn<-nrow(en1)+1
    en1[nn,2:3]<-colSums(en[grep("^TE",en$TE_family,invert=TRUE),2:3])
    en1[nn,1]<-"other"
    print(cor.test(en1$Bp,en1$ACR_in_TE.Bp ))
    write.table(en1,file=paste0("ACR.panID.",g,".txt"),row.names=FALSE,sep="\t" )
}
   
save(ACRte.A2,ACRte.D5,ACRte.F1a,ACRte.F1d,ACRte.AD1a,ACRte.AD1d,"Qregulation3/ACRsC.rdata")



## assemble together
setwd("~/Downloads")
setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/7.Questions/TErelated")
## Input panTE ids
y<-read.table("panTEids.txt",sep="\t")
head(y)
#               id classL  classS
#      TE_00000000    DNA     DTA
#      TE_00000001    DNA     DTC
#  TE_00000002_LTR    LTR   Gypsy
#  TE_00000003_LTR    LTR   Gypsy
#  TE_00000004_LTR    LTR unknown
#  TE_00000005_LTR    LTR   Gypsy
rownames(y)<-y$id
table(y$classS)
#   Copia      DTA      DTC      DTH      DTM      DTT    Gypsy Helitron  unknown
#    3095      881      402      144     1390       41    16028       12     6064
dim(y) # 28057     4


class=y$classS
class=gsub("Helitron","nonTIR/Helitron",class)
class=gsub("DTC","TIR/CACTA",class)
class=gsub("DTM","TIR/Mutator",class)
class=gsub("DTH","TIR/PIF_Harbinger",class)
class=gsub("DTT","TIR/Tc1_Mariner",class)
class=gsub("DTA","TIR/hAT",class)
class=gsub("Gypsy","LTR/Gypsy",class)
class=gsub("Copia","LTR/Copia",class)
class=gsub("unknown","LTR/unknown",class)
y$class=class
z<-y[,c("id","class")]


for(g in c("A2","D5","F1a","F1d","AD1a","AD1d"))
{
    print(g)
    
    # read panTE bp and TE-derived ACRs
    x<-read.table(paste0("ACR.panID.",g,".txt"),sep="\t",header=T)
    rownames(x)<-x$TE_family
    # Now, use the phyper() function to conduct the hypergeometric test, enriched TEid
    x$enrichP<-NA
    pop_size=sum(x$Bp)
    pop_successes=sum(x$ACR_in_TE.Bp)
    for(i in 1:nrow(x) ){
        sample_successes =x$ACR_in_TE.Bp[i]
        sample_size=x$Bp[i]
        x$enrichP[i] <- phyper(q = sample_successes, m = pop_successes, n = pop_size - pop_successes, k = sample_size, lower.tail = FALSE)
            }
    x$enrichP.bh<-p.adjust(x$enrichP,"BH")
    # get enriched only
    print(dim(xx<-x[x$enrichP.bh<0.05,]))
    
    # check top10 enriched panTEs
    print(xx[order(xx$ACR_in_TE.Bp,decreasing=T),][1:10,])
    z[,g]<-xx[y$id,"ACR_in_TE.Bp"]
}
head(z)
#    id       class   A2   D5  F1a  F1d AD1a AD1d
# 1     TE_00000000     TIR/hAT   NA  152   NA   91   39   NA
# 2     TE_00000001   TIR/CACTA 3739 1268 2725 2527 1300 1755
# 3 TE_00000002_LTR   LTR/Gypsy   NA   NA   NA   NA   NA   NA
# 4 TE_00000003_LTR   LTR/Gypsy  462   NA  496   NA  359   NA
# 5 TE_00000004_LTR LTR/unknown   NA   NA   NA   NA  174   NA
# 6 TE_00000005_LTR   LTR/Gypsy   NA   NA   NA   NA   NA   NA
z[is.na(z)]=0

# filter IDs with significant contribution to ACRs
m<-z[,3:8]
m[is.na(m)]=0
m[m>0]=1
s<-which(rowSums(m)>0)

# plot panTE quantitative contribution to ACR bp
use<-as.matrix(z[s,3:8])
dim(use) # 8680    6
rownames(use)<-z$id[s]
pdf("te-ACR.bp.pdf")
boxplot(log2(use+1), main="TE-derived ACRs: bp by panTE id, log2")
heatmap(log2(as.matrix(use)+1),scale="none",main="Freq in ACR, log2",cexRow=0.2)
dev.off()


m<-m[s,]
colSums(m)
#   A2   D5  F1a  F1d AD1a AD1d
# 4270 2484 4615 2685 3545 2622
TEc<-y$class[s]
library(ggplot2)
library(ComplexUpset)
# library(UpSetR)
# upset(m,set=ncol(m),order.by="freq",decreasing=TRUE)
pdf("te-ACR.pav.pdf",width=11,height=7)
upset(
    m,
    colnames(m),
    annotations = list(
        'Percentage'=(
            ggplot(mapping=aes(fill=TEc))
            + geom_bar(stat='count', position='fill')
            + scale_y_continuous(labels=scales::percent_format())
            + ylab('TE superfamily')
            + scale_fill_brewer(palette="Paired")
        )
    ),
    width_ratio=0.2
)
upset(
    m,
    colnames(m),
    annotations = list(
        'Percentage'=(
            ggplot(mapping=aes(fill=TEc))
            + geom_bar(stat='count', position='fill')
            + scale_y_continuous(labels=scales::percent_format())
            + ylab('TE superfamily')
            + scale_fill_brewer(palette="Paired") + theme(legend.position="none")
        )
    ),
    width_ratio=0.15,
    min_size=50
)
dev.off()


rownames(use)  #8680
load('~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/4.RNAseqAnalysis/IndivGNT2021rRNA/expression.IndivGNTrRNA.rdata')->l;l
# "A2"   "D5"   "F1"   "AD1"  "A2f"  "D5f"  "F1f"  "AD1f"

# check order
summary(D5f$tpmTE[,1:2]==A2f$tpmTE[,1:2])
summary(D5f$tpmTE[,1:2]==AD1f$tpmTE[,1:2])
summary(D5f$tpmTE[,1:2]==F1f$tpmTE[,1:2])
idx<-D5f$tpmTE[,1:2]
tpm<-cbind(A2f$tpmTE[,4:6], D5f$tpmTE[,4:5], F1f$tpmTE[,c(4:6,8:10)], AD1f$tpmTE[,c(4:6,8:10)])
table(filter<-rowSums(tpm)>0) #7056
colSums(tpm)
# SD5-A2-S1      SD5-A2-S4      SD5-A2-S5      SD5-D5-S1      SD5-D5-S3 SD5-A2xD5-S1.x SD5-A2xD5-S2.x
# 26676.138      23413.417      22952.990      63300.319      41726.233       9330.054       9060.534
# SD5-A2xD5-S3.x SD5-A2xD5-S1.y SD5-A2xD5-S2.y SD5-A2xD5-S3.y SD5-Maxxa-S1.x SD5-Maxxa-S2.x SD5-Maxxa-S3.x
#  13171.480      11438.332      11015.529      15114.518      19821.401      23096.534      25623.108
# SD5-Maxxa-S1.y SD5-Maxxa-S2.y SD5-Maxxa-S3.y
# 18623.288      20457.506      22907.616
idx<-idx[filter,]
idx$class= gsub("DNA|MITE","TIR",idx$Classification)
table(idx$class)
#   LTR/Copia    LTR/Gypsy  LTR/unknown      TIR/DTA      TIR/DTC      TIR/DTH      TIR/DTM
#         953         3143         1407          399          263           85          759
#   TIR/DTT TIR/Helitron
#          26          10

# how many express TEs are enriched ACR sources?
table(idx$Name %in% rownames(use))
# FALSE  TRUE
# 2423  4622
table(rownames(use) %in% idx$Name  )


tpm<-tpm[filter,]
info = data.frame(sample=names(tpm), species0 = gsub("SD5-|-S.*","",names(tpm)), species=c(rep("A2",3), rep("D5",2),rep("F1",6),rep("AD1",6)), ploidy=c(rep("diploid",5), rep("F1",6),rep("AD1",6)), subgenome=c(rep("A",3), rep("D",2),rep(c("A","D"),2,each=3)) )

## draw complex heapmap check relationships between different patterns
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)

# head annotation
ha = HeatmapAnnotation(df=info[,c("ploidy","subgenome")], col = list(ploidy=structure(brewer.pal(3, "Set1"), names = c("diploid","F1","AD1")), subgenome=structure(c("brown","purple"), names = c("A","D"))))
# Row annotation
ra = rowAnnotation(class=idx$class, col = list( class=structure(brewer.pal(9, "Paired"),names = levels(factor(idx$class)))))


# heatmap of log2tpm, not scaled
mat<- as.matrix(log2(tpm+1))
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
nkm=4
hm<-Heatmap(mat, name = "TE expression",col = hmcol,top_annotation = ha, right_annotation=ra, row_km=nkm, cluster_columns = FALSE,column_title = "log2(TPM+1)", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = FALSE)
pdf("plotTE.tpm.Fig.pdf")
par(mfrow=c(1,1))
draw(hm)
c<-row_order(hm)
pcol= structure(brewer.pal(9, "Paired"),names = levels(factor(idx$class)))
par(mfrow=c(3,4))
for(i in 1:nkm)
{
    pie(table(idx$class[c[[i]]]), col=pcol,clockwise=TRUE,main=i,labels="")
    barplot(as.matrix(rowsum(mat[c[[i]],], group=idx$class[c[[i]]])), col=pcol,las=2, main="log2(TPM+1)")
}
dev.off()
