# https://github.com/gthar/nucleServ/blob/master/bin/nucleR.R

samp
thread=8
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/Nucleosome")
# what portions of genomes are covered by numcleosomes? Well-positioned or fuzzy?
res=c("sample","genome", "N.Fuzzy","N.uncertain","N.well","N.well150","N.fuzzy150", "bp.Fuzzy","bp.uncertain","bp.well","bp.well150","bp.fuzzy150")
for(sample in c("McH","FcH","A6Hn","DcH")){
    print(sample)
    df=fread(paste0(sample,"_q20.nucleosome.gff"),select=c("V1","V9"))
    chr =df$V1
    width =as.numeric(gsub("width=|;.*","",df$V9))
    class =gsub(".*=","",df$V9)
    # a secondary classficiation: Well-positioned (W) as width<150 & score_w>0.4 & score_h>0.6; Otherwise fuzzy (F)
    w150=(width<=150 & class=="W")
    # by chromosome
    res = cbind(t(table(class,by=chr)),t(table(w150,by=chr))[,c("TRUE","FALSE")],t(tapply(width, list(class,chr), sum)), t(tapply(width, list(w150,chr), sum))[,c("TRUE","FALSE")])
    res=rbind(res,c(table(class), table(w150)[c("TRUE","FALSE")],tapply(width, class, sum),tapply(width, w150, sum)[c("TRUE","FALSE")]))
    # F1 and AD1
    if(sample =="McH"|sample=="FcH"){
        genome =gsub("[0-9]","",df$V1)
        # summary by subgenome
        tbl = cbind(t(table(class,by=genome)), t(table(w150,by=genome)[c("TRUE","FALSE"),]),t(tapply(width, list(class,genome), sum)),t(tapply(width, list(w150,genome), sum)[c("TRUE","FALSE"),]))
        res = rbind(res,tbl)
    }
    write.table(res,file=paste0(sample,".sumTbl.nuclosome.txt"),sep="\t")
}
res

library(reshape2)
library(ggplot2)
library(RColorBrewer)
options(stringAsFactor=F)

## nuclesome classification
###########################
cz=read.table("chr.size.txt",sep="\t",header=T)
cz[,-1]=apply(cz[,-1],2,as.numeric)
nA=read.table("A6Hn.sumTbl.nuclosome.txt",header=T,sep="\t"); nA$genome="A"; nA$ploidy="Diploid";
nD=read.table("DcH.sumTbl.nuclosome.txt",header=T,sep="\t"); nD$genome="D"; nD$ploidy="Diploid";
nF=read.table("FcH.sumTbl.nuclosome.txt",header=T,sep="\t"); nF$genome=c(rep(c("A","D"),each=13),NA,"A","D");
nM=read.table("McH.sumTbl.nuclosome.txt",header=T,sep="\t"); nM$genome=c(rep(c("A","D"),each=13),NA,"A","D"); nM$ploidy="AD1"
nA$size=cz$A2
nD$size=cz$D5
nF$ploidy="F1"; nF$size = c(cz$A2[1:13],cz$D5[1:13], cz$A2[14]+cz$D5[14], cz$A2[14], cz$D5[14])
nM$size = c(cz$AD1At[1:13],cz$AD1Dt[1:13], cz$AD1At[14]+cz$AD1Dt[14], cz$AD1At[14], cz$AD1Dt[14])
n=rbind(nA[1:13,],nD[1:13],nF[1:26,],nM[1:26,])
n$NC.W = n$"TRUE..1"/n$size
n$NC.F = n$"FALSE..1"/n$size
n$NC=n$NC.F+n$NC.W
# ANOVA test
a= aov(NC~ploidy+genome, n)
TukeyHSD(a)
# $ploidy
#                     diff           lwr         upr     p adj
# Diploid-AD1  0.021535013  0.0076874044 0.035382622 0.0011075
# F1-AD1       0.014358111  0.0003804749 0.028335747 0.0427341
# F1-Diploid  -0.007176902 -0.0210245111 0.006670707 0.4340266
# $genome
#           diff        lwr        upr p adj
# D-A 0.08461524 0.07516663 0.09406385     0
with(n,aggregate(NC, list(ploidy,genome),function(x)c(mean(x),sd(x)/sqrt(13))))
#  Group.1 Group.2        mean         error
# 1     AD1       A 0.801037967 0.003045319
# 2 Diploid       A 0.782255535 0.005001122
# 3      F1       A 0.788015529 0.005652702
# 4     AD1       D 0.841592288 0.003544322
# 5 Diploid       D 0.899116561 0.001615257
# 6      F1       D 0.883330949 0.001164685
df=n
df$ploidy=factor(df$ploidy, levels=c("Diploid","F1","AD1"))
# plot NC
p=ggplot(df,aes(x=ploidy, y=NC,fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Dark2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6)) + ylab ("Nucleosome Coverage")
pW=ggplot(df,aes(x=ploidy, y=NC.W,fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Dark2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6)) + ylab ("Nucleosome Coverage")+ggtitle("Well")
pF=ggplot(df,aes(x=ploidy, y=NC.F,fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Dark2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6)) + ylab ("Nucleosome Coverage")+ggtitle("Fuzzy")

## NRL summary
##############
library(reshape2)
x<-read.table("estimateNRL.txt",sep="\t",header=T)
df=melt(x[-1,-1])
df$genome = gsub(".*_","",df$variable)
df$ploidy = factor(gsub("_.*","",df$variable),levels=c("Diploid","F1","AD1"))
pNRL = ggplot(df,aes(x=ploidy, y=value,fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Dark2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6))  + ylab("Nucleosome Repeat Length (bp)")
a= aov(value~ploidy+genome, df)
TukeyHSD(a)
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#Fit: aov(formula = NRL ~ ploidy + genome, data = df[df$class == "Estimation", ])
#
#$ploidy
#                 diff        lwr      upr     p adj
#F1-Diploid  0.2423077 -0.6816578 1.166273 0.8055725
#AD1-Diploid 3.1450769  2.2119176 4.078236 0.0000000
#AD1-F1      2.9027692  1.9696099 3.835929 0.0000000
#$genome
#          diff       lwr        upr     p adj
#D-A -0.8669514 -1.499531 -0.2343718 0.0079027


## NUC prediction
A<-read.table("predict.A.txt",sep="\t",header=T); A$genome="A";A$ploidy="Diploid"
D<-read.table("predict.D.txt",sep="\t",header=T); D$genome="D";D$ploidy="Diploid"
M<-read.table("predict.M.txt",sep="\t",header=T); M$genome=rep(c("A","D"),each=13); M$ploidy="AD1"
df=rbind(A,D,M)
df$ploidy = factor(df$ploidy,levels=c("Diploid","F1","AD1"))
#NRL
pNRL.p = ggplot(df,aes(x=ploidy, y=NRL.mean ,fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Set2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6))  + ylab("Nucleosome Repeat Length (bp)")
pNC.p = ggplot(df,aes(x=ploidy, y=NucCov ,fill=genome)) + geom_boxplot() + scale_color_grey()+ scale_fill_brewer(palette="Set2") + theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust=0.6))  + ylab("Nucleosome Coverage")




library(gridExtra)
pdf("nucleosome.pdf")
grid.arrange(pNRL + theme(legend.position = "none")+xlab(""),
             p + theme(legend.position = "none")+xlab(""),
             pNRL.p+ theme(legend.position = "none")+xlab(""),
             pNC.p + theme(legend.position = "none")+xlab(""),
             nrow = 2)
pNRL
p
pNRL.p
pNC.p
pF
pW
dev.off()
