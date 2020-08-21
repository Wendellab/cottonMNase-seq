####### DNS-seq pipline ########
## See detailed notes "iSeg_analysis_notes"

## Pre-installed programs
# module load R
# R

library(travis)
options(threads=8)
options(verbose=T)

## examine MAD and SD across genome, within "iseg/"
fileL =c(list.files("mappingA_WHU",pattern="MAD&SD.txt", full.names=TRUE), list.files("mappingD",pattern="MAD&SD.txt", full.names=TRUE), list.files("mappingF2020",pattern="MAD&SD.txt", full.names=TRUE), list.files("mappingM_UTX",pattern="MAD&SD.txt", full.names=TRUE))
res = data.frame(chr=c(1:13,"all"))
for(file in fileL)
{
    x = read.table(file, header=TRUE,sep="\t")$MAD
    if(length(x)>14){
        x = cbind(x[c(1:13,27)],x[14:27])
        colnames(x) = paste0(gsub(".*/|_.*","",file),c(".A",".D"))
        res = cbind(res,x)
    }else{
        res = cbind(res,x)
        colnames(res)[ncol(res)] = gsub(".*/|_.*","",file)
    }
}
res
write.table(res, "isegRes/sumMAD.txt", sep="\t",row.names=FALSE)
pdf("isegRes/sumMAD.pdf")
boxplot(res[,-1], las=2)
dev.off()

## quantile normalization for pooled data, expect A
bgs = c(list.files("mappingA_WHU",pattern="..Dn.*_qnorm.bg",full.names=TRUE),
list.files("mappingD",pattern=".cD.*_qnorm.bg",full.names=TRUE),
list.files("mappingF2020",pattern=".cD.*_qnorm.bg",full.names=TRUE),
list.files("mappingM_UTX",pattern=".cD.*_qnorm.bg",full.names=TRUE))
bgs
# unqual lines normalization?
qbgs=bgQuantileNorm(bgs)
qbgs
# mv files into iseg/
cmd = paste0("mv ", qbgs, " isegRes/")
sapply(cmd,system)

## Check MAD and SD again before running iSeg
# define input
fileL = list.files("isegRes", pattern= "..D.*_qnorm.bg", full.names=TRUE)
fileL
# system("mkdir ../iseg")
# calculate SD and MAD
library(data.table)
getMADnSD<-function(df_bg)
{
    val<-rep(df_bg$cov, df_bg$end-df_bg$start)
    MAD = median(abs(val-median(val)))
    SD = 1.4826 * MAD
    res<-c(MAD,SD)
    names(res)=c("MAD","SD")
    return(res)
}
for(input in fileL){
    bg<-fread(input)
    names(bg)<-c("chr","start","end","cov")
    chrs<-unique(bg$chr)
    MS<-c("MAD","SD")
    for(i in chrs){
        MS<-rbind(MS,getMADnSD(bg[chr==i]))
    }
    MS<-rbind(MS,getMADnSD(bg))
    MS<-cbind(c("chr",chrs,"all"),MS)
    write.table(MS, file=paste0(gsub("_.*","",input),"_qnorm_MAD&SD.txt"), row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}
fileL = list.files("isegRes",pattern="_qnorm_MAD&SD.txt",full.names=TRUE)
res = data.frame(chr=c(1:13,"all"))
for(file in fileL)
{
    x = read.table(file, header=TRUE,sep="\t")$MAD
    if(length(x)>14){
        x = cbind(x[c(1:13,27)],x[14:27])
        colnames(x) = paste0(gsub("isegRes/|_.*","",file),c(".A",".D"))
        res = cbind(res,x)
    }else{
        res = cbind(res,x)
        colnames(res)[ncol(res)] = gsub("isegRes/|_.*","",file)
    }
}
res
write.table(res, "isegRes/sumMAD.qnorm.txt", sep="\t",row.names=FALSE)
pdf("isegRes/sumMAD.qnorm.pdf")
boxplot(res[,-1], las=2)
dev.off()


## now run iSegv190207 v1.3.2 on bigram
# define input and output
setwd("isegRes")
isegDir = "iseg_v1.3.4_062020"
cat("# run iseg\nmodule load boost/1.69.0-py3-openmpi3-yc4whx2\nexport PATH=$PATH:/work/LAS/jfw-lab/hugj2006/tools/isegv190624/iSeg", file="runiSeg.sh",sep="\n")
cat(paste0("mkdir ",isegDir), file="runiSeg.sh",sep="\n",append=TRUE)
fileL = list.files(".", pattern= "..D.*_qnorm.bg", full.names=TRUE)
fileL
maxwl =100
minwl =10
for(input in fileL){
    outputDir = gsub("_.*","",input)
    cmd1<-paste0("iSeg -cz -ctp -bc 4.0-4.5-5.0-5.5-6.0-6.5-7.0 -minwl ",minwl," -maxwl ",maxwl," -d ",input," -of ",outputDir, " >",isegDir,"/",outputDir,".manifest.txt")
    cmd2 <- paste0("mv ",outputDir,"/* ",isegDir,"/; rm -r ",outputDir)
    cat(cmd1,cmd2, file="runiSeg.sh",sep="\n", append=TRUE)
}
system("cat runiSeg.sh")
system("bash runiSeg.sh")

q("no")
# module load python
# module load py-numpy
# module load py-pandas
# module load py-openpyxl
# python3 ~/iSegv190207/iSeg.py -cz -ctp -bc 1.0-2.0-3.0-4.0-5.0-6.0 -minwl 10 -maxwl 100 -d iseg/ScD_q20_chr.size_w20_qnorm_qnorm.bg -of ScD >isegv1.3.2/ScD.manifest.txt
# mv ScD/* isegv1.3.2/ ; rm -r ScD




#############################################################################

## now run iSeg v180809b
#  ~/iSegv180809b/build/iSeg -cz -ctp -nc 13 -bc 1.0-2.0-3.0 -minwl 10 -maxwl 100 -d iseg/A6Dn_q20_chr.size_w20_qnorm_qnorm.bg -of A6Dn >iseg/A6Dn.manifest.txt

# i cannot get the wraper work yet
# module load python
# module load py-numpy
# module load py-pandas
# python ~/iSegv180809b/build/iSeg.py -cz -ctp -nc 13 -bc 1.0-2.0-3.0 -minwl 10 -maxwl 100 -d iseg/A6Dn_q20_chr.size_w20_qnorm_qnorm.bg -of A6Dn >iseg/A6Dn.manifest.txt

