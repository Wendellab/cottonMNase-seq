s####### DNS-seq pipline ########
## See detailed notes "iSeg_analysis_notes"

## Pre-installed programs
# module load R
# R

library(travis)
options(threads=8)
options(verbose=T)

## examine MAD and SD across genome, within "iseg/"
fileL = list.files("iseg",pattern="MDS&SD.txt",full.names=TRUE)
res = data.frame(chr=c(1:13,"all"))
for(file in fileL)
{
    x = read.table(file, header=TRUE,sep="\t")$MAD
    if(length(x)>14){
        x = cbind(x[c(1:13,27)],x[14:27])
        colnames(x) = paste0(gsub("iseg/|_.*","",file),c(".A",".D"))
        res = cbind(res,x)
    }else{
        res = cbind(res,x)
        colnames(res)[ncol(res)] = gsub("iseg/|_.*","",file)
    }
}
res
write.table(res, "iseg/sumMDS.txt", sep="\t",row.names=FALSE)
pdf("iseg/sumMDS.pdf")
boxplot(res[,-1], las=2)
dev.off()

## quantile normalization for pooled data, expect A
bgs = c(list.files("mappingA_new",pattern="..Dn.*_qnorm.bg",full.names=TRUE),
list.files("mappingD",pattern=".cD.*_qnorm.bg",full.names=TRUE),
list.files("mappingF",pattern=".cD.*_qnorm.bg",full.names=TRUE),
list.files("mappingMnew",pattern=".cD.*_qnorm.bg",full.names=TRUE))
bgs
# unqual lines normalization?
qbgs=bgQuantileNorm(bgs)
qbgs
# mv files into iseg/
cmd = paste0("mv ", qbgs, " iseg/")
sapply(cmd,system)

## check MDS again
# Check MAD and SD before running iSeg
# define input
fileL = list.files("iseg", pattern= "..D.*_qnorm.bg", full.names=TRUE)
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
    write.table(MS, file=paste0(gsub("_.*","",input),"_MDS&SD.txt"), row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}
fileL = list.files("iseg",pattern="MDS&SD.txt",full.names=TRUE)
res = data.frame(chr=c(1:13,"all"))
for(file in fileL)
{
    x = read.table(file, header=TRUE,sep="\t")$MAD
    if(length(x)>14){
        x = cbind(x[c(1:13,27)],x[14:27])
        colnames(x) = paste0(gsub("iseg/|_.*","",file),c(".A",".D"))
        res = cbind(res,x)
    }else{
        res = cbind(res,x)
        colnames(res)[ncol(res)] = gsub("iseg/|_.*","",file)
    }
}
res
write.table(res, "iseg/sumMDS.quantile.txt", sep="\t",row.names=FALSE)
pdf("iseg/sumMDS.quantile.pdf")
boxplot(res[,-1], las=2)
dev.off()

## now run iSegv190207 v1.3.2
# define input and output
system("mkdir isegv1.3.2")
fileL = list.files("iseg", pattern= "..D.*_qnorm.bg", full.names=TRUE)
fileL
maxwl =100
minwl =10
for(input in fileL){
    outputDir = gsub("iseg/|_.*","",input)
    # run iseg for BC 1, 2, 3
    cmd<-paste0("~/iSegv190207/iSeg -cz -ctp -bc 1.0-2.0-3.0-4.0-5.0-6.0 -minwl ",minwl," -maxwl ",maxwl," -d ",input," -of ",outputDir, " >isegv1.3.2/",outputDir,".manifest.txt")
    message(cmd)
    system(cmd)
    system(paste0("mv ",outputDir,"/* isegv1.3.2/ ; rm -r ",outputDir))
}

q()
n
# module load python
# module load py-numpy
# module load py-pandas
# module load py-openpyxl
# python3 ~/iSegv190207/iSeg.py -cz -ctp -bc 1.0-2.0-3.0-4.0-5.0-6.0 -minwl 10 -maxwl 100 -d iseg/ScD_q20_chr.size_w20_qnorm_qnorm.bg -of ScD >isegv1.3.2/ScD.manifest.txt
# mv ScD/* isegv1.3.2/ ; rm -r ScD

#############################################################################
## now run iSeg v180809b
# define input and output
fileL = list.files("iseg", pattern= "..D.*_qnorm.bg", full.names=TRUE)
fileL
# system("mkdir ../iseg")
# run iseg
# set parameters
maxwl =100
minwl =10
nchr=c(13,13,26,26)
names(nchr) = fileL
for(input in fileL){
    outputDir = gsub("iseg/|_.*","",input)
    # run iseg for BC 1, 2, 3
    cmd<-paste0("~/iSegv180809b/build/iSeg -cz -ctp -nc ",nchr[input]," -bc 4.0-5.0-6.0 -minwl ",minwl," -maxwl ",maxwl," -d ",input," -of ",outputDir, " >iseg/",outputDir,".manifest.txt")
    message(cmd)
    system(cmd)
    system(paste0("mv ",outputDir,"/* iseg/ ; rm -r ",outputDir))
}
# i cannot get the wraper work yet
# module load python
# module load py-numpy
# module load py-pandas
# python ~/iSegv180809b/build/iSeg.py -cz -ctp -nc 13 -bc 1.0-2.0-3.0 -minwl 10 -maxwl 100 -d iseg/A6Dn_q20_chr.size_w20_qnorm_qnorm.bg -of A6Dn >iseg/A6Dn.manifest.txt
# calculate

#############################################################################
# before v180809
# output dir
system("mkdir iseg")
# set parameters
maxwl =100
minwl =10
nchr=c(13,13, 26, 26)
# define input and output
fileL = grep("D.*_qnorm.bg",list.files(),value=TRUE)
fileL
for(input in fileL){
    outputM = paste0("iseg/",gsub("_.*","",input),"_bc")
    # run iseg for BC 1, 2, 3
    for(bc in 1:3)
    {
        cmd<-paste0("~/iSegv180808sp3/iseg -cz -ctp -nc ",nchr," -bc ",bc," -minwl ",minwl," -maxwl ",maxwl," -d ",input," -of ",outputM,bc,".txt >>",outputM,bc,".log")
        message(cmd)
        cat(cmd,file=paste0(outputM,bc,".log"),sep="\n")
        system(cmd)
    }
}

