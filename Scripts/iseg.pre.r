####### DNS-seq pipline ########
## See detailed notes "iSeg_analysis_notes"

## Pre-installed programs
# * kent source utilities - "bedGraphToBigWig"
# export PATH="$PATH:~/kent"
# module load R
# module load bedtools2
# R

library(travis)
options(threads=8)
options(verbose=T)

# start with bed files from H and L, get coverage for 20bp windows, quantile normalize, then get difference L-H
allbeds=files("*.bed")
# make sure to use only unparsed data
upbeds = files("*q20.bed")


# make 20bp windows genome wide
chromsizes="chr.size.txt" 
w=bedtoolsMakeWindows(chromsizes,20,genome=T)
w # chr.size_w20.bed

# select beds that are unparsed for iSeg
bgs=bedtoolsCoverage(upbeds,w)
bgs
qbgs=bgQuantileNorm(bgs)
qbgs

# now prepare normlized w20 coverge files for L-H
l=files("*L*qnorm.bg")
h=files("*H*qnorm.bg")
# check to see if light and heavy are paired properly
data.frame(l,h)

# calculate difference between light and heavy
dbgs=bgOps(l,"difference",h,pattern="L",replacement="D")
dbgs
# dbgs files will be used for iSeg

# convert bedGraph files to bigWig
bws=bedGraphToBigWig(c(l,h,dbgs),chromsizes)
# The resulting bigWig files can be used on most genome browsers
bws

# Check MAD and SD before running iSeg
# define input
fileL = grep("..D.*_qnorm.bg",list.files(),value=TRUE)
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
    write.table(MS, file= gsub("_.*","_MAD&SD.txt",input), row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}


q()
n
