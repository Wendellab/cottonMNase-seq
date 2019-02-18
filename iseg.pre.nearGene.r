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

# mappiable ref region, faster with bg
library(data.table)
getEffectiveRefSize.bg=function(bedGraphFile){
    rpm <- fread(bedGraphFile, select = 4)$V4;
    start <- fread(bedGraphFile, select = 2)$V2;
    end <- fread(bedGraphFile, select = 3)$V3;
    ts <- end-start
    cs <- sum(ts[rpm!=0])
    return(cs)
}
# calculate SD and MAD
getMAD<-function(bg, by.chr=FALSE)
{
    df<-fread(bg)
    names(df)<-c("chr","start","end","cov")
    val<-rep(df$cov, df$end-df$start)
    MAD = median(abs(val-median(val)))
    # SD = 1.4826 * MAD
    res<-MAD
    names(res)="all"
    if(by.chr==TRUE)
    {
        chrs = unique(df$chr)
        for(i in chrs){
            df_chr = df[chr==i]
            val<-rep(df_chr$cov, df_chr$end-df_chr$start)
            MAD = median(abs(val-median(val)))
            names(MAD)=i
            res = c(res,MAD)
        }
    }
    return(res)
}
getSD <- function(bg)
{
    val<-fread(bg,select=4)$V4
    return(sd(val))
}



# coverage with reads per million covering 1MB
genomes = c("AD1","F1","A2","D5")
tags=c("M","F","A","D")
dirs = c("mappingMnew","mappingF","mappingA_new","mappingD")
for(i in 1:4)
{
    genome = genomes[i]
    dir=dirs[i]
    setwd(paste0("/work/LAS/jfw-lab/hugj2006/cottonLeaf/",dir))
    
    # start with bed files from H and L, get coverage for20bp windows, quantile normalize, then get difference L-H
    # allbeds=files("*.bed")
    # make sure to use only unparsed data
    # upbeds = files("*q20.bed")
    
    # subset reads to gene regions
    for(i in upbeds)
    {system(paste0("bedtools window -a ",i," -b ../refGenomes/",genome,".og.bed -l 3000 -r 3000 | cut -f1-6 |uniq >",gsub("q20","q20.genic",i)))}
    genbeds = list.files(pattern="q20.genic.bed")
    print(genbeds)
    
    # make 20bp windows genome wide
    chromsizes="chr.size.txt"
    # w=bedtoolsMakeWindows(chromsizes,20,genome=T)
    w="chr.size_w20.bed"
    print(w)
    
    # subset 20 bp windowns to gene regions, some redundant due to overlapped regions between genes
    system(paste0("bedtools window -a ",w," -b ../refGenomes/",genome,".og.bed -l 3000 -r 3000 | cut -f1-3 |uniq >",gsub("w20","w20.genic",w)))
    genw = gsub("w20","w20.genic",w)
    print(genw)
    
    # calculate genomce coverge, NOT rpm
    bgs=bedtoolsCoverage(genbeds,genw,scalar = 1)
    print(bgs)
    
    # get some numbers
    readN = filelines(genbeds, threads = getOption("threads", 1L))
    effSize = sapply(bgs,getEffectiveRefSize.bg)
    refSize=filelines(genw)*20
    print(data.frame(readN,effSize, refSize))
    write.table(data.frame(readN,effSize, refSize), file=paste0("../iseg/", genomes[i],"_genic_par.txt"), row.names=T,col.names=T,quote=FALSE,sep="\t")
    
    # deption correction: effRefSize/totalReadsreadN
    cmd = paste0("awk '{print $1,$2  ,$3 ,$4* ",effSize/readN," }' OFS='\t' ",bgs," > ",gsub(".bg$",".depth.bg",bgs))
    cmdRun(cmd,threads = getOption("threads", 1L))
    sbgs = gsub(".bg$",".depth.bg",bgs)
    # quantiel normalize
    qsbgs=bgQuantileNorm(sbgs)
    print(qsbgs)
    
    # examine SD, see what correction and qnorm did
    print(data.frame(sapply(bgs,getSD), sapply(sbgs,getSD), sapply(qsbgs,getSD) ))
    
    # now prepare normlized w20 coverge files for L-H
    l=grep("L",qsbgs,value=TRUE)
    h=grep("H",qsbgs,value=TRUE)
    # check to see if light and heavy are paired properly
    print(data.frame(l,h))
    
    # calculate difference between light and heavy
    dbgs=bgOps(l,"difference",h,pattern="L",replacement="D")
    print(dbgs)
    # dbgs files will be used for iSeg
    
    #list.files(pattern="genic")
    #file.remove(l)
    
    # convert bedGraph files to bigWig
    bws=bedGraphToBigWig(c(l,h,dbgs),chromsizes)
    # The resulting bigWig files can be used on most genome browsers
    print(bws)
    
    # Check MAD and SD before running iSeg
    # define input
    # dbgs = grep("..D.*genic.depth_qnorm.bg",list.files(),value=TRUE)
    res=cbind(sapply(dbgs,function(x)getMAD(x, by.chr=TRUE)))
    colnames(res) = gsub("_.*","",dbgs)
    print(res)
    write.table(res, file=paste0("../iseg/", genome,"_genic_MAD.txt"), row.names=T,col.names=T,quote=FALSE,sep="\t")
}

# summarize MAD and plot
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/")
# examine all MADs
fileL = paste0("iseg/",genomes,"_genic_MAD.txt")
res = data.frame(chr=c("all",1:13))
for(file in fileL)
{
    x = read.table(file, header=TRUE,sep="\t")
    if(nrow(x)>14){
        y = cbind(x[1:14,],x[c(1,15:27),])
        names(y) = paste(names(y), c(rep("A",ncol(x)),rep("D",ncol(x)) ) , sep=".")
        res = cbind(res,y)
    }else{
        res = cbind(res,x)
    }
}
res
pdf("iseg/sumMAD.genic.pdf")
cols = ifelse(grepl("^A|A$",names(res)[-1],),"pink","skyblue")
cols[grep("^A",names(res)[-1])]="red"
cols[grep("^D",names(res)[-1])]="blue"
boxplot(res[2:14,-1], las=2, col=cols)
boxplot(res[2:14,grep(".c|n",names(res))], las=2, col=cols[grep(".c|n",names(res)[-1])])
dev.off()


# summarize refSize and effSize
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/")
# examine all MADs
fileL = paste0("iseg/",genomes,"_genic_par.txt")
par =list()
for(file in fileL)
{
    x = read.table(file, header=TRUE,sep="\t")
    perc = x$effSize/x$refSize
    par[[file]]=perc
}
par

q()
n