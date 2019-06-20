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

#######################################
# make working dir and gather input BED
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/")
system("mkdir isegRes/MOA")
setwd("isegRes/MOA")
genome=c("D","A","F","M")

#######################################
# remove reads over mt region
system("bedtools intersect -v -a ../../mappingD/DcL_q20.bed -b ../../ATAC/blacklist.bed >DcL_q20.bed")
system("wc -l ../../mappingD/DcL_q20.bed") #36048546
system("wc -l DcL_q20.bed") #35963776

#######################################
# prep inpu for coverage
beds = c("DcL_q20.bed", "../../mappingA_new/A6Ln_q20.bed","../../mappingF/FcL_q20.bed","../../mappingMnew/McL_q20.bed")
chromsizes=c("../../mappingD/chr.size.txt", "../../mappingA_new/chr.size.txt","../../mappingF/chr.size.txt","../../mappingMnew/chr.size.txt")
names(chromsizes)=genome
w=c("../../mappingD/chr.size_w20.bed", "../../mappingA_new/chr.size_w20.bed","../../mappingF/chr.size_w20.bed","../../mappingMnew/chr.size_w20.bed")
names(w)=genome

#######################################
# parse fragments by length
bpl=bedParseLengths(beds,c(0,130,260))
names(bpl)=genome

#######################################
# calculate genome coverage and convert to bw for IGV visualization
bg=bedtoolsGenomeCov(bpl$D,chromsizes["D"]); bedGraphToBigWig(bg,chromsizes["D"])
bg=bedtoolsGenomeCov(bpl$A,chromsizes["A"]); bedGraphToBigWig(bg,chromsizes["A"])
bg=bedtoolsGenomeCov(bpl$F,chromsizes["F"]); bedGraphToBigWig(bg,chromsizes["F"])
bg=bedtoolsGenomeCov(bpl$M,chromsizes["M"]); bedGraphToBigWig(bg,chromsizes["M"])

#######################################
# make bedgraph for 0-130bp coverage over 20bp windown.
A=bedtoolsCoverage(bpl$A[1],w["A"])
D=bedtoolsCoverage(bpl$D[1],w["D"])
F=bedtoolsCoverage(bpl$F[1],w["F"])
M=bedtoolsCoverage(bpl$M[1],w["M"])
bgs=c(A,D,F,M)

#######################################
# quantile normalization to obtain consistent MAD for iSeg
qbgs=bgQuantileNorm(bgs)
#####I am not sure if it is the best practice, because the nonzero region significantly changed before and after qnorm.

#######################################
# Check MAD and SD before running iSeg
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
for(input in qbgs){
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


#######################################
## run and summarize iSeg results
# plot function
log10hist <- function(x, main = filePath, xlab="",breaks = "Sturges")
{
    h<-hist(x,plot=F,breaks = breaks);
    h$counts[h$counts>0] <- log10(h$counts[h$counts>0])
    plot(h, ylab='log10(Freq)', xlab =xlab,main=main)
}
# summary function: generate summary table, plot each seg file by chr, convert output to BED
isegSummary = function(inputBGs,isegDir){
    print(inputBGs)
    print(isegDir)
    library(data.table)
    library(genomation)
    library(GenomicRanges)
    # chromosome length and non-zero region
    sl=list()
    nonzero = list()
    for(bg in inputBGs)
    {
        flag = gsub(".*/|_.*","",bg)
        dt=fread(bg,sep="\t")
        sl[[flag]] = tapply(dt$V3,dt$V1,max)
        dt$width = dt$V3-dt$V2
        select= (dt$V4!=0 )
        x = aggregate(dt$width[select],by=list(dt$V1[select]),sum)
        nonzero[[flag]] =  x$x
        names(nonzero[[flag]]) =x$Group.1
    }
    # chr length
    print(sl)
    # nonzero region length
    print(nonzero)
    # get stats
    resL = list.files(isegDir,pattern="0.txt")
    ## loop import and report stats for each flag
    flags = unique(gsub(".*/|_.*","",resL))
    for(flag in flags){
        stat = grep("SD",readLines(paste0(isegDir,"/",flag,".manifest.txt")),value=TRUE )
        print(flag)
        print(stat)
        SD=as.numeric(gsub(".* ","",stat))
        # collect files under same flag
        flagFiles = resL[grep(flag,resL)]
        # prep result table
        res = data.frame(ChrName=names(sl[[flag]]), Length = sl[[flag]], Chromosome = 1:length(sl[[flag]]))
        res.nz = data.frame(ChrName=names(nonzero[[flag]]), Length = nonzero[[flag]], Chromosome = 1:length(nonzero[[flag]]))
        res.peak = data.frame(ChrName=names(sl[[flag]]), Length = sl[[flag]], Chromosome = 1:length(sl[[flag]]))
        # loop each file
        for(file in flagFiles){
            bc = gsub(".*bc|[.].*","",file)
            # read into grange
            filePath = paste0(isegDir,"/",file)
            gr = readGeneric(file=filePath, chr = 1, start = 2, end = 3, strand = NULL, meta.col=list(height=4,tstat=5,pvalue=6), sep=" ")
            gr$pvalue=fread(filePath,sep=" ",select=6)$V6
            gr$dns = ifelse(gr$height>0,"Pos","Neg")
            ## Region bps called as peaks
            sRegion = data.frame(tapply(width(gr), list(seqnames(gr), gr$dns),sum))
            names(sRegion) =paste0(names(sRegion),".bc",bc)
            ## Proportion of regions called
            pRegion = sweep(sRegion,1, nonzero[[flag]],"/")
            ## reduce range to count peaks
            s =reduce(gr[gr$dns=="Pos"]); s$dns = "Pos"
            r =reduce(gr[gr$dns=="Neg"]); if(length(r)>0){r$dns = "Neg"}
            rgr = c(s,r)
            ## number of peaks
            nPeaks = data.frame(tapply(as.character(seqnames(rgr)), list(seqnames(rgr), rgr$dns),length))
            names(nPeaks) =paste0(names(nPeaks),".bc",bc)
            ## make report table
            res = cbind(res,sRegion)
            res.nz = cbind(res.nz,pRegion)
            res.peak = cbind(res.peak,nPeaks)     
            ## plot segment width, using reduced
            pdf(paste0(isegDir,"/plotDis.",gsub(".txt","",file),".pdf"))
            # hist of segment width
            log10hist(width(rgr), main = "Segment Width", xlab="bp")
            par(mfrow=c(4,4 ))
            rgrl <- split(rgr, seqnames(rgr))
            for(i in names(rgrl)) log10hist(width(rgrl[[i]]),main=i,xlab="bp")
            ## hist of height
            par(mfrow=c(1,1))
            log10hist(gr$height/SD, main = "Segment Biological Cutoff", xlab="BC", breaks=100)
            par(mfrow=c(4,4))
            grl <- split(gr, seqnames(gr))
            for(i in names(grl)) log10hist(grl[[i]]$height/SD,main=i,xlab="BC", breaks=100)
            # Making karyograms and circos plots
            # par(mfrow=c(1,1))
            # autoplot(seqinfo(gr)) + layout_karyogram(gr, geom="point",aes(x=start,y=height,color=dns))
            dev.off() 
            ## output BED
            # sort rgr
            srgr = sort(rgr)
            # to bed
            df <- data.frame(seqnames=seqnames(srgr), starts=start(srgr), ends=end(srgr), names="Peak", scores=0, strands="+", starts=start(srgr), ends=end(srgr), color = ifelse(srgr$dns=="Pos", "20,20,255","255,20,20"))
            write.table(df, file=paste0(isegDir,"/",gsub(".txt","",file),".Fus.bed"), quote=F, sep="\t", row.names=F, col.names=F)
        } 
        cname = c(names(res)[1:3],sort(names(res[-(1:3)])))
        out=paste0(isegDir,"/sumStat.",flag,".txt")
        ## write results
        write.table("# segment bp",file=out,row.names=FALSE, sep="\t",quote=FALSE)
        write.table(res[,cname],file=out, row.names=FALSE, sep="\t",quote=FALSE, append=TRUE)
        write.table("# segment bp %",file=out, row.names=FALSE, sep="\t",quote=FALSE, append=TRUE)
        write.table(res.nz[,cname], file=out, row.names=FALSE, sep="\t",quote=FALSE, append=TRUE)
        write.table("# segment number",file=out, row.names=FALSE, sep="\t",quote=FALSE, append=TRUE)
        write.table(res.peak[,cname],file=out, row.names=FALSE, sep="\t",quote=FALSE, append=TRUE)
    }
}

#######################################
## now run iSegv190207 v1.3.2 on bigram
## Qnorm input, MAD = 0.0474014 based on MAD calculation
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/isegRes/MOA")
fileL = list.files( pattern= "qnorm.bg", full.names=TRUE)
fileL
maxwl =100
minwl =10
isegDir = "iseg_v1.3.2_qnorm"
cat("# run iseg", file="runiSeg.sh",sep="\n")
cat(paste0("mkdir ",isegDir), file="runiSeg.sh",sep="\n",append=TRUE)
for(input in fileL){
    outputDir = gsub(".*/|_.*","",input)
    # run iseg for BC 1, 2, 3
    cmd1<-paste0("~/iSegv190207/iSeg -cz -ctp -bc 1.0-2.0-3.0-4.0-5.0-6.0 -minwl ",minwl," -maxwl ",maxwl," -d ",input," -of ",outputDir, " >",isegDir,"/",outputDir,".manifest.txt")
    cmd2 <- paste0("mv ",outputDir,"/* ",isegDir,"/; rm -r ",outputDir)
    cat(cmd1,cmd2, file="runiSeg.sh",sep="\n", append=TRUE)
}
system("cat runiSeg.sh")
system("bash runiSeg.sh")
isegSummary(inputBGs=fileL ,isegDir = "iseg_v1.3.2_qnorm")


############ input without qnorm 
fileL = list.files( pattern= "w20.bg", full.names=TRUE)
fileL
maxwl =100
minwl =10
isegDir = "iseg_v1.3.2_041219"
cat("# run iseg", file="runiSeg.sh",sep="\n")
cat(paste0("mkdir ",isegDir), file="runiSeg.sh",sep="\n",append=TRUE)
for(input in fileL){
    outputDir = gsub(".*/|_.*","",input)
    # run iseg for BC 1, 2, 3
    cmd1<-paste0("~/iSegv190207/iSeg -cz -ctp -bc 1.0-2.0-3.0-4.0-5.0-6.0 -minwl ",minwl," -maxwl ",maxwl," -d ",input," -of ",outputDir, " >",isegDir,"/",outputDir,".manifest.txt")
    cmd2 <- paste0("mv ",outputDir,"/* ",isegDir,"/; rm -r ",outputDir)
    cat(cmd1,cmd2, file="runiSeg.sh",sep="\n", append=TRUE)
}
system("cat runiSeg.sh")
system("bash runiSeg.sh")
isegSummary(inputBGs=fileL ,isegDir = "iseg_v1.3.2_041219")
# need higher BC for A genome
system("~/iSegv190207/iSeg -cz -ctp -bc 7.0-8.0 -minwl 10 -maxwl 100 -d ./A6Ln_q20_000-130_chr.size_w20.bg -of A6Ln >>iseg_v1.3.2_041219/A6Ln.manifest.txt")
isegSummary(inputBGs="./A6Ln_q20_000-130_chr.size_w20.bg" ,isegDir = "A6Ln")
