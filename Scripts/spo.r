####### DNS-seq pipline ########
## See detailed notes "iSeg_analysis_notes"

## Pre-installed programs
# * kent source utilities - "bedGraphToBigWig"
# export PATH="$PATH:/work/LAS/jfw-lab/hugj2006/tools/kent"
# module load R
# module load bedtools2
# R

library(travis)
options(threads=8)
options(verbose=T)

#######################################
# make working dir and gather input BED
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/")
system("mkdir isegRes/SPO")
setwd("isegRes/SPO")
genome=c("D","A","F","M")

#######################################
# remove reads over mt region
system("bedtools intersect -v -a ../../mappingD/DcL_q20.bed -b ../../ATAC/blacklist.bed >DcL_q20.bed")
system("wc -l ../../mappingD/DcL_q20.bed") #36048546
system("wc -l DcL_q20.bed") #35963776

#######################################
# prep inpu for coverage
beds = c("DcL_q20.bed", "../../mappingA_WHU/A6Ln_q20.bed","../../mappingF2020/FcL_q20.bed","../../mappingM_UTX/McL_q20.bed")
chromsizes=c("../../mappingD/chr.size.txt", "../../mappingA_WHU/chr.size.txt","../../mappingF2020/chr.size.txt","../../mappingM_UTX/chr.size.txt")
names(chromsizes)=genome
w20=c("../../mappingD/chr.size_w20.bed", "../../mappingA_WHU/chr.size_w20.bed","../../mappingF2020/chr.size_w20.bed","../../mappingM_UTX/chr.size_w20.bed")
names(w20)=genome

#######################################
# parse fragments by length
bpl=bedParseLengths(beds,c(0,130,260))
# awk '{ l=$3-$2;  if(l >=0 && l <130){ print $0 > "DcL_q20_000-130.bed" }if(l >=130 && l <260){ print $0 > "DcL_q20_130-260.bed" } }' OFS='    ' DcL_q20.bed
# awk '{ l=$3-$2;  if(l >=0 && l <130){ print $0 > "A6Ln_q20_000-130.bed" }if(l >=130 && l <260){ print $0 > "A6Ln_q20_130-260.bed" } }' OFS='    ' ../../mappingA_new/A6Ln_q20.bed
# awk '{ l=$3-$2;  if(l >=0 && l <130){ print $0 > "FcL_q20_000-130.bed" }if(l >=130 && l <260){ print $0 > "FcL_q20_130-260.bed" } }' OFS='    ' ../../mappingF/FcL_q20.bed
# awk '{ l=$3-$2;  if(l >=0 && l <130){ print $0 > "McL_q20_000-130.bed" }if(l >=130 && l <260){ print $0 > "McL_q20_130-260.bed" } }' OFS='    ' ../../mappingMnew/McL_q20.bed
names(bpl)=genome
bpl
# $D
# [1] "DcL_q20_000-130.bed" "DcL_q20_130-260.bed"
# $A
# [1] "A6Ln_q20_000-130.bed" "A6Ln_q20_130-260.bed"
# $F
# [1] "FcL_q20_000-130.bed" "FcL_q20_130-260.bed"
# $M
# [1] "McL_q20_000-130.bed" "McL_q20_130-260.bed"

#######################################
# calculate genome coverage and convert to bw for IGV visualization
bg=bedtoolsGenomeCov(bpl$D,chromsizes["D"]); bedGraphToBigWig(bg,chromsizes["D"])
bg=bedtoolsGenomeCov(bpl$A,chromsizes["A"]); bedGraphToBigWig(bg,chromsizes["A"])
bg=bedtoolsGenomeCov(bpl$F,chromsizes["F"]); bedGraphToBigWig(bg,chromsizes["F"])
bg=bedtoolsGenomeCov(bpl$M,chromsizes["M"]); bedGraphToBigWig(bg,chromsizes["M"])

#######################################
# make 21 bp window by sliding 5 bp
w=c()
w["D"]=bedtoolsMakeWindows(chromsizes["D"],21, 5, genome=T,outnames="Dchr.size_w21s5.bed")
w["A"]=bedtoolsMakeWindows(chromsizes["A"],21, 5, genome=T,outnames="Achr.size_w21s5.bed")
w["F"]=bedtoolsMakeWindows(chromsizes["F"],21, 5, genome=T,outnames="Fchr.size_w21s5.bed")
w["M"]=bedtoolsMakeWindows(chromsizes["M"],21, 5, genome=T,outnames="Mchr.size_w21s5.bed")
# bedtools makewindows -w 21 -s 5 -g ../../mappingD/chr.size.txt | awk 'x != $3; {x=$3}' > Dchr.size_w21s5.bed
# bedtools makewindows -w 21 -s 5 -g ../../mappingA_new/chr.size.txt | awk 'x != $3; {x=$3}' > Achr.size_w21s5.bed
# bedtools makewindows -w 21 -s 5 -g ../../mappingF/chr.size.txt | awk 'x != $3; {x=$3}' > Fchr.size_w21s5.bed
# bedtools makewindows -w 21 -s 5 -g ../../mappingMnew/chr.size.txt | awk 'x != $3; {x=$3}' > Mchr.size_w21s5.bed
w
#                     D                     A                     F                    M
#"Dchr.size_w21s5.bed" "Achr.size_w21s5.bed" "Fchr.size_w21s5.bed" "Mchr.size_w21s5.bed"

########################################
# get midpoint of small light fragment, remmeber to sort
bedMidPoint=function(bed){
    out = gsub(".bed",".mid.bed",bed)
    cmd=paste0("cat ",bed,"|awk '{ l=($3+$2)/2; print $1,int(l),int(l)+1,$4,$5,\"+\" }' OFS='\t' | sortBed -i >",out)
    print(cmd)
    system(cmd)
    return(out)
}
bmps = c()
bmps["A"] = bedMidPoint(bpl$A[1])
bmps["D"] = bedMidPoint(bpl$D[1])
bmps["F"] = bedMidPoint(bpl$F[1])
bmps["M"] = bedMidPoint(bpl$M[1])
# cat A6Ln_q20_000-130.bed|awk '{ l=($3+$2)/2; print $1,int(l),int(l)+1,$4,$5,"+" }' OFS='\t' | sortBed -i >A6Ln_q20_000-130.mid.bed
# cat DcL_q20_000-130.bed|awk '{ l=($3+$2)/2; print $1,int(l),int(l)+1,$4,$5,"+" }' OFS='\t' | sortBed -i >DcL_q20_000-130.mid.bed
# cat FcL_q20_000-130.bed|awk '{ l=($3+$2)/2; print $1,int(l),int(l)+1,$4,$5,"+" }' OFS='\t' | sortBed -i >FcL_q20_000-130.mid.bed
# cat McL_q20_000-130.bed|awk '{ l=($3+$2)/2; print $1,int(l),int(l)+1,$4,$5,"+" }' OFS='\t' | sortBed -i >McL_q20_000-130.mid.bed
bmps
#                        A                          D                         F                          M
#"A6Ln_q20_000-130.mid.bed"  "DcL_q20_000-130.mid.bed"  "FcL_q20_000-130.mid.bed"  "McL_q20_000-130.mid.bed"

#######################################
# intersect w20s5 windows with midpoint, then make coverge bedgraph for MOA frenters
A = bedtoolsIntersect(w["A"], bmps["A"],extraargs = "-wa")
D = bedtoolsIntersect(w["D"], bmps["D"],extraargs = "-wa")
F = bedtoolsIntersect(w["F"], bmps["F"],extraargs = "-wa")
M = bedtoolsIntersect(w["M"], bmps["M"],extraargs = "-wa")
# bedtools intersect -a Achr.size_w21s5.bed -b A6Ln_q20_000-130.mid.bed -wa > Achr.size_w21s5_x_A6Ln_q20_000-130.mid.bed
bgA=bedtoolsCoverage(A,w20["A"])
bgD=bedtoolsCoverage(D,w20["D"])
bgF=bedtoolsCoverage(F,w20["F"])
bgM=bedtoolsCoverage(M,w20["M"])
# bedtools coverage -a Achr.size_w21s5_x_A6Ln_q20_000-130.mid.bed -b ../../mappingA_new/chr.size_w20.bed > Achr.size_w21s5_x_A6Ln_q20_000-130.mid_chr.size_w20.bg
bgs= c(bgA,bgD,bgM,bgF)

#######################################
# make bigwig for MOA frenters viualization
bedGraphToBigWig(bgA,chromsizes["A"])
bedGraphToBigWig(bgD,chromsizes["D"])
bedGraphToBigWig(bgF,chromsizes["F"])
bedGraphToBigWig(bgM,chromsizes["M"])
# bedGraphToBigWig Achr.size_w21s5_x_A6Ln_q20_000-130.mid_chr.size_w20.bg ../../mappingA_new/chr.size.txt Achr.size_w21s5_x_A6Ln_q20_000-130.mid_chr.size_w20.bw
# bedGraphToBigWig Dchr.size_w21s5_x_DcL_q20_000-130.mid_chr.size_w20.bg ../../mappingD/chr.size.txt Dchr.size_w21s5_x_DcL_q20_000-130.mid_chr.size_w20.bw
# bedGraphToBigWig Fchr.size_w21s5_x_FcL_q20_000-130.mid_chr.size_w20.bg ../../mappingF/chr.size.txt Fchr.size_w21s5_x_FcL_q20_000-130.mid_chr.size_w20.bw
# bedGraphToBigWig Mchr.size_w21s5_x_McL_q20_000-130.mid_chr.size_w20.bg ../../mappingMnew/chr.size.txt Mchr.size_w21s5_x_McL_q20_000-130.mid_chr.size_w20.bw

#######################################
# quantile normalization to obtain consistent MAD for iSeg
# qbgs=bgQuantileNorm(bgs)
# Not a good practice, because the nonzero region can be significantly changed before and after qnorm.

#######################################
# Check MAD and SD before running iSeg
library(data.table)
getMADnSD<-function(df_bg)
{
    val<-rep(df_bg$cov, df_bg$end-df_bg$start)
    MAD = median(abs(val-median(val)))
    SD = 1.4826 * MAD
    #aSD = sd(val)
    u=mean(val)
    aSD=sqrt(sum((val - u)^2)/length(val))
    res<-c(MAD,SD,aSD)
    names(res)=c("MAD","SD","sd")
    return(res)
}
for(input in bgs){
    bg<-fread(input)
    names(bg)<-c("chr","start","end","cov")
    chrs<-unique(bg$chr)
    MS<-c("MAD","SD","sd")
    for(i in chrs){
        MS<-rbind(MS,getMADnSD(bg[chr==i]))
    }
    MS<-rbind(MS,getMADnSD(bg))
    MS<-cbind(c("chr",chrs,"all"),MS)
    write.table(MS, file= gsub("chr.*","_MAD&SD.txt",input), row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
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
isegSummary = function(inputBGs,isegDir,BGsub=".*/|_.*", SDs=NULL){
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
        flag = gsub(BGsub,"",bg)
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
    resL = list.files(isegDir,pattern="bc.*txt")
    ## loop import and report stats for each flag
    flags = unique(gsub(".*/|_bc.*","",resL))
    for(flag in flags){
        # collect files under same flag
        flagFiles = resL[grep(flag,resL)]
        if(is.null(SDs))
        {
            stat = grep("SD",readLines(paste0(isegDir,"/",flag,".manifest.txt")),value=TRUE )
            SD=as.numeric(gsub(".* ","",stat))
            print(flag)
            print(stat)
        }else{SD=SDs[gsub("_.*","",flag)]}
        # prep result table
        flag=gsub("_.*","",flag)
        res = data.frame(ChrName=names(sl[[flag]]), Length = sl[[flag]], Chromosome = 1:length(sl[[flag]]))
        res.nz = data.frame(ChrName=names(nonzero[[flag]]), Length = nonzero[[flag]], Chromosome = 1:length(nonzero[[flag]]))
        res.peak = data.frame(ChrName=names(sl[[flag]]), Length = sl[[flag]], Chromosome = 1:length(sl[[flag]]))
        # loop each file
        for(file in flagFiles){
            bc = gsub(".*bc|[.]txt","",file)
            # read into grange
            filePath = paste0(isegDir,"/",file)
            gr = readGeneric(file=filePath, chr = 1, start = 2, end = 3, strand = NULL, meta.col=list(height=4,tstat=5,pvalue=6), sep=" ")
            gr$pvalue=fread(filePath,sep=" ",select=6)$V6
            gr$dns = ifelse(gr$height>0,"Pos","Neg")
            ## reduce range to count peaks
            s =reduce(gr[gr$dns=="Pos"]); s$dns = "Pos"
            r =reduce(gr[gr$dns=="Neg"]); if(length(r)>0){r$dns = "Neg"}
            rgr = c(s,r)
            ## Region bps called as peaks
            sRegion = data.frame(tapply(width(rgr), list(seqnames(rgr), rgr$dns),sum))
            names(sRegion) =paste0(names(sRegion),".bc",bc)
            ## Proportion of regions called
            pRegion = sweep(sRegion,1, nonzero[[flag]],"/")
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
            df <- data.frame(seqnames=seqnames(srgr), starts=start(srgr)-1, ends=end(srgr), names="Peak", scores=0, strands="+", starts=start(srgr), ends=end(srgr), color = ifelse(srgr$dns=="Pos", "20,20,255","255,20,20"))
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
## now run iSeg
setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/isegRes/SPO")
# MADs
sds=c()
sds["A"] = as.numeric(gsub(".*\t","",readLines("A_MAD&SD.txt")[15]))
sds["D"] = as.numeric(gsub(".*\t","",readLines("D_MAD&SD.txt")[15]))
sds["F"] = as.numeric(gsub(".*\t","",readLines("F_MAD&SD.txt")[28]))
sds["M"] = as.numeric(gsub(".*\t","",readLines("M_MAD&SD.txt")[28]))
bgs = list.files( pattern= "mid_chr.size_w20.bg")
names(bgs) =gsub("chr.*","",bgs)
maxwl =100
minwl =10
isegDir = "iseg_v1.3.4_060520"
cat("# run iseg\nmodule load boost/1.69.0-py3-openmpi3-yc4whx2\nexport PATH=$PATH:/work/LAS/jfw-lab/hugj2006/tools/isegv190624/iSeg", file="runiSeg.sh",sep="\n")
cat(paste0("mkdir ",isegDir), file="runiSeg.sh",sep="\n",append=TRUE)
genome=c("D","A","F","M")
for(i in genome){
    input= bgs[i]
    sd=sds[i]
    outputDir = gsub("chr.*","_frenters",input)
    # run iseg for BC 1, 2, 3
    cmd1<-paste0("iSeg -cz -ctp -bc 4.0-4.5-5.0-5.5-6.0-6.5 -minwl ",minwl," -maxwl ",maxwl," -sd ",sd," -d ",input," -of ",outputDir, " >",isegDir,"/",outputDir,".manifest.txt")
    cmd2 <- paste0("mv ",outputDir,"/* ",isegDir,"/; rm -r ",outputDir)
    cat(cmd1,cmd2, file="runiSeg.sh",sep="\n", append=TRUE)
}
system("cat runiSeg.sh")
system("bash runiSeg.sh")
isegSummary(inputBGs=bgs ,isegDir = isegDir, BGsub="chr.*", SDs=sds)
## inspect "plotDis.*.pdf" for each BC of each sample
## inspect "sumStat.*.txt" for each sample
q("no")
