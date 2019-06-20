###############################################################
################ Aggregation and Visulization #################
###############################################################
# exploit aggregated profiles over genomic features (gene body, TSS, TTS, repeats, etc.)
# make comparison across genomes and ploidy

# use speedy
# module load py-deeptools
# module load r/3.5.0-py2-ufvuwmm
# module load bedtools2
# cd /work/LAS/jfw-lab/hugj2006/cottonLeaf
# R

#################
## Prep BigWig ##
#################
## BAM to BigWig
# I was debating whether to use
cfile= "aggregation/command.bam2bw.sh"
cat("## BAM to Bigwig\nmodule load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
dirs = c("mappingA_new","mappingD","mappingF","mappingMnew")
# effSize = c()
for(dir in dirs)
{
    print(dir)
    bams=list.files(dir, pattern="*bam$",full=T)
    for(bam in bams){
        # bam to bigwig
        out = paste0("aggregation/", gsub(".*/|.bam","",bam))
           # 1. full range
        cmd = paste0("bamCoverage -b ", bam," -o ",out,".full.bw --binSize 20 --normalizeUsingRPKM")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
           # 2. size 0 to 130
        cmd = paste0("bamCoverage -b ", bam," -o ",out,".0-130.bw --binSize 20 --normalizeUsingRPKM --minFragmentLength 0 --maxFragmentLength 130")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
           # 3. size 130 to 260
        cmd = paste0("bamCoverage -b ", bam," -o ",out,".130-260.bw --binSize 20 --normalizeUsingRPKM --minFragmentLength 130 --maxFragmentLength 260")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
           # 4. nucleosome
        cmd = paste0("bamCoverage -b ", bam," -o ",out,".mnase.bw --binSize 1 --normalizeUsingRPKM --MNase") # --MNase sets centerReads to be 3
        cat(cmd, file =cfile,sep="\n", append=TRUE)
    }
}
system(paste0("cat ",cfile))
system(paste0("bash ",cfile))


## Full occupancy profile: BigWig D = L-H
light = list.files("aggregation",pattern="L.*full.bw",full=T)
heavy = list.files("aggregation",pattern="H.*full.bw",full=T)
diff = gsub("H","D",heavy)
# check matching
df=data.frame(light,heavy,diff)
df
cat("### D = L-H", file ="aggregation/command.getD.sh",sep="\n")
for(i in 1:nrow(df)){
    cmd = paste0("bigwigCompare -b1 ",df$light[i]," -b2 ",df$heavy[i]," -o ",df$diff[i]," --ratio subtract")
    cat(cmd, file ="aggregation/command.getD.sh",sep="\n", append=TRUE)
}
system("cat aggregation/command.getD.sh")
system("bash aggregation/command.getD.sh")

## Merge replicates
cat("### Merge replicates\nmkdir aggregation/combo", file ="aggregation/command.mergeRep.sh",sep="\n")
for(flag in c("D","M","F")){
    path=list.files("aggregation",pattern=paste0("^",flag,".*bw"),full=T)
    type=gsub(".*/..|.bw","",path)
    reps = split(path,type)
    for(t in 1:length(reps))
    {
        b1=reps[[t]][1]
        b2=reps[[t]][2]
        bout=paste0("aggregation/combo/",flag,"c",names(reps)[t],".bw")
        cmd = paste0("bigwigCompare -b1 ",b1," -b2 ",b2," -o ",bout," --ratio mean")
        cat(cmd, file ="aggregation/command.mergeRep.sh",sep="\n", append=TRUE)
    }
}
system("cat aggregation/command.mergeRep.sh")
system("bash aggregation/command.mergeRep.sh")

## partition FcD_q20.full.bw to A and D chromosomes, each qnorm with diploid tracks
#export PATH="$PATH:/work/LAS/jfw-lab/hugj2006/tools/kent/"
system("bash aggregation/command.parseFcD.sh")
library(travis)
setwd("aggregation/parseF/")
bgs =c("AcD_q20.full.bg","FcD_q20.full.chrA.bg","DcD_q20.full.bg","FcD_q20.full.chrD.bg")
qbgs=bgQuantileNorm(bgs)
bws=bedGraphToBigWig(c("FcD_q20.full.chrD_qnorm.bg","DcD_q20.full_qnorm.bg"),chromsizes="../../mappingD/chr.size.txt")
bws
bws=bedGraphToBigWig(c("FcD_q20.full.chrA_qnorm.bg","AcD_q20.full_qnorm.bg"),chromsizes="../../mappingA_new/chr.size.txt")
bws
system("bigwigCompare -b1 FcD_q20.full.chrA_qnorm.bw -b2 AcD_q20.full_qnorm.bw --binSize 20 --ratio log2 -o FAtvsA2.log2ratio.bw")
system("bigwigCompare -b1 FcD_q20.full.chrD_qnorm.bw -b2 DcD_q20.full_qnorm.bw --binSize 20 --ratio log2 -o FDtvsD5.log2ratio.bw")



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
        id = gr$gene_id
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

ref = data.frame(genome = c("A2", "D5","F1","AD1"),
pattern=c("^A","^D","^F","^M"),
txDB =  list.files("refGenomes",pattern="txdb",full.names =T),
subgenome =c(1,1,2,2))

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
    
    if(i %in% 1:2){
        writeGRtoBED(genes.gr, group=NULL, out=paste0("refGenomes/", ref$genome[i], ".gene"))
    }
    if(i %in% 3:4){
        writeGRtoBED(genes.gr[grep("^A",seqnames(genes.gr)),], group=NULL,out=paste0("refGenomes/", ref$genome[i],".gene.A"))
        writeGRtoBED(genes.gr[grep("^D",seqnames(genes.gr)),], group=NULL,out=paste0("refGenomes/", ref$genome[i],".gene.D"))
    }
}
list.files("refGenomes",pattern="gene.*bed")
# "A2.gene.bed"    "AD1.gene.A.bed" "AD1.gene.D.bed" "D5.gene.bed"
# "F1.gene.A.bed"  "F1.gene.D.bed"

################################################
## prepare gene bed files for ortho-homoelogs ##
################################################

library(GenomicRanges)
library(GenomicFeatures)
# FUN - prepare quantile BED files
writeGRtoBED = function(gr, group="quantile",out="")
{
    df <- data.frame(seqnames=seqnames(gr), starts=start(gr)-1, ends=end(gr), names=c(rep(".", length(gr))),scores=c(rep(".", length(gr))),strands=strand(gr))
    if(!is.null(group)){
        id = gr$gene_id
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
# /lss/research/jfw-lab/Projects/MNase-seq/orthohomoeologQuadruplets101218.txt
ogQ<-read.table("orthohomoeologQuadruplets101218.txt", sep="\t", header=TRUE)
head(ogQ)

txdb=loadDb("refGenomes/txdb.A2du.sqlite")
genes.gr = genes(txdb)
genes.gr = genes.gr[genes.gr$gene_id %in% ogQ$A2,]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/A2.og")

txdb=loadDb("refGenomes/txdb.D5.sqlite")
genes.gr = genes(txdb)
genes.gr = genes.gr[genes.gr$gene_id %in% ogQ$D5,]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/D5.og")

txdb=loadDb("refGenomes/txdb.TM1saski.sqlite")
genes.gr = genes(txdb)
genes.gr = genes.gr[genes.gr$gene_id %in% c(as.character(ogQ$At),as.character(ogQ$Dt)),]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/AD1.og")

txdb=loadDb("refGenomes/txdb.F1.sqlite")
genes.gr = genes(txdb)
genes.gr = genes.gr[genes.gr$gene_id %in% c(as.character(ogQ$A2),as.character(ogQ$D5)),]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/F1.og")

list.files("refGenomes",pattern="*og.bed")
# "A2.og.bed"  "AD1.og.bed" "D5.og.bed"  "F1.og.bed"

#############################
## Aggregation per species ##
#############################
library(GenomicRanges)
#library(RColorBrewer)
#library(ggbio)
library(GenomicFeatures)
library(genomation)
options(scipen=999)

## FUN
# extract genomic regions from txDB
getTSS=function(txdb, include=NULL, exclude.pattern=NULL, be=1000, af=1000){
    tss =  genes(txdb)
    if(!is.null(exclude.pattern)) {tss = tss[grep(exclude.pattern,tss$gene_id,invert=T)]}
    if(!is.null(include)) {tss = tss[tss$gene_id %in% include]}
    # end of the + strand genes must be equalized to start pos
    end(tss[strand(tss)=="+",])  =start(tss[strand(tss)=="+",])
    # start of the - strand genes must be equalized to end pos
    start(tss[strand(tss)=="-",])=end(tss[strand(tss)=="-",])
    # define flanking region
    start(tss)=start(tss)-be
    end(tss)  = end(tss)+ af
    # remove duplicated TSSes ie alternative transcripts
    # this keeps the first instance and removes duplicates
    tss=tss[!duplicated(tss),]
    seqlevels(tss)=unique(as.character(seqnames(tss)))
    return(tss)
}
getTTS=function(txdb, include=NULL, exclude.pattern=NULL, be=1000, af=1000){
    tts =  genes(txdb)
    if(!is.null(exclude.pattern)) {tts = tts[grep(exclude.pattern,tts$gene_id,invert=T)]}
    if(!is.null(include)) {tts = tts[tts$gene_id %in% include]}
    # start of the + strand genes must be equalized to end pos
    start(tts[strand(tts)=="+",])  =end(tts[strand(tts)=="+",])
    # end of the - strand genes must be equalized to start pos
    end(tts[strand(tts)=="-",])=start(tts[strand(tts)=="-",])
    # define flanking region
    start(tts)=start(tts)-be
    end(tts)  = end(tts)+ af
    # remove duplicated TSSes ie alternative transcripts
    # this keeps the first instance and removes duplicates
    tts=tts[!duplicated(tts),]
    seqlevels(tts)=unique(as.character(seqnames(tts)))
    return(tts)
}
plotTSSnTTS=function(bw.files,tss.gr,tts.gr,tag,title)
{
    #TSS
    sml = ScoreMatrixList(target = bw.files, windows = tss.gr, bin.num=50, type="bigWig",strand.aware=TRUE)
    names(sml) = tag
    plotMeta(sml,  xcoords = c(-1000, 1000), profile.names=names(sml), ylab="RPKM",dispersion = "se", xlab="TSS", main=title)
    abline(v=0,lty=2)
    #TTS
    sml = ScoreMatrixList(target = bw.files, windows = tts.gr, bin.num=50, type="bigWig",strand.aware=TRUE)
    names(sml) = tag
    plotMeta(sml,  xcoords = c(-1000, 1000), profile.names=names(sml), ylab="RPKM",dispersion = "se", xlab="TTS", main=title)
    abline(v=0,lty=2)
}

## Prepare genomic annotation
ref = data.frame(genome = c("A2", "D5","F1","AD1"), pattern=c("^A","^D","^F","^M"), txDB =  list.files("refGenomes",pattern="txdb",full.names =T), subgenome =c(1,1,2,2), nuc = list.files("nucleosome",pattern="*coverage.bw", full=TRUE) )

# write gene annotation to bed
pdf("aggregation/plotPerGenome.pdf")
for(i in 1:4)
{
    # genomic features
    print(ref[i,])
    txdb = loadDb(as.character(ref$txDB[i]))
    tss.gr =  getTSS(txdb)
    tts.gr =  getTTS(txdb)
    #print(tss.gr)
    #print(tts.gr)
    
    # list all BWs
    bws = list.files("aggregation",pattern=as.character(ref$pattern[i]),full=T)
    
    # full occupancy
    title= paste0(ref[i,],": nucleosome occupancy")
    bw.files = grep("full",bws,value=TRUE)
    print(bw.files)
    tag = gsub("aggregation/|_.*","",bw.files)
    plotTSSnTTS(bw.files,tss.gr,tts.gr,tag,title)
    
    # small size fractioned occupancy
    title= paste0(ref[i,],": 0-130bp fragment occupancy")
    bw.files = grep("0-130",bws,value=TRUE)
    print(bw.files)
    tag = gsub("aggregation/|_.*","",bw.files)
    plotTSSnTTS(bw.files,tss.gr,tts.gr,tag,title)
    
    # big size fractioned occupancy
    title= paste0(ref[i,],": 130-260bp fragment occupancy")
    bw.files = grep("130-260",bws,value=TRUE)
    print(bw.files)
    tag = gsub("aggregation/|_.*","",bw.files)
    plotTSSnTTS(bw.files,tss.gr,tts.gr,tag,title)
    
    # nucleosome
    title= paste0(ref[i,],": 3bp nucleosome centers")
    bw.files = grep("mnase",bws,value=TRUE)
    print(bw.files)
    tag = gsub("aggregation/|_.*","",bw.files)
    plotTSSnTTS(bw.files,tss.gr,tts.gr,tag,title)
}
dev.off()
