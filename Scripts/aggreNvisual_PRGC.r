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


cat("## BAM to Bigwig\nmodule load py-deeptools\nmodule load bedtools2", file="aggregation/command.bam2bw.sh",sep="\n")
dirs = c("mappingA_new","mappingD","mappingF","mappingMnew")
# effSize = c()
for(dir in dirs)
{
    print(dir)
    bams=list.files(dir, pattern="*bam$",full=T)
    for(bam in bams){
        
        # get effective size
        cmd = paste0("bedtools genomecov -ibam ", bam," -d |awk '$3>0 {print $3}' |wc -l >",dir," /effectiveSize.txt")
        message(cmd); cat(cmd, file ="aggregation/command.log.txt",sep="\n", append=TRUE); system(cmd)
        size = read.table("effectiveSize.txt")[1,1]
        #names(size) = bam
        #effSize = c(effSize,size)
        
        # bam to bigwig
        out = paste0("aggregation/", gsub(".*/|.bam","",bam))
           # 1. full range
        cmd = paste0("bamCoverage -b ", bam," -o ",out,".full.bw --binSize 20 --normalizeTo1x ",size)
        message(cmd); cat(cmd, file ="aggregation/command.log.txt",sep="\n", append=TRUE); system(cmd)
           # 2. size 0 to 130
        cmd = paste0("bamCoverage -b ", bam," -o ",out,".0-130.bw --binSize 20 --normalizeTo1x ",size," --minFragmentLength 0 --maxFragmentLength 130")
        message(cmd); cat(cmd, file ="aggregation/command.log.txt",sep="\n", append=TRUE); system(cmd)
           # 3. size 130 to 260
        cmd = paste0("bamCoverage -b ", bam," -o ",out,".130-260.bw --binSize 20 --normalizeTo1x ",size," --minFragmentLength 130 --maxFragmentLength 260")
        message(cmd); cat(cmd, file ="aggregation/command.log.txt",sep="\n", append=TRUE); system(cmd)
           # 4. nucleosome
        cmd = paste0("bamCoverage -b ", bam," -o ",out,".mnase.bw --binSize 1 --normalizeTo1x ",size," --MNase") # --MNase sets centerReads to be 3
        message(cmd); cat(cmd, file ="aggregation/command.log.txt",sep="\n", append=TRUE); system(cmd)
    }
}
system("mv aggregation/command.log.txt aggregation/bam2bw.sh")
system("cat runDeeptools.DNS.sh")
system("bash runDeeptools.DNS.sh")


## Full occupancy profile: BigWig D = L-H
light = list.files("aggregation",pattern="L.*full.bw",full=T)
heavy = list.files("aggregation",pattern="H.*full.bw",full=T)
diff = gsub("H","D",heavy)
# check matching
df=data.frame(light,heavy,diff)
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
    plotMeta(sml,  xcoords = c(-1000, 1000), profile.names=names(sml), ylab="reads per genomic content (1x normalization)",dispersion = "se", xlab="TSS", main=title)
    abline(v=0,lty=2)
    #TTS
    sml = ScoreMatrixList(target = bw.files, windows = tts.gr, bin.num=50, type="bigWig",strand.aware=TRUE)
    names(sml) = tag
    plotMeta(sml,  xcoords = c(-1000, 1000), profile.names=names(sml), ylab="reads per genomic content (1x normalization)",dispersion = "se", xlab="TTS", main=title)
    abline(v=0,lty=2)
}

## Prepare genomic annotation
ref = data.frame(genome = c("A2", "D5","F1","AD1"),
pattern=c("^A","^D","^F","^M"),
txDB =  list.files("refGenomes",pattern="txdb",full.names =T),
subgenome =c(1,1,2,2),
nuc = list.files("nucleosome",pattern="*coverage.bw", full=TRUE)  # nucleosome positions
)

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


##############################
## Comparing across genomes ##
##############################


A2 vs D5, At vs Dt -
A

## Genomation
library(genomation)
for(i in 1:4)
{
    i=2
    print(ref[i,])
    txdb = loadDb(as.character(ref$txDB[i]))
    tss.gr =  getTSS(txdb)
    bw.files = as.character(as.matrix(ref[i,c("nuc","occD","occH","occL")]))
    sml = ScoreMatrixList(target = bw.files, windows = tss.gr, bin.num=50, type="bigWig",strand.aware=TRUE)
    names(sml) = c("nuc","occD","occH","occL")
    # sml.scaled<-scaleScoreMatrixList(sml,column=FALSE,row=TRUE)
    # multiHeatMatrix(sml.scaled, xcoords = c(-1000, 1000))
    # plotMeta(sml.scaled, xcoords = c(-1000, 1000))
    
    multiHeatMatrix(sml,winsorize=c(1,99), order=TRUE,legend.name="tpm",xlab="TSS", xcoords=-1500:1500, grid=FALSE, cex.legend=0.6,cex.lab=0.8, cex.axis=0.6)
    plotMeta(sml,  xcoords = c(-1000, 1000), profile.names=names(sml), ylab="nucleosome occupancy score (tpm)",dispersion = "se", xlab="TSS", main=paste0(ref$genome[i],": Meta-profile"))
    abline(v=0,lty=2)
    plotMeta(sml.scaled,  xcoords = c(-1000, 1000), profile.names=names(sml), ylab="nucleosome occupancy score (tpm)",dispersion = "se", xlab="TSS", main=paste0(ref$genome[i],": Meta-profile, scaled"))
    abline(v=0,lty=2)
    
    
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








sml = ScoreMatrixList(target = "McD_q20_chr.size_w20_qnorm_qnorm.bw", windows = promoter, type="bigWig",strand.aware=TRUE)
# names used by profile.names and matrix.main
names(sml)<-gsub("_u.*|_q20","",names(sml))
sml.scaled<-scaleScoreMatrixList(sml,column=FALSE,row=TRUE)
multiHeatMatrix(sml,winsorize=c(1,99), order=TRUE,legend.name="tpm",xlab="TSS",
xcoords=-1500:1500, grid=FALSE,
cex.legend=0.6,cex.lab=0.8,
cex.axis=0.6)
plotMeta(sml,  xcoords = c(-1500,1500), profile.names=names(sml),
ylab="nucleosome occupancy score (tpm)",dispersion = "se",
xlab="bases around TSS", main=paste0(r,": Meta-profile"))

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

ogQ<-read.table("/lss/research/jfw-lab/Projects/MNase-seq/orthohomoeologQuadruplets101218.txt", sep="\t", header=TRUE)
head(ogQ)

txdb=loadDb("refGenomes/txdb.A2du.sqlite")
genes.gr = genes(txdb)
genes.gr = genes.gr[genes.gr$gene_id %in% ogQ$A2,]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/A2.og")

txdb=loadDb("refGenomes/txdb.D5.sqlite")
genes.gr = genes(txdb)
genes.gr = genes.gr[genes.gr$gene_id %in% ogQ$D5,]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/D5.og")

txdb=loadDb("refGenomes/txdb.TM1saski.sqlit")
genes.gr = genes(txdb)
genes.gr = genes.gr[genes.gr$gene_id %in% c(as.character(ogQ$At),as.character(ogQ$Dt)),]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/AD1.og")

txdb=loadDb("refGenomes/txdb.F1.sqlit")
genes.gr = genes(txdb)
genes.gr = genes.gr[genes.gr$gene_id %in% c(as.character(ogQ$A2),as.character(ogQ$D5)),]
writeGRtoBED(genes.gr, group=NULL, out="refGenomes/F1.og")


#########################################################################################
## Inspect ortho-homoeolog quadruplets for the presence of MSFs in 1kb promoter region ##
#########################################################################################

l = load("iseg/genomicAnnotation.byType.rdata")



## collect gene containing MSFs in 1kb promoter
getGenes=function(csAnno, feature="Promoter (<=1kb)")
{
    anno = csAnno@anno
    genes = anno$geneId[anno$annotation==feature]
    return(genes)
}
geneL = lapply(MSFpeakAnnoList, getGenes)
str(geneL)
# List of 6
#  $ A2    : chr [1:12009] "Ga01G0007" "Ga01G0012" "Ga01G0017" "Ga01G0018" ...
#  $ D5    : chr [1:10128] "Gorai.001G001000" "Gorai.001G001700" "Gorai.001G002500" "Gorai.001G002600" ...
#  $ F1.At : chr [1:8452] "Ga01G0003" "Ga01G0005" "Ga01G0005" "Ga01G0005" ...
#  $ F1.Dt : chr [1:10035] "Gorai.001G000500" "Gorai.001G001300" "Gorai.001G002400" "Gorai.001G002900" ...
#  $ AD1.At: chr [1:17351] "Gohir.A01G000100" "Gohir.A01G000300" "Gohir.A01G000300" "Gohir.A01G000500" ...
#  $ AD1.Dt: chr [1:16626] "Gohir.D01G000400" "Gohir.D01G000700" "Gohir.D01G000800" "Gohir.D01G000900" ...

## load 20362 ortholog quadruplets derived from Justin's orthoMCL results version3_081618
ogQ<-read.table("/lss/research/jfw-lab/Projects/MNase-seq/orthohomoeologQuadruplets101218.txt", sep="\t", header=TRUE)
head(ogQ)

# correspondece
og6 = data.frame(A2=ogQ$A2, D5=ogQ$D5, F1.At=ogQ$A2, F1.Dt=ogQ$D5, AD1.At=ogQ$At, AD1.Dt=ogQ$Dt)
rownames(og6)=ogQ$groupID
head(og6)
# within in each quadruplets, which contains MSFs and which doesn't
paMSF = data.frame( A2=ogQ$A2 %in% geneL$A2, D5=ogQ$D5 %in% geneL$D5, F1.At=ogQ$A2 %in% geneL$F1.At, F1.Dt=ogQ$D5 %in% geneL$F1.Dt, AD1.At=ogQ$At %in% geneL$AD1.At, AD1.Dt=ogQ$Dt %in% geneL$AD1.Dt)
rownames(paMSF) = ogQ$groupID
head(paMSF)
# convert to 6 digit string code
paMSFn = apply(apply(paMSF,2,as.numeric), 1, function(x)paste0(x,collapse=""))
# examine the pattern of p/a
ddply(as.data.frame(paMSF[,-1]),.(A2,D5,F1.At,F1.Dt,AD1.At,AD1.Dt),nrow)
table(paMSFn)
# gene expression, use RPKM, log2mean
l=load("RNAseq/Ranalysis/expression.rdata")
l # "A2"  "D5"  "F1"  "AD1"
expr = data.frame(A2=A2$rpkm[og6$A2,"log2mean"], D5=D5$rpkm[og6$D5,"log2mean"], F1.At=F1$rpkm[og6$F1.At,"log2mean"], F1.Dt=F1$rpkm[og6$F1.Dt,"log2mean"], AD1.At=AD1$rpkm[og6$AD1.At,"log2mean"],AD1.Dt=F1$rpkm[og6$AD1.Dt,"log2mean"])
rownames(expr) = ogQ$groupID
# save
info = "og6 - ortho-homoeolog gene group index, paMSF - presence of 1kb promoter MSFs, paMSFn - 6 digit code for presence of 1kb promoter MSFs, expr - log2 mean rpkm+1"
save(info, og6, paMSF, paMSFn, expr, file="RNAseq/Ranalysis/expNuc.rdata")

## heatmap
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
# simple heatmap of total expression
select=1:500
mat = as.matrix(expr[select,])
mat_scaled = t(apply(mat, 1, scale))
heatmap.2(mat, scale = "none", Colv=FALSE, dendrogram="row",col = bluered(100), trace = "none", density.info = "none", margins=c(8,5))
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(mat, scale = "none", Colv=FALSE, dendrogram="row",col = hmcol, trace = "none", margins=c(8,5))
# ComplexHeatmap
ha = HeatmapAnnotation(df = data.frame(ploidy =c("diploid","diploid","F1","F1","AD1","AD1"),genome=c("A","D","A","D","A","D") ), col = list(ploidy=structure(brewer.pal(3, "Set2"), names = c("diploid","F1","AD1")),genome=c("A"="pink","D"="royalblue")))
#     present absence profiles, use default redBlue color
pa = as.matrix(apply(paMSF[select,],2,as.numeric))
hm= Heatmap(pa, name="MSF", col = c("0" = "white", "1" = "blue"),  cluster_column = FALSE, top_annotation = ha, top_annotation_height = unit(4, "mm"),column_title = "Promoter MSFs", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = T)
#     expression log2rpkm profiles, use default redBlue color
he.scaled=Heatmap(mat_scaled, name = "scaled log2rpkm", cluster_column = FALSE, top_annotation_height = unit(4, "mm"),column_title = "Scaled Expression", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = TRUE)
#     expression log2rpkm profiles, use default redBlue color
he=Heatmap(mat, name = "log2rpkm",col = hmcol, cluster_column = FALSE, top_annotation = ha,  top_annotation_height = unit(4, "mm"),column_title = "Expression", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = TRUE)
he+he.scaled+hm
pdf("expressionByMSFpa.pdf")
df = as.data.frame(table(paMSFn))
for(i in 1:nrow(df)){
    select=(paMSFn==df$paMSFn[i])
    mat = as.matrix(expr[select,])
    mat_scaled = t(apply(mat, 1, scale))
    colnames(mat_scaled)=colnames(mat)
    pa = as.matrix(apply(paMSF[select,],2,as.numeric))
    
    #     present absence profiles, use default redBlue color
    hm= Heatmap(pa, name="MSF", col = c("0" = "white", "1" = "blue"),  cluster_column = FALSE, top_annotation = ha, top_annotation_height = unit(4, "mm"),column_title = paste0(df$Freq[i]," genes"), column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = T)
    #     expression log2rpkm profiles, use default redBlue color
    he.scaled=Heatmap(mat_scaled, name = "scaled log2rpkm", cluster_column = FALSE, top_annotation = ha, top_annotation_height = unit(4, "mm"),column_title = "Scaled Expression", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = TRUE)
    #     expression log2rpkm profiles, use default redBlue color
    he=Heatmap(mat, name = "log2rpkm",col = hmcol, cluster_column = FALSE, top_annotation = ha, top_annotation_height = unit(4, "mm"),column_title = "Expression", column_title_gp = gpar(fontsize = 10), row_title_gp = gpar(fontsize = 10), show_row_names = FALSE, show_column_names = TRUE)
    draw(he+he.scaled+hm)
}
dev.off()

## Differential MSFs presence absence NOT in association with differential expression
makeTransparent = function(..., alpha=0.5) {
    if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
    alpha = floor(255*alpha)
    newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
    .makeTransparent = function(col, alpha) {
        rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
    }
    newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
    return(newColor)
}
plot(expr[,c("A2","D5")], pch=20, col =rgb(0,0,0,0.4))
points(expr[apply(paMSF[,c("A2","D5")],1,function(x)x[1]!=x[2]),c("A2","D5")], col=rgb(1,0,0,0.4),pch=20)
plot(expr[,c("F1.At","F1.Dt")], pch=20, col =rgb(0,0,0,0.4))
points(expr[apply(paMSF[,c("F1.At","F1.Dt")],1,function(x)x[1]!=x[2]),c("F1.At","F1.Dt")], col=rgb(1,0,0,0.4),pch=20)
plot(expr[,c("AD1.At","AD1.Dt")], pch=20, col =rgb(0,0,0,0.4))
points(expr[apply(paMSF[,c("AD1.At","AD1.Dt")],1,function(x)x[1]!=x[2]),c("AD1.At","AD1.Dt")], col=rgb(1,0,0,0.4),pch=20)


############################
# AD1 differential expression
At = log2(as.matrix(AD1$rpkm[as.character(og6$AD1.At),1:3])+1)
Dt = log2(as.matrix(AD1$rpkm[as.character(og6$AD1.Dt),1:3])+1)
count = cbind(At, Dt)
names(count) = paste0(c(rep("A",ncol(At)),rep("D",ncol(Dt))),gsub(".*-","",names(count)))
lfc = rowMeans(At)-rowMeans(Dt)
pval=apply(count,1,function(x){t.test(x[1:ncol(At)],x[(ncol(At)+1):ncol(count)])$p.value})
qval=p.adjust(diff$pval,"BH")
df=cbind(og6[,c("AD1.At","AD1.Dt")],data.frame(lfc,pval,qval))
head(df)

Abias=which(df$lfc>0 & df$qval<0.05 & !is.na(df$qval))
Dbias=which(df$lfc<0 & df$qval<0.05 & !is.na(df$qval))
unbias = setdiff(1:nrow(df),c(Abias,Dbias))
length(unbias) # 16042
length(Abias) # 2010
length(Dbias) # 2310

plot(lfc,-log2(qval), pch=20, col =rgb(0,0,0,0.4))
points(lfc[Abias],-log2(qval[Abias]), col=rgb(1,0,0,0.4),pch=20)
points(lfc[Dbias],-log2(qval[Dbias]), col=rgb(0,1,0,0.4),pch=20)

i=4
bw = "McD_q20_chr.size_w20_qnorm_qnorm.bw"
txdb = loadDb("txdb.TM1saski.sqlite")
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
TSS.all = getTSS(txdb, exclude.pattern="Gohir.1Z")
TSS.Abias.A =  getTSS(txdb, exclude.pattern="Gohir.1Z", include=df$AD1.A[Abias])
TSS.Abias.D =  getTSS(txdb, exclude.pattern="Gohir.1Z", include=df$AD1.D[Abias])
TSS.Dbias.A =  getTSS(txdb, exclude.pattern="Gohir.1Z", include=df$AD1.A[Dbias])
TSS.Dbias.D =  getTSS(txdb, exclude.pattern="Gohir.1Z", include=df$AD1.D[Dbias])
TSS.unbias.A =  getTSS(txdb, exclude.pattern="Gohir.1Z", include=df$AD1.A[unbias])
TSS.unbias.D =  getTSS(txdb, exclude.pattern="Gohir.1Z", include=df$AD1.D[unbias])


# extract target data over pre-defined windowns, multiple file
plotMetaPair = function(bw, windowA, windowD, overlay=TRUE, main="Meta-profile")
{
    sml = list(A = ScoreMatrix(target = bw, windows = windowA, type="bigWig",strand.aware=TRUE), Dt = ScoreMatrix(target = bw, windows = windowD, type="bigWig",strand.aware=TRUE) )
    names(sml) =c("A","D")
    sml= new("ScoreMatrixList",sml)
    if(overlay){
        plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="differential nuclease sensitivity score (tpm)",dispersion = "se", xlab="bases around TSS", main=main)
        abline(v=0,lty=2)
    }else{
        par(mfrow=c(1,length(sml)))
        plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="differential nuclease sensitivity score (tpm)",dispersion = "se", xlab="bases around TSS", main=main, overlay=FALSE)
    }
}

plotMetaPair(bw, windowA=TSS.Abias.A, windowD=TSS.Abias.D, main="A bias expression: AT > Dt")
plotMetaPair(bw, windowA=TSS.Dbias.A, windowD=TSS.Dbias.D, main="D bias expression: At < Dt")
plotMetaPair(bw, windowA=TSS.unbias.A, windowD=TSS.unbias.D, main="Equal expression: At = Dt")
# Abiased expression A>D, shows D profile higher in A than D

############################
# F1 differential expression
At = log2(as.matrix(F1$rpkm[as.character(og6$F1.At),1:3])+1)
Dt = log2(as.matrix(F1$rpkm[as.character(og6$F1.Dt),1:3])+1)
count = cbind(At, Dt)
names(count) = paste0(c(rep("A",ncol(At)),rep("D",ncol(Dt))),gsub(".*-","",names(count)))
lfc = rowMeans(At)-rowMeans(Dt)
pval=apply(count,1,function(x){t.test(x[1:ncol(At)],x[(ncol(At)+1):ncol(count)])$p.value})
qval=p.adjust(diff$pval,"BH")
df=cbind(og6[,c("F1.At","F1.Dt")],data.frame(lfc,pval,qval))
head(df)

Abias=which(df$lfc>0 & df$qval<0.05 & !is.na(df$qval))
Dbias=which(df$lfc<0 & df$qval<0.05 & !is.na(df$qval))
unbias = setdiff(1:nrow(df),c(Abias,Dbias))
length(Abias) # 159
length(Dbias) # 134
length(unbias) # 20069


plot(lfc,-log2(qval), pch=20, col =rgb(0,0,0,0.4))
points(lfc[Abias],-log2(qval[Abias]), col=rgb(1,0,0,0.4),pch=20)
points(lfc[Dbias],-log2(qval[Dbias]), col=rgb(0,1,0,0.4),pch=20)

bw = "FcD_q20_chr.size_w20_qnorm_qnorm.bw"
txdb = loadDb("txdb.F1.sqlite")
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
TSS.all = getTSS(txdb, exclude.pattern="Gorai.N")
TSS.Abias.A =  getTSS(txdb, exclude.pattern="Gorai.N", include=df$F1.A[Abias])
TSS.Abias.D =  getTSS(txdb, exclude.pattern="Gorai.N", include=df$F1.D[Abias])
TSS.Dbias.A =  getTSS(txdb, exclude.pattern="Gorai.N", include=df$F1.A[Dbias])
TSS.Dbias.D =  getTSS(txdb, exclude.pattern="Gorai.N", include=df$F1.D[Dbias])
TSS.unbias.A =  getTSS(txdb, exclude.pattern="Gorai.N", include=df$F1.A[unbias])
TSS.unbias.D =  getTSS(txdb, exclude.pattern="Gorai.N", include=df$F1.D[unbias])


# extract target data over pre-defined windowns, multiple file
plotMetaPair = function(bw, windowA, windowD, overlay=TRUE, main="Meta-profile")
{
    sml = list(A = ScoreMatrix(target = bw, windows = windowA, type="bigWig",strand.aware=TRUE), Dt = ScoreMatrix(target = bw, windows = windowD, type="bigWig",strand.aware=TRUE) )
    names(sml) =c("A","D")
    sml= new("ScoreMatrixList",sml)
    if(overlay){
        plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="differential nuclease sensitivity score (tpm)",dispersion = "se", xlab="bases around TSS", main=main)
        abline(v=0,lty=2)
    }else{
        par(mfrow=c(1,length(sml)))
        plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="differential nuclease sensitivity score (tpm)",dispersion = "se", xlab="bases around TSS", main=main, overlay=FALSE)
    }
}

plotMetaPair(bw, windowA=TSS.Abias.A, windowD=TSS.Abias.D, main="A bias expression: AT > Dt")
dev.new()
plotMetaPair(bw, windowA=TSS.Dbias.A, windowD=TSS.Dbias.D, main="D bias expression: At < Dt")
dev.new()
plotMetaPair(bw, windowA=TSS.unbias.A, windowD=TSS.unbias.D, main="Equal expression: At = Dt")
# Abiased expression A>D, shows D profile higher in A than D

############################
# Diploid differential expression
At = log2(as.matrix(A2$rpkm[as.character(og6$A2),1:3])+1)
Dt = log2(as.matrix(D5$rpkm[as.character(og6$D5),1:3])+1)
count = cbind(At, Dt)
names(count) = paste0(c(rep("A",ncol(At)),rep("D",ncol(Dt))),gsub(".*-","",names(count)))
lfc = rowMeans(At)-rowMeans(Dt)
pval=apply(count,1,function(x){t.test(x[1:ncol(At)],x[(ncol(At)+1):ncol(count)])$p.value})
qval=p.adjust(diff$pval,"BH")
df=cbind(og6[,c("A2","D5")],data.frame(lfc,pval,qval))
head(df)

Abias=which(df$lfc>0 & df$qval<0.05 & !is.na(df$qval))
Dbias=which(df$lfc<0 & df$qval<0.05 & !is.na(df$qval))
unbias = setdiff(1:nrow(df),c(Abias,Dbias))
length(Abias) # 136
length(Dbias) # 156
length(unbias) # 20070


plot(lfc,-log2(qval), pch=20, col =rgb(0,0,0,0.4))
points(lfc[Abias],-log2(qval[Abias]), col=rgb(1,0,0,0.4),pch=20)
points(lfc[Dbias],-log2(qval[Dbias]), col=rgb(0,1,0,0.4),pch=20)

bw = "ScD_q20_chr.size_w20_qnorm_qnorm.bw"
txdb = loadDb("txdb.F1.sqlite")
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
TSS.all = getTSS(txdb, exclude.pattern="Gorai.N")
TSS.Abias.A =  getTSS(txdb, exclude.pattern="Gorai.N", include=df$A2[Abias])
TSS.Abias.D =  getTSS(txdb, exclude.pattern="Gorai.N", include=df$D5[Abias])
TSS.Dbias.A =  getTSS(txdb, exclude.pattern="Gorai.N", include=df$A2[Dbias])
TSS.Dbias.D =  getTSS(txdb, exclude.pattern="Gorai.N", include=df$D5[Dbias])
TSS.unbias.A =  getTSS(txdb, exclude.pattern="Gorai.N", include=df$A2[unbias])
TSS.unbias.D =  getTSS(txdb, exclude.pattern="Gorai.N", include=df$D5[unbias])


# extract target data over pre-defined windowns, multiple file
plotMetaPair = function(bw, windowA, windowD, overlay=TRUE, main="Meta-profile")
{
    sml = list(A = ScoreMatrix(target = bw, windows = windowA, type="bigWig",strand.aware=TRUE), Dt = ScoreMatrix(target = bw, windows = windowD, type="bigWig",strand.aware=TRUE) )
    names(sml) =c("A","D")
    sml= new("ScoreMatrixList",sml)
    if(overlay){
        plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="differential nuclease sensitivity score (tpm)",dispersion = "se", xlab="bases around TSS", main=main)
        abline(v=0,lty=2)
    }else{
        par(mfrow=c(1,length(sml)))
        plotMeta(sml,  xcoords = c(-1000,1000), profile.names=names(sml), ylab="differential nuclease sensitivity score (tpm)",dispersion = "se", xlab="bases around TSS", main=main, overlay=FALSE)
    }
}

plotMetaPair(bw, windowA=TSS.Abias.A, windowD=TSS.Abias.D, main="A bias expression: AT > Dt")
dev.new()
plotMetaPair(bw, windowA=TSS.Dbias.A, windowD=TSS.Dbias.D, main="D bias expression: At < Dt")
dev.new()
plotMetaPair(bw, windowA=TSS.unbias.A, windowD=TSS.unbias.D, main="Equal expression: At = Dt")
# Abiased expression A>D, shows D profile higher in A than D
