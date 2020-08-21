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

#####################
## Assemble BigWig ##
#####################

setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/")
# check prepprocessed BigWig files
# --- nucleosome occupancy (NO)
bw = c(list.files("mappingA_new", pattern="*bw$",full=T),list.files("mappingD", pattern="*bw$",full=T), list.files("mappingF", pattern="*bw$",full=T),list.files("mappingMnew", pattern="*bw$",full=T))
df1 = data.frame(bigwig=bw, genome=substr(gsub(".*[/]","",bw),1,1), combo=grepl("cH|cD|cL|A6.n",bw))
df1$feature=ifelse(grepl("D_|Dn_",df1$bigwig),"DNS",ifelse(grepl("L_|Ln_",df1$bigwig),"NO_L","NO_H"))
## --- subnucleosomal particle occupancy (SPO)
bw=list.files("isegRes/SPO", pattern="*bw$",full=T)
df2 = data.frame(bigwig=bw, genome=substr(gsub(".*[/]","",bw),1,1), combo=grepl("cH|cD|cL|A6.n",bw))
df2$feature=ifelse(grepl("w21s5",df2$bigwig),"SPOfcenter","SPO")
## --- nucleosome positioning (NP)
bw=list.files("nucleosome", pattern="*bw$",full=T)
df3 = data.frame(bigwig=bw, genome=substr(gsub(".*[/]","",bw),1,1), combo=grepl("cH|cD|cL|A6.n",bw), feature="NP")
## --- ATAC
bw=list.files("ATAC/peaks", pattern="*bw$",full=T)
df4=data.frame(bigwig=bw, genome="D", combo=TRUE, feature="ATAC")
## assemble
df=rbind(df1,df2,df3,df4)

###########################
## Deeptools correlation ##
###########################
cfile= "Qregulation/command.correlateFeatures.sh"
cat("module load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
for(genome in c("A","D","F","M")){
    select = (df$genome==genome&df$feature!="NP")
    bws=as.character(df$bigwig[select])
    cmd=paste0("multiBigwigSummary bins -b ",paste(bws,collapse=" ")," -out results.npz\nplotCorrelation -in results.npz -o Qregulation/plotHeatmap_pearson.",genome,".pdf --corMethod pearson --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix Qregulation/plotHeatmap_pearson.",genome,".txt\nplotCorrelation -in results.npz -o Qregulation/plotHeatmap_spearman.",genome,".pdf --corMethod spearman --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix Qregulation/plotHeatmap_spearman.",genome,".txt\nplotPCA -in results.npz -o Qregulation/plotPCA.",genome,".pdf\n")
    cat(cmd, file =cfile,sep="\n", append=TRUE)
}
system(paste0("cat ",cfile))
system(paste0("bash ",cfile))

cfile= "Qregulation/command.correlateFeatures.combo.sh"
cat("module load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
for(genome in c("A","D","F","M")){
    select = (df$genome==genome&df$feature!="NP"&df$combo==TRUE)
    labels = paste(as.character(df$feature[select]),collapse=" ")
    bws=as.character(df$bigwig[select])
    cmd=paste0("multiBigwigSummary bins -b ",paste(bws,collapse=" ")," -l ",labels," -out results.npz\nplotCorrelation -in results.npz -o Qregulation/plotHeatmap_pearson.",genome,".combo.pdf --corMethod pearson --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix Qregulation/plotHeatmap_pearson.",genome,".combo.txt\nplotCorrelation -in results.npz -o Qregulation/plotHeatmap_spearman.",genome,".combo.pdf --corMethod spearman --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix Qregulation/plotHeatmap_spearman.",genome,".combo.txt\nplotPCA -in results.npz -o Qregulation/plotPCA.",genome,".combo.pdf\n")
    cat(cmd, file =cfile,sep="\n", append=TRUE)
}
system(paste0("cat ",cfile))
system(paste0("bash ",cfile))

##########################
## Aggregation on genes ##
##########################
refBeds = c("refGenomes/A2du.gene.bed", "refGenomes/D5.gene.bed", "refGenomes/F1.gene.bed", "refGenomes/TM1saski.gene.bed")
names(refBeds) = c("A","D","F","M")
cfile= "Qregulation/command.aggregateGenes.sh"
cat("module load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
for(genome in c("A","D","F","M")){
    bed=refBeds[genome]
    select = (df$genome==genome&df$combo==TRUE)

    
    # plot H, L, D together
    use = df[select & df$feature %in% c("DNS","NO_H","NO_L"), ]
    bws=paste(use$bigwig,collapse=" ")
    labels=paste(use$feature,collapse=" ")
    cmd=paste0("computeMatrix scale-regions -S ",bws," -R ",bed," --samplesLabel ",labels," -o scaleGeneBody.gz --regionBodyLength 3000 -b 1500 -a 1500 --skipZeros\nplotHeatmap -m scaleGeneBody.gz -out Qregulation/profileGB.",genome,".DHL1.png\nplotProfile -m scaleGeneBody.gz --perGroup -out Qregulation/profileGB.",genome,".DHL2.png")
    cat(cmd, file =cfile,sep="\n", append=TRUE)
    
    # plot each profile for others
    use = df[select & !df$feature %in% c("DNS","NO_H","NO_L"), ]
    bws=paste(use$bigwig,collapse=" ")
    labels=paste(use$feature,collapse=" ")
    cmd=paste0("computeMatrix scale-regions -S ",paste(bws)," -R ",bed," --samplesLabel ",labels," -o scaleGeneBody.gz --regionBodyLength 3000 -b 1500 -a 1500 --skipZeros\nplotProfile -m scaleGeneBody.gz --numPlotsPerRow 2 -out Qregulation/profileGB.",genome,".each.png\n")
    cat(cmd, file =cfile,sep="\n", append=TRUE)
}
system(paste0("cat ",cfile))
system(paste0("bash ",cfile))

########################
## Aggregation on TEs ##
########################
------------

## BAM to BigWig

# nucleosomal occupancy 1-4
cfile= "aggregation/command.bam2bw.sh"
cat("## BAM to Bigwig\nmodule load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
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

# subnucleosomal occupancy 5
cfile= "aggregation/command.bam2bw.moa.sh"
cat("## BAM to Bigwig\nmodule load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
dirs = c("mappingA_new","mappingD","mappingF","mappingMnew")
# effSize = c()
for(dir in dirs)
{
    print(dir)
    bams=list.files(dir, pattern="L_q20.bam$|Ln_q20.bam$",full=T)
    for(bam in bams){
        # bam to bigwig
        out = paste0("aggregation/", gsub(".*/|.bam","",bam))
        # 5. small light frenters
        cmd = paste0("bamCoverage -b ", bam," -o ",out,".moa.bw --binSize 1 --smoothLength 5 --MNase --normalizeUsingRPKM --minFragmentLength 0 --maxFragmentLength 130")
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
