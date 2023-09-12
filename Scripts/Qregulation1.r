########################################################
################ Regulation question 1 #################
########################################################
# How are various chromatin features compared with each other and in associated with gene expression?

# module load py-deeptools
# module load r/3.5.0-py2-ufvuwmm
# module load bedtools2
# cd /work/LAS/jfw-lab/hugj2006/cottonLeaf
# R
sbatch=function(cfile, name="JOB", time = "7-00:00:00", N=1, n=8, partition="whatever", mem="200G",email="hugj2006@iastate.edu"){
    slurmF=paste0(cfile,".slurm")
    cat("#!/bin/bash\n#SBATCH --constraint=AVX2\n", file=slurmF,sep="\n")
    cat(paste0("#SBATCH -t ",time,"   # walltime\n#SBATCH -N ",N,"   # number of nodes in this job\n#SBATCH -n ",n,"   # total number of processor cores in this job; each node has 272 cores, Nova has 36\n#SBATCH --mem=",mem," # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G\n#SBATCH --partition=",partition," # use sinfo to get queue names\n#SBATCH -J ",name,"   # job name\n#SBATCH --output=job%j.txt\n#SBATCH --mail-type=FAIL,END\n#SBATCH --mail-user=",email,"   # email address\n"), file=slurmF, sep="\n",append=TRUE)
    cat(paste0("## LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE\necho 'Start'\ndate\n"), file=slurmF,sep="\n",append=TRUE)
    system(paste0("cat ",slurmF," ",cfile, " >>",slurmF))
    cat(paste0("\necho 'END'\ndate\n"), file=slurmF,sep="\n",append=TRUE)
    system(paste0("cat ",slurmF))
    system(paste0("sbatch ",slurmF))
}

#####################
## Assemble BigWig ##
#####################

setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/")
# check prepprocessed BigWig files
# --- nucleosome occupancy (NO)
bw = c(list.files("mappingA_WHU", pattern="*bw$",full=T),list.files("mappingD", pattern="*bw$",full=T), list.files("mappingF2020", pattern="*bw$",full=T),list.files("mappingM_UTX", pattern="*bw$",full=T))
df1 = data.frame(bigwig=bw, genome=substr(gsub(".*[/]","",bw),1,1), combo=grepl("cH|cD|cL|A6.n",bw))
df1$feature=ifelse(grepl("D_|Dn_",df1$bigwig),"DNS",ifelse(grepl("L_|Ln_",df1$bigwig),"NO_L","NO_H"))
## --- subnucleosomal particle occupancy (SPO)
bw=list.files("isegRes/SPO", pattern="000-130.*bw$",full=T)
df2 = data.frame(bigwig=bw, genome=substr(gsub(".*[/]","",bw),1,1), combo=grepl("cH|cD|cL|A6.n",bw))
df2$feature=ifelse(grepl("w21s5",df2$bigwig),"SPOfcenter","SPO")
## --- nucleosome positioning (NP)
bw=list.files("Nucleosome", pattern="*bw$",full=T)
df3 = data.frame(bigwig=bw, genome=substr(gsub(".*[/]","",bw),1,1), combo=grepl("cH|cD|cL|A6.n",bw), feature="NP")
## --- ATAC
bw=list.files("ATAC/mapping", pattern="*filter.bw$",full=T)
df4=data.frame(bigwig=bw, genome="D", combo=FALSE, feature=c("ATAC.R25","ATAC.R26"))
## --- DNase-seq
bw=c(list.files("DNase/D5", pattern="*filter.bw$",full=T),list.files("DNase/A2", pattern="*bw$",full=T),list.files("DNase/AD1", pattern="*bw$",full=T))
df5=data.frame(bigwig=bw, genome=c("D","D","A","A","M","M","M","M"), combo=FALSE, feature=c("DNase.R1","DNase.R2","DNase.R1","DNase.R2","DNase.R1","DNase.R2","DNase.fiber.R1","DNase.fiber.R2"))

## assemble
df=rbind(df1,df2,df3,df4,df5)
save(df,file="Qregulation/bws.rdata")

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
sbatch(cfile, name="correlateFeatures", time = "06:00:00", N=1, n=8)
#system(paste0("bash ",cfile))

cfile= "Qregulation/command.correlateFeatures2.sh"
cat("module load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
for(genome in c("A","D","F","M")){
    select = (df$genome==genome&df$feature!="NP"&df$feature!="NO_H"&df$feature!="NO_L")
    bws=as.character(df$bigwig[select])
    cmd=paste0("multiBigwigSummary bins -b ",paste(bws,collapse=" ")," -out results.npz\nplotCorrelation -in results.npz -o Qregulation/corr_ACR/plotHeatmap_pearson.",genome,".pdf --corMethod pearson --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix Qregulation/corr_ACR/plotHeatmap_pearson.",genome,".txt\nplotCorrelation -in results.npz -o Qregulation/corr_ACR/plotHeatmap_spearman.",genome,".pdf --corMethod spearman --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix Qregulation/corr_ACR/plotHeatmap_spearman.",genome,".txt\nplotPCA -in results.npz -o Qregulation/corr_ACR/plotPCA.",genome,".pdf\n")
    cat(cmd, file =cfile,sep="\n", append=TRUE)
}
system(paste0("cat ",cfile))
sbatch(cfile, name="correlateFeatures", time = "06:00:00", N=1, n=8)


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
sbatch(cfile, name="correlateCombo", time = "03:00:00", N=1, n=8)
#system(paste0("bash ",cfile))

##########################
## Aggregation on genes ##
##########################
load("Qregulation/bws.rdata")
g=c("A","D","F","M")
refBeds = c("refGenomes/A2WHU.gene.bed", "refGenomes/D5.gene.bed", "refGenomes/F2020.gene.bed", "refGenomes/AD1utx.gene.bed")
refQBeds =c("refGenomes/bedQ/A2WHU.Q*", "refGenomes/bedQ/D5.Q*", "refGenomes/bedQ/F2020.Q*", "refGenomes/bedQ/AD1utx.Q*")
refQBeds_YL =c("refGenomes/bedQ_YL/A2WHU.Q*", "refGenomes/bedQ_YL/D5.Q*", "", "refGenomes/bedQ_YL/AD1utx.Q*")
names(refBeds) = g
names(refQBeds) = g
names(refQBeds_YL) = g

# Plot over gene body (-1.5, 3, 1.5)
cfile= "Qregulation/command.aggregateGenes.GB.sh"
cat("module load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
for(genome in g){
    bed=refBeds[genome]
    bedQ=refQBeds[genome]
    select = (df$genome==genome& (df$combo==TRUE |grepl("ATAC",df$feature)| grepl("DNase",df$feature)))
    use0=df[select,]
    
    for(each in use0$feature){
        use = df[select & df$feature==each, ]
        bws=use$bigwig
        labels=use$feature
        cmd=paste0("computeMatrix scale-regions -S ",bws," -R ",bedQ," --samplesLabel ",labels," -o scaleGeneBody.gz --regionBodyLength 3000 -b 1500 -a 1500 --skipZeros\nplotProfile -m scaleGeneBody.gz --regionsLabel Q0 Q1 Q2 Q3 Q4 -out Qregulation/plotGB.",each,".",genome,".png")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
    }
}
system(paste0("cat ",cfile))
#system(paste0("bash ",cfile))
sbatch(cfile, name="JOB", time = "06:00:00", N=1, n=8)

cfile= "Qregulation/command.aggregateGenes_YL.GB.sh"
cat("module load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
for(genome in c("A","D","M")){
    bed=refBeds[genome]
    bedQ=refQBeds_YL[genome]
    select = (df$genome==genome& (grepl("ATAC",df$feature)| grepl("DNase",df$feature)))
    use0=df[select,]
    
    for(each in use0$feature){
        use = df[select & df$feature==each, ]
        bws=use$bigwig
        labels=use$feature
        cmd=paste0("computeMatrix scale-regions -S ",bws," -R ",bedQ," --samplesLabel ",labels," -o scaleGeneBody.gz --regionBodyLength 3000 -b 1500 -a 1500 --skipZeros\nplotProfile -m scaleGeneBody.gz --regionsLabel Q0 Q1 Q2 Q3 Q4 -out Qregulation/plotGB.",each,".",genome,".YL.png")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
    }
}
system(paste0("cat ",cfile))
#system(paste0("bash ",cfile))
sbatch(cfile, name="JOB", time = "06:00:00", N=1, n=8)

# Plot over TSS and TES (-1.5, 1.5)
cfile= "Qregulation/command.aggregateGenes.Points.sh"
cat("module load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
for(genome in g){
    bed=refBeds[genome]
    bedQ=refQBeds[genome]
    select = (df$genome==genome&df$combo==TRUE)
    use0=df[select,]
    # c("NO_L","NO_H","DNS","SPO","SPOfcenter","NP")
    for(each in use0$feature){
        use = df[select & df$feature==each, ]
        bws=use$bigwig
        labels=use$feature
        cmd=paste0("computeMatrix reference-point --referencePoint TSS -S ",bws," -R ",bed," --samplesLabel ",labels," -o point.gz -b 1500 -a 1500 --skipZeros\nplotProfile -m point.gz -out Qregulation/plotTSS.",each,".",genome,".png")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
        cmd=paste0("computeMatrix reference-point --referencePoint TES -S ",bws," -R ",bed," --samplesLabel ",labels," -o point.gz -b 1500 -a 1500 --skipZeros\nplotProfile -m point.gz -out Qregulation/plotTES.",each,".",genome,".png")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
        cmd=paste0("computeMatrix reference-point --referencePoint TSS -S ",bws," -R ",bedQ," --samplesLabel ",labels," -o point.gz -b 1500 -a 1500 --skipZeros\nplotProfile -m point.gz  --regionsLabel Q0 Q1 Q2 Q3 Q4 -out Qregulation/plotTSS.",each,".",genome,".quantile.png")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
        cmd=paste0("computeMatrix reference-point --referencePoint TES -S ",bws," -R ",bedQ," --samplesLabel ",labels," -o point.gz -b 1500 -a 1500 --skipZeros\nplotProfile -m point.gz  --regionsLabel Q0 Q1 Q2 Q3 Q4 -out Qregulation/plotTES.",each,".",genome,".quantile.png")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
    }
}
system(paste0("cat ",cfile))
#system(paste0("bash ",cfile))
sbatch(cfile, name="JOB", time = "06:00:00", N=1, n=8)
    
------------


########################
## Aggregation on TEs ##
########################
load("Qregulation/bws.rdata")
g=c("A","D","F","M")
refBeds = c("refGenomes/A2WHU.gypsy.bed", "refGenomes/D5.gypsy.bed", "refGenomes/F2020.gypsy.bed", "refGenomes/AD1utx.gypsy.bed")
refIBeds = c("refGenomes/A2WHU.gypsy.intact.bed", "refGenomes/D5.gypsy.intact.bed", "refGenomes/F2020.gypsy.intact.bed", "refGenomes/AD1utx.gypsy.intact.bed")
names(refBeds) = g
names(refIBeds) = g

# Plot over gene body (-1.5, 3, 1.5)
cfile= "Qregulation/command.aggregateGypsy.GB.sh"
cat("module load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
for(genome in g){
    bed=refBeds[genome]
    bedI=refIBeds[genome]
    select = (df$genome==genome& (df$combo==TRUE |grepl("ATAC",df$feature)| grepl("DNase",df$feature)))
    use0=df[select,]
    
    for(each in use0$feature){
        use = df[select & df$feature==each, ]
        bws=use$bigwig
        labels=use$feature
        cmd=paste0("computeMatrix scale-regions -S ",bws," -R ",bed," ",bedI," --samplesLabel ",labels," -o scaleGeneBody.gz --regionBodyLength 3000 -b 1500 -a 1500 --skipZeros\nplotProfile -m scaleGeneBody.gz -out Qregulation/plotGypsy.",each,".",genome,".png")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
    }
}
system(paste0("cat ",cfile))
#system(paste0("bash ",cfile))
sbatch(cfile, name="JOB", time = "06:00:00", N=1, n=8)


# Plot over TSS and TES (-1.5, 1.5)
cfile= "Qregulation/command.aggregateGypsy.Points.sh"
cat("module load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
for(genome in g){
    bed=refBeds[genome]
    bedI=refIBeds[genome]
    select = (df$genome==genome&df$combo==TRUE)
    use0=df[select,]
    # c("NO_L","NO_H","DNS","SPO","SPOfcenter","NP")
    for(each in use0$feature){
        use = df[select & df$feature==each, ]
        bws=use$bigwig
        labels=use$feature
        cmd=paste0("computeMatrix reference-point --referencePoint TSS -S ",bws," -R ",bed," ",bedI," --samplesLabel ",labels," -o point.gz -b 1500 -a 1500 --skipZeros\nplotProfile -m point.gz -out Qregulation/plotGypsySS.",each,".",genome,".png")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
        cmd=paste0("computeMatrix reference-point --referencePoint TES -S ",bws," -R ",bed," ",bedI," --samplesLabel ",labels," -o point.gz -b 1500 -a 1500 --skipZeros\nplotProfile -m point.gz -out Qregulation/plotGypsyES.",each,".",genome,".png")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
    }
}
system(paste0("cat ",cfile))
#system(paste0("bash ",cfile))
sbatch(cfile, name="JOB", time = "06:00:00", N=1, n=8)
    
