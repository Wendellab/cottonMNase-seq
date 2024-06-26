 ########################################################
################ Regulation question 2 #################
########################################################
# Do different approaches - DNS, SPO, and ATAC provide a conserved view of chromatin assessiblity?

# module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core
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

###########################
## Deeptools correlation ##
###########################
# Qregulation1 already generate correlations between ATAC, SPO, SPfrenter, and DNS


###################
## Assemble ACRs ##
###################
library(rtracklayer)
library(GenomicAlignments)
library(GenomicFeatures)
library(plot.matrix)
library(gplots)

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")

##------ Han et al. PNAS 2023: DNase
peaks0 = import("DNase/D5/genrich.narrowPeak", format="narrowPeak")
peaks0 #59763
DHS.Han2023=peaks0
peaks0 = import("DNase/D5/MvG.intersect.bed", format="bed",extraCols = extraCols_narrowPeak)
peaks0=peaks0[width(peaks0)>5,]
peaks0 #18083
DHS.Han2023.MivG=peaks0
peaks0 = import("DNase/D5/MvG.combine.bed", format="bed")
peaks0=peaks0[width(peaks0)>5,]
peaks0 #51943
DHS.Han2023.McvG=peaks0
peaks0 = import("DNase/D5/macs2.intersect.bed", format="bed",extraCols = extraCols_narrowPeak)
peaks0=peaks0[width(peaks0)>5,]
peaks0 #27697
DHS.Han2023.Mi=peaks0
peaks0 = import("DNase/D5/macs2.combine.bed", format="bed")
peaks0 #104391
DHS.Han2023.Mc=peaks0


####------ Wang et al. NP 2018: DNase
# DNase-seq DHS macs2, combine
peaks0 = import("DNase/D5w/genrich.narrowPeak", format="narrowPeak")
peaks0 #2059
DHS.Wang2018=peaks0
#peaks0 = import("Qregulation1/D5.dnase.macs2.combine.bed", format="bed")
#peaks0 #38108
peaks0 = import("DNase/D5w/macs2.combine.bed", format="bed")
peaks0 #40488
DHS.Wang2018.Mc=peaks0
# DNase-seq DHS macs2, intersect
peaks0 = import("DNase/D5w/macs2.intersect.bed", format="bed",extraCols = extraCols_narrowPeak)
peaks0=peaks0[width(peaks0)>5,]
peaks0 #13664
#peaks0 = import("Qregulation1/D5.dnase.macs2.intersect.bed", format="bed",extraCols = extraCols_narrowPeak)
#peaks0=peaks0[width(peaks0)>5,]
#peaks0 #12931
DHS.Wang2018.Mi=peaks0

####------ You et al. The Crop Journal 2022: ATAC
peaks0 = import("ATAC/You/genrich.narrowPeak", format="narrowPeak")
peaks0 #59763
ATAC.You2022=peaks0

## This study by Udall: ATAC
# ACR homer. combined
peaks0 = import("Qregulation1/homer.combine.bed", format="bed")
peaks0 #17102
ATAC.homer.c=peaks0
# ACR homer.intersect
peaks0 = import("Qregulation1/homer.intersect.bed", format="bed")
peaks0 #7479
ATAC.homer.i=peaks0
# ACR macs2, combine
peaks0 = import("Qregulation1/macs2.combine.bed", format="bed")
peaks0 #16195
ATAC.macs2.c=peaks0
# ACR macs2, intersect
peaks0 = import("Qregulation1/macs2.intersect.bed", format="bed",extraCols = extraCols_narrowPeak)
peaks0=peaks0[width(peaks0)>5,]
peaks0 #6167
ATAC.macs2.i=peaks0
# ACRs genrich
peaks0 = import("Qregulation1/genrich.narrowPeak", format="narrowPeak")
peaks0 #11687
ATAC.genrich=peaks0

# MSFs
bl =GRanges("Chr01", IRanges(23162000, 23820000))
# blacklist region
peaks0 = import("Qregulation1/DcD_bc4.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]
h=findOverlaps(peaks,bl); length(h) # 907 in black region
MSF.bc4.0=peaks[-queryHits(h)]; length(MSF.bc4.0)# 147291

peaks0 = import("Qregulation1/DcD_bc4.5.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]; peaks
h=findOverlaps(peaks,bl); length(h)
MSF.bc4.5=peaks[-queryHits(h)]; length(MSF.bc4.5)# 114730

peaks0 = import("Qregulation1/DcD_bc5.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]; peaks
h=findOverlaps(peaks,bl); length(h)
MSF.bc5.0=peaks[-queryHits(h)]; length(MSF.bc5.0)# 93820

peaks0 = import("Qregulation1/DcD_bc5.5.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]; peaks
h=findOverlaps(peaks,bl); length(h)
MSF.bc5.5=peaks[-queryHits(h)]; length(MSF.bc5.5)# 73930

peaks0 = import("Qregulation1/DcD_bc6.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]; peaks # 63282
h=findOverlaps(peaks,bl); length(h) # 1241 in black region
MSF.bc6.0=peaks[-queryHits(h)];  length(MSF.bc6.0)# 62041

peaks0 = import("Qregulation1/DcD_bc6.5.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]; peaks # 48287
h=findOverlaps(peaks,bl); length(h) # 1238 in black region
MSF.bc6.5=peaks[-queryHits(h)];  length(MSF.bc6.5)# 47049

peaks = import("Qregulation1/Dc_peaks.narrowPeak", format="narrowPeak") #q=0.05
peaks #48502
h=findOverlaps(peaks,bl); length(h) # 555 in black region
MSF.macs2=peaks[-queryHits(h)]; MSF.macs2 #47947

peaks = import("Qregulation1/Dc1_peaks.narrowPeak", format="narrowPeak") #q=0.01
peaks #12185
h=findOverlaps(peaks,bl); length(h) # 728 in black region
MSF.macs2.q1=peaks[-queryHits(h)]; MSF.macs2.q1 # 11457

# SPO Fr
peaks0 = import("Qregulation1/D_frenters_bc4.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"] # 290381
SPO.bc4.0=peaks
peaks0 = import("Qregulation1/D_frenters_bc4.5.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]
SPO.bc4.5=peaks # 231420
peaks0 = import("Qregulation1/D_frenters_bc5.0.Fus.bed",format="bed")
peaks = peaks0[peaks0$itemRgb=="#1414FF"]
SPO.bc5.0=peaks # 140268

#####################
## ACRs comparison ##
#####################

# compare characteristics: bp, %bp, peak numer, peak width, annotation with respect to genes and TFs
res=data.frame(type=c("MSF.bc5.0","MSF.bc5.5","MSF.bc6.0","MSF.bc6.5","MSF.macs2","MSF.macs2.q1","SPO.bc4.0","SPO.bc4.5","SPO.bc5.0","ATAC.homer.c","ATAC.homer.i","ATAC.macs2.c","ATAC.macs2.i","ATAC.genrich","ATAC.You2022","DHS.Wang2018","DHS.Wang2018.Mc","DHS.Wang2018.Mi","DHS.Han2023","DHS.Han2023.Mc","DHS.Han2023.Mi","DHS.Han2023.McvG","DHS.Han2023.MivG"), stringsAsFactors =FALSE)
for(i in 1:nrow(res))  # loop
{
    peaks = get(res$type[i])
    res$number[i] = length(peaks)
    res$bp[i] = sum(width(peaks))
    res$bpPerc[i] = paste0(round(sum(width(peaks))/749228090*100,2),"%")
    res$min.width[i] = min(width(peaks))
    res$mean.width[i] = mean(width(peaks))
    res$median.width[i] = median(width(peaks))
    res$max.width[i] = max(width(peaks))
}
res
write.table(res,file="Qregulation1/compareACRs2024.txt", sep="\t",row.names=F)
resC=res

# Pairwise overlaps - percentage of query found overlaped with subject
resM = matrix(nrow=nrow(res),ncol=nrow(res))
type=res$type
rownames(resM)=type
colnames(resM)=type
for(r in 1:nrow(res))  # row as query
for(c in 1:nrow(res))  # col as subject
{
    if(r!=c) # leave diagnal
    {
        q=get(type[r])
        s=get(type[c])
        o=countOverlaps(q,s)
        resM[r,c] = length(which(o>0))/length(q)
    }
}
resM
write.table(resM,file="Qregulation1/compareACRs_overlap2024.txt", sep="\t",row.names=F)

library(corrplot)
pdf("Qregulation1/compareACRs_overlap2024.pdf")
corrplot(resM, is.corr = FALSE, col.lim = c(0, 1),addCoef.col = 'brown',number.cex=0.5,tl.col = 'black',diag = FALSE)
use= c("MSF.bc6.0","SPO.bc5.0" ,"ATAC.genrich","ATAC.You2022","DHS.Wang2018","DHS.Han2023")
resM1<-resM[use,use]
corrplot(resM1, is.corr = FALSE, col.lim = c(0, 1),addCoef.col = 'grey50',number.cex=1,tl.col = 'black',diag = FALSE)
dev.off()


####################
## ACR annotation ##
####################
library(rtracklayer)
library(GenomicAlignments)
library(GenomicFeatures)
library(dplyr)

library(plot.matrix)
library(gplots)
library(ggplot2)
library(ggridges)
library(ggthemes)
library(gridExtra)
library(viridis)
library(UpSetR)
library(RColorBrewer)

source("FUN.acr.r") # ann

## input D genome annotation
txdb <- loadDb("refGenomes/txdb.D5.sqlite")
gns= genes(txdb)
gns= gns[grep("Chr",seqnames(gns))]# 37223
# export gene and 2kb flanking region
flank = gns
start(flank) = start(gns) - 2000
end(flank) = end(gns) + 2000
start(flank)[start(flank)<1] =1
export(flank, file("temp/gene_with_2k_flanking.bed"), format ="bed")
# TE annotation
tes = import("refGenomes/G.raimondii_JGI_221_v2.0.assembly.fasta.EDTA.TEanno.gff", format="gff")
seqlevels(tes)=gsub("D5_","",seqlevels(tes))
tes=tes[grep("^Chr",seqnames(tes))]
table(tes$type)
#LTR/Gypsy      DNA/DTC DNA/Helitron      DNA/DTM      DNA/DTA  LTR/unknown
#    90406        31226        26634        85649        22295        65079
#LTR/Copia      DNA/DTT      DNA/DTH     MITE/DTA     MITE/DTM     MITE/DTH
#    35421        11727        15172        11463         7784         1915
# MITE/DTT     MITE/DTC
#       57          157

# annotation and plots
resC
for(i in 1:nrow(resC))  # loop
{
    print(flag<-resC$type[i])
    peaks = get(flag)
    peaks = peaks[width(peaks)>1,]
        
    # annotate ACRs first with TE
    peaksTE=annotateACRs_TE(peaks,tes,upsetplot=TRUE)
    p.te.upset = peaksTE$plot
    # then with respect to nearby genes
    df = annotateACRs(peaksTE$gr,gns, distance=2000)
    d_df = df[df$type=="Distal",c(1:3)]
    p_df = df[df$type=="Proximal",c(1:3)]
    g_df = df[df$type=="Genic",c(1:3)]
    write.table(d_df, file="temp/dACRs.bed",sep="\t",row.names=FALSE,col.names=FALSE, quote=FALSE)
    write.table(g_df, file="temp/gACRs.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    write.table(p_df, file="temp/pACRs.bed",sep="\t",row.names=FALSE,col.names=FALSE,quote=FALSE)
    # df3k = annotateACRs(peaks,gns, distance=3000)

    # to test whether ACR percentage in TEs is higher or lower than expected, generate null distribution by shuffling peak resgions 100 times
    res=peaksTE$enrich
    gc_permutation =data.frame(TE_family=res$TE_family)
    rownames(gc_permutation) =res$TE_family
    for(pi in 1:100){
        # shuffle peaks with BEDtools
        print(pi)
        system("bedtools shuffle -i temp/dACRs.bed -excl temp/gene_with_2k_flanking.bed -g mappingD/chr.size.txt >temp/shuffle.bed; bedtools shuffle -i temp/pACRs.bed -incl temp/gene_with_2k_flanking.bed -g mappingD/chr.size.txt >>temp/shuffle.bed; bedtools shuffle -i temp/gACRs.bed -incl temp/gene_with_2k_flanking.bed -g mappingD/chr.size.txt >>temp/shuffle.bed")
        # import peaks and get TE overlap
        pp = import("temp/shuffle.bed", format="bed")
        gc_permutation[,pi] = annotateACRs_TE(pp,tes,upsetplot=FALSE)$en[res$TE_family,"ACR_perc"]
    }
    pp$type= c(rep("Distal",nrow(d_df)),rep("Proximal",nrow(p_df)),rep("Genic",nrow(g_df)))
    res$null_perc_mean = rowMeans(gc_permutation)
    res$enrichment = log2(res$ACR_perc/res$null_perc_mean)
    empiricalP = function(x,y){sum(x>y)/length(y)}
    res$Probablity = 0
    for(i in 1:nrow(res)){
        res$Probablity[i] = empiricalP(res$ACR_perc[i],gc_permutation[i,])
    }
    res_TE =res
    # plot heatmap
    #labeledHeatmap(xx,xLabels=1:2,yLabels=res$TE_family,colors = brewer.pal(n=10,name="PRGn"))
    res0=res
    res0$enrichment[res$enrichment==-Inf] = min(res0$enrichment[!is.infinite(res0$enrichment)])
    p.te = ggplot(data = res0, aes(x = 1, y = TE_family)) + geom_tile(aes(fill = enrichment),color="white") + scale_fill_gradient2(low="purple", high="darkgreen", guide="colorbar") +theme_bw()
    ## ACRs are depleted from TEs in general

    # categorization: ACR size
    res_Gene = aggregate(width(peaks),by=list(df$type),sum)
    names(res_Gene)=c("type","Bp")
    #  Group.1       x
    #   Distal 1108508
    #    Genic 1330851
    # Proximal 1059869
    
    #### pie plot of ACR numbers by category
    # pie(as.numeric(table(df$type) , labels = c("dACRs","gACRs","pACRs"), border="white", col = c("#999999", "#E69F00", "#56B4E9"))
    Prop = as.data.frame(table(df$type))
    names(Prop)=c("type","number")
    Prop$perc = round(as.numeric(Prop$number)/nrow(df)*100,1)
    mycols=few_pal()(3)
    # Add label position
    Prop <- Prop %>%
      arrange(desc(type)) %>%
      mutate(lab.ypos = cumsum(perc) - 0.5*perc)
    Prop$label = paste0(Prop$type,"\n", Prop$perc,"%")
    # pie & donut
    p.pie = ggplot(Prop, aes(x = "", y = perc, fill = type)) +
      geom_bar(stat = "identity", color = "white") +
      coord_polar("y", start = 0) +
      geom_text(aes(y = lab.ypos, label = label), color = "white",size=6) +
      scale_fill_manual(values = mycols) +
      theme_void()
    p.donut = ggplot(Prop, aes(x = 2, y = perc, fill = type)) +
      geom_bar(stat = "identity", color = "white") +
      coord_polar("y", start = 0) +
      geom_text(aes(y = lab.ypos, label = label), color = "white",size=6) +
      scale_fill_manual(values = mycols) +
      theme_void()+xlim(0.5, 2.5)
    p.donut.smallfont = ggplot(Prop, aes(x = 2, y = perc, fill = type)) +
      geom_bar(stat = "identity", color = "white") +
      coord_polar("y", start = 0) +
      geom_text(aes(y = lab.ypos, label = label), color = "white",size=3) +
      scale_fill_manual(values = mycols) +
      theme_void()+xlim(0.5, 2.5)
    res_Gene=merge(res_Gene,Prop, by="type")

    ### distribution of ACR lengths
    # hist(width(peaks), breaks=100, xlim =c(0,1000), freq=F)
    df$width2 = df$width
    df$width2[df$width>500] = 501
    p.width = ggplot(df, aes(x = width2))+ geom_histogram(aes(color = type, fill = type),alpha=0.4, position = "identity") + scale_fill_manual(values = mycols) +  scale_color_manual(values = mycols) + theme_classic() +theme(axis.text = element_text(size=15), legend.text = element_text(size=15))

    ### density of distance to nearest genes
    #meanweight
    density_peak=function(xx){
        x=log10(xx)
        d=density(x)
        i = which.max(d$y)
        peak = d$x[i]
        #h = hist(x, breaks=500)
        #i = which.max(h$counts)
        #peak = h$mids[i]
        return(10^peak)
    }
    md <- df[df$type!="Genic",] %>%
      group_by(type) %>%
      summarise(d.submit = density_peak(distance2nearest))
    res_Gene = merge(res_Gene, as.data.frame(md),by="type",all.x=T)
    require(scales)
    p.distance = ggplot(df, aes(x = distance2nearest))+ geom_density(aes( color=type, fill = type),alpha=0.4) + scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), labels=trans_format("log10", math_format(10^.x))) + scale_fill_manual(values = mycols[c(1,3)]) +  scale_color_manual(values = mycols[c(1,3)]) + theme_classic() + geom_vline(aes(xintercept = d.submit, color = type),data = md, linetype = "dashed")+theme(axis.text = element_text(size=15), legend.text = element_text(size=15))
    #geom_label(aes(x=c(1311), y = 1.05, label = c("1.3 Kb"))) +
    #geom_label(aes(x=c(37661), y = 1.05, label = c("37.7 Mb")))

    ### get CG distribution in peaks and control
    write.table(df[,1:7], file="temp/peaks.bed",sep="\t",row.names=FALSE,col.names=FALSE, quote=FALSE)
    system("bedtools nuc -fi refGenomes/Dgenome2_13.fasta -bed temp/peaks.bed >temp/peaksCG.bed")
    df$gc_content=read.table("temp/peaksCG.bed",sep="\t")$V9
    system("bedtools nuc -fi refGenomes/Dgenome2_13.fasta -bed temp/shuffle.bed >temp/shuffleCG.bed")
    pp$gc_content=read.table("temp/shuffleCG.bed",sep="\t")$V5
    pp$group ="control"
    df$group ="observed"
    GC = rbind(df[,c("type","gc_content","group")],as.data.frame(pp)[,c("type","gc_content","group")])
    gc_anova = anova(f<-aov(gc_content~type+group,data=GC)) #  ***
    res_GC = rbind(TukeyHSD(f)$type,TukeyHSD(f)$group)
    p.gc=ggplot(GC, aes(x=type, y=gc_content, color=type, fill=type)) + geom_violin() + geom_boxplot(width=0.1, color="grey", alpha=0.2) + theme_bw() + theme(legend.position="none",axis.text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("GC content") + scale_fill_manual(values = mycols) +  scale_color_manual(values = mycols) + facet_grid(. ~ group)

    ### plot CG distribution: peaks vs control
    pdf(paste0("Qregulation1/peakPlots2024",flag,".pdf"))
    textplot(res_Gene)
    textplot(res_TE)
    textplot(res_GC)
    print(p.pie)
    print(p.donut)
    print(p.distance)
    print(p.width)
    print(p.gc)
    print(p.te.upset)
        
    p.donut.smallfont = ggplot(Prop, aes(x = 2, y = perc, fill = type)) +
      geom_bar(stat = "identity", color = "white") +
      coord_polar("y", start = 0) +
      geom_text(aes(y = lab.ypos, label = label), color = "white",size=2.5, fontface = "bold") +
      scale_fill_manual(values = mycols) + xlim(0.5, 2.5)+
      theme_void()+ theme(legend.position = "none")
    grid.arrange( p.donut.smallfont,
    p.distance+ theme(legend.position = "none",axis.text = element_text(size=8)),
    p.width+ theme(legend.position = "none",axis.text = element_text(size=8)),
    p.gc+ theme(legend.position = "none",axis.text = element_text(size=8)),
    p.te + theme(axis.text = element_text(size=8)),
    nrow = 2, ncol =3, layout_matrix = rbind(c(1,2,5), c(3,4,5)))
    dev.off()
    save(df, pp, gc_permutation, res_TE, res_Gene, res_GC, gc_anova, file=paste0("Qregulation1/peakRes2024",flag,".rdata"))
}

resC<-read.table("Qregulation1/compareACRs2024.txt", sep="\t",header=T)
flags<-resC$type
# flags<-c("MSF.bc6.5","SPO.bc5.0","ATAC.genrich","ATAC.You2022","DHS.Wang2018","DHS.Han2023")
rm(res0)
rm(d0)
for(flag in flags){
    load(paste0("Qregulation1/peakRes2024",flag,".rdata"))
    res_Gene$sample = flag
    if(exists("res0")){res0<-rbind(res0,res_Gene)}else{res0<-res_Gene}
    df$sample = flag
    if(exists("d0")){d0=rbind(d0,df[,c("width","type","distance2nearest","sample","gc_content")])}else
    {d0<-df[,c("width","type","distance2nearest","sample","gc_content")]}
}
res0$method = gsub("[.].*","",res0$sample)
d0$method = gsub("[.].*","",d0$sample)

use=c("MSF.bc6.0","SPO.bc5.0","ATAC.genrich","ATAC.You2022","DHS.Wang2018","DHS.Han2023")
res<-res0[res0$sample %in% use,]
d  <-  d0[d0$sample %in% use,]

# overlapping with plot
theme(axis.text.x = element_text(angle = 45, vjust=0.5))

#
mycols=few_pal()(3)
p.bar= ggplot(res, aes(x = sample, y = perc, fill = type)) + labs(title="Categorization",x="", y = "percentage") + geom_bar(stat = "identity", color = "white") + geom_text(aes(y = lab.ypos, label = gsub(".*\n","",label)), color = "white",size=3) + scale_fill_manual(values = mycols) + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust=0.5))
p.w= ggplot(d, aes(x=sample, y=width, fill=method)) + labs(title="Width",x="", y = "bp")+ geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, notch=T) + theme_minimal() + ylim(0,quantile(d$width,0.99)) + scale_fill_brewer(palette="Greys")+ theme(axis.text.x = element_text(angle = 45, vjust=0.5))
p.d= ggplot(d, aes(x=sample, y=distance2nearest, fill=method)) + labs(title="Distance to nearest gene",x="", y = "bp")+ geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, notch=T) + theme_minimal() + scale_y_continuous(trans="log2")+ scale_fill_brewer(palette="Greys")+ theme(axis.text.x = element_text(angle = 45, vjust=0.5))
p.g= ggplot(d, aes(x=sample, y=gc_content, fill=method)) + labs(title="GC content",x="", y = "percentage")+ geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=1, notch=T) + theme_minimal() + scale_fill_brewer(palette="Greys")+ theme(axis.text.x = element_text(angle = 45, vjust=0.5))+ theme(axis.text.x = element_text(angle = 45, vjust=0.5))
library(grid)
blank <- grid.rect(gp=gpar(col="white"))
#
pdf("Qregulation1/comparACRs2024.pdf")
grid.arrange( p.w + theme(legend.position = "none",axis.text.x = element_blank()),
              p.g+theme(legend.position = "none",axis.text.x = element_blank()),
              p.bar+theme(legend.position = "bottom",axis.text = element_text(size=8)) ,
              p.d+theme(legend.position = "none",axis.text = element_text(size=8)),
              blank, nrow = 3, ncol =2,layout_matrix = rbind(c(1,5), c(2,4),c(3,5),c(3,5))
)
p.w
p.g
p.d
p.bar
dev.off()



############################
## Aggregate ACR on genes ##
############################

library(ChIPseeker)
library(gridExtra)
library(dplyr)
source("FUN.aggreNvisual.r")
source("chipseeker.utilities.R")
source("plotTagMatrix.r")

resC

# peak lists
resC$type

# gene windows
tss = getPoint(gns,"start",be=3000, af=1000)
tes = getPoint(gns,"end",be=1000, af=3000)

# aggregate
grl=list(MSF.bc6.0,SPO.bc5.0,ATAC.genrich,ATAC.You2022,DHS.Wang2018,DHS.Han2023)
names(grl) = c("MSF.bc6.0","SPO.bc5.0","ATAC.genrich","ATAC.You2022", "DHS.Wang2018","DHS.Han2023")
pdf("Qregulation1/ACR_over_genes2024.pdf")
tagMatrixListS = lapply(grl, getTagMatrix, windows=tss)
tagMatrixListT = lapply(grl, getTagMatrix, windows=tes)
p1=plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS")
p2=plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TES")
ylim= range(c(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range))
grid.arrange(p1+ theme(legend.position = 'bottom') +  scale_y_continuous(limits=ylim), p2+ theme(legend.position = 'bottom')+  scale_y_continuous(limits=ylim), nrow = 1)
plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS",facet="row")
plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TES",facet="row")
dev.off()

### compare ATAC-seq ACRs
grl=list(ATAC.homer.c,ATAC.macs2.c,ATAC.genrich)
names(grl) = c("ATAC.homer.c","ATAC.macs2.c","ATAC.genrich")
# aggregate
pdf("Qregulation/ACRs_D/ACR_over_genes.atac.pdf")
tagMatrixListS = lapply(grl, getTagMatrix, windows=tss)
tagMatrixListT = lapply(grl, getTagMatrix, windows=tes)
p1=plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS")
p2=plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TES")
ylim= range(c(ggplot_build(p1)$layout$panel_scales_y[[1]]$range$range,ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range))
grid.arrange(p1+ theme(legend.position = 'bottom') +  scale_y_continuous(limits=ylim), p2+ theme(legend.position = 'bottom')+  scale_y_continuous(limits=ylim), nrow = 1)
plotAvgProf.internal(tagMatrixListS, xlim=c(-3000, 1000), origin_label = "TSS",facet="row")
plotAvgProf.internal(tagMatrixListT, xlim=c(-1000, 3000), origin_label = "TES",facet="row")
dev.off()

    

#######################
## Cross Aggregation ##
#######################

# peak lists
resC$type
# ACR homer
peakl=list(MSF.bc6,MSF.bc6.5,SPO.bc4.5,SPO.bc5,ATAC.homer.c,ATAC.homer.i,ATAC.genrich,ATAC.macs2.c,ATAC.macs2.i, DHS.macs2.c,DHS.macs2.i)
names(peakl)=resC$type

# ATAC summit by macs2.i
ATAC.summit <- ATAC.macs2.i
s =start(ATAC.summit)+ATAC.summit$peak-1
end(ATAC.summit)=s
start(ATAC.summit)=s
# DHS summit by mascs2.i
DHS.summit <- DHS.macs2.i
s =start(DHS.summit)+DHS.summit$peak-1
end(DHS.summit)=s
start(DHS.summit)=s
# window list
grl=list(MSF.bc6.5,SPO.bc5,ATAC.summit, DHS.summit)
names(grl) = c("MSF.bc6.5", "SPO.bc5","ATAC.summit", "DHS.summit")

# ACRs over summit/midpoint windows
for(i in 1:length(grl))  # loop
{
    print(flag<-names(grl)[i])
    peaks = get(flag)
    # get center window
    mid=round((start(peaks)+end(peaks))/2)
    cc = peaks
    start(cc)=mid-2000
    end(cc)=mid+2000
    tagMatrixListC = lapply(peakl, getTagMatrix, windows=cc)
    
    pdf(paste0("Qregulation/ACRs_D/ACR_over_",flag,".pdf"))
    print(plotAvgProf.internal(tagMatrixListC, xlim=c(-2000, 2000), origin_label = "Center"))
    print(plotAvgProf.internal(tagMatrixListC, xlim=c(-2000, 2000), origin_label = "Center",facet="row"))
    for(j in 1:length(peakl))
    {
        tagHeatmap(tagMatrixListC[[j]], xlim=c(-2000, 2000), color="red", title=names(tagMatrixListC)[j])
    }
    dev.off()
}

# Plot D BWs over ATAC summit
export(ATAC.summit,"Qregulation/ACRs_D/ATAC.summit.bed")
export(DHS.summit,"Qregulation/ACRs_D/DHS.summit.bed")
load("Qregulation/bws.rdata")
genome="D"
select = (df$genome==genome & df$feature!="NP" & (df$combo==TRUE| grepl("ATAC",df$feature)|grepl("DNase",df$feature)))
use=df[select,]
cfile= "Qregulation/command.aggregateSummit.sh"
cat("module load py-deeptools\nmodule load bedtools2", file=cfile,sep="\n")
for(bed in c("Qregulation/ACRs_D/ATAC.summit.bed","Qregulation/ACRs_D/DHS.summit.bed")){
    flag=gsub("Qregulation/ACRs_D/|.summit.bed","",bed)
    for (each in use$feature){
        bws=use$bigwig[use$feature==each]
        labels=each
        cmd=paste0("computeMatrix reference-point --referencePoint TSS -S ",bws," -R ",bed," --samplesLabel ",labels," -o summit.gz -b 1500 -a 1500 --skipZeros\nplotHeatmap -m summit.gz --refPointLabel summit -out Qregulation/plot",flag,"Summit.",each,".",genome,".png")
        cat(cmd, file =cfile,sep="\n", append=TRUE)
    }
}
system(paste0("cat ",cfile))
#system(paste0("bash ",cfile))
sbatch(cfile, name="JOB", time = "26:00:00", N=1, n=8)
# Profiles show peak but heatmaps reveal little signal enrichment

# examine nucleosome positioning over ACRs. Null hypothesis is that ATAC is more depleted with nucleosomes
nucleosome = import("Nucleosome/DcH_q20.nucleosome.gff", format="gff")
grl=list(MSF.bc6.5,SPO.bc5,ATAC.summit, DHS.summit)
names(grl) = c("MSF.bc6.5", "SPO.bc5","ATAC.summit", "DHS.summit")
tagMatrixList=list()
pdf("Qregulation/ACRs_D/NP_over_ACRs.pdf")
for(i in 1:length(grl))  # loop
{
    print(flag<-names(grl)[i])
    peaks = get(flag)
    # get center window
    mid=round((start(peaks)+end(peaks))/2)
    cc = peaks
    start(cc)=mid
    end(cc)=mid
    export(cc,paste0("Qregulation/",flag,".mid.bed"),format="bed")
    start(cc)=mid-2000
    end(cc)=mid+2000
    tagMatrix = getTagMatrix(nucleosome, windows=cc)
    tagHeatmap(tagMatrix, xlim=c(-2000, 2000), color="red", title=flag)
    tagMatrixList[[i]] = tagMatrix
}
names(tagMatrixList)=names(grl)
print(plotAvgProf.internal(tagMatrixList, xlim=c(-2000, 2000), origin_label = "Center"))
print(plotAvgProf.internal(tagMatrixList, xlim=c(-2000, 2000), origin_label = "Center",facet="row"))
dev.off()

# plot NP over summit/midpoint windows
load("Qregulation/bws.rdata")
genome="D"
select = (df$genome==genome&df$combo==TRUE)
use = df[select & df$feature=="NP", ]
bws=use$bigwig
bedL =list.files("Qregulation",pattern="mid.bed",full=T)

cfile= "Qregulation/command.aggregateNP.sh"
cat("module load py-deeptools", file=cfile,sep="\n")
cmd=paste0("computeMatrix reference-point --referencePoint center -S ",bws," -R ",paste(bedL,collapse=" ")," -o summit.gz -b 1500 -a 1500 --skipZeros\nplotHeatmap -m summit.gz --refPointLabel summit -out Qregulation/ACRs_D/plotCenter.NP.png\nplotProfile -m summit.gz --refPointLabel summit --perGroup -out Qregulation/ACRs_D/plotCenter.NP.profile.png")
cat(cmd, file =cfile,sep="\n", append=TRUE)
system(paste0("cat ",cfile))
#system(paste0("bash ",cfile))
sbatch(cfile, name="JOB", time = "06:00:00", N=1, n=8)
# Profiles show peak but heatmaps reveal little signal enrichment


## Aggregation on TEs ##
########################


