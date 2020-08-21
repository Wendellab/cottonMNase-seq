


###-------Deeptools visualzation: ATAC BW on TSS
# Convert .bam file to a normalized bigwig using the bamCoverage tool in deepTools:
# check aggregation over TSS
cat("module load py-deeptools", file="runDeeptools.sh",sep="\n")
bamL <- list.files("mapping",pattern="*q20.bam$", full=T)
for(each in bamL){
    cmd <- paste0("bamCoverage -b ",each," -o ",gsub("bam","bw",each)," -bs=1 --normalizeUsingRPKM -p=6")
    cat(cmd, file="runDeeptools.sh",sep="\n", append=TRUE)
}
cmd="computeMatrix reference-point -S mapping/*bw -R ../refGenomes/D5.gene.bed -o mapping/TSS.gz --referencePoint TSS -b 1000 -a 1000 --skipZeros; plotHeatmap -m mapping/TSS.gz -out mapping/atac_TSS_Heatmap.png"
cat(cmd, file="runDeeptools.sh",sep="\n", append=TRUE)
system("cat runDeeptools.sh")
system("bash runDeeptools.sh")

###-------Deeptools visualzation: ATAC BW on Peak centers detected by homer and macs2
#####visualize reads on detected Homer peaks
cat("module load py-deeptools", file="runDeeptools.homer.sh",sep="\n")
cmd1="computeMatrix reference-point -S mapping/*bw -R peaks/*merged.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_homer.png --colorMap=Blues --refPointLabel=center"
cmd2="computeMatrix reference-point -S mapping/*bw -R peaks/*merged.filtered.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_homer.filter.png --colorMap=Blues --refPointLabel=center"
cat(cmd1,cmd2, file="runDeeptools.homer.sh",sep="\n", append=TRUE)
system("cat runDeeptools.homer.sh")
system("bash runDeeptools.homer.sh")
#####visualize reads on detected MACS2 peaks
cat("module load py-deeptools", file="runDeeptools.macs2.sh",sep="\n")
cmd1="computeMatrix reference-point -S mapping/*bw -R peaks/*narrowPeak -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_macs2.png --colorMap=Blues --refPointLabel=center"
cmd2="computeMatrix reference-point -S mapping/*bw -R peaks/*narrowPeak.filtered.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_macs2.filter.png --colorMap=Blues --refPointLabel=center"
cat(cmd1,cmd2, file="runDeeptools.macs2.sh",sep="\n", append=TRUE)
system("cat runDeeptools.macs2.sh")
system("bash runDeeptools.macs2.sh")

###-------Deeptools visualzation: ATAC BW on MSF and subnucleosmal centers of DNS-MNase-seq
#####visualize reads on detected Homer peaks
system("grep "20,20,255" ../isegRes/isegv1.3.2_032519/DcD_bc6.0.Fus.bed  >DcD_bc6.0.Fus.MSF.bed")
system("grep "255,20,20" ../isegRes/isegv1.3.2_032519/DcD_bc6.0.Fus.bed  >DcD_bc6.0.Fus.MRF.bed")
cat("module load py-deeptools", file="runDeeptools.MSF&L0-130.sh",sep="\n")
cmd1="computeMatrix reference-point -S mapping/*bw -R ../isegRes/MOA/iseg_v1.3.2_041219/DcL_bc4.0.Fus.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_L130.png --colorMap=Blues --refPointLabel=center"
cmd2="computeMatrix reference-point -S mapping/*bw -R DcD_bc6.0.Fus.MSF.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_msf.png --colorMap=Blues --refPointLabel=center"
cmd3="computeMatrix reference-point -S mapping/*bw -R DcD_bc6.0.Fus.MRF.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_mrf.png --colorMap=Blues --refPointLabel=center"
cat(cmd1,cmd2, file="runDeeptools.MSF&L0-130.sh",sep="\n", append=TRUE)
system("cat runDeeptools.MSF\&L0-130.sh")
system("bash runDeeptools.MSF\&L0-130.sh")

###-------Deeptools visualzation: MNase-seq BW on TSS peak centers (intersect between reps)
######check DNS profiles
cat("module load py-deeptools", file="runDeeptools.DNS.sh",sep="\n")
cmd1="computeMatrix reference-point -S ../aggregation/combo/D*full.bw -R peaks/*res.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.dnsFull_atacPeaks.png --colorMap=Blues --refPointLabel=center"
cmd2="computeMatrix reference-point -S ../aggregation/combo/D*0-130*.bw -R peaks/*res.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.dnsSizeS_atacPeaks.png --colorMap=Blues --refPointLabel=center"
cmd3="computeMatrix reference-point -S ../aggregation/combo/D*130-260*.bw -R peaks/*res.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.dnsSizeL_atacPeaks.png --colorMap=Blues --refPointLabel=center"
cmd4="computeMatrix reference-point -S ../aggregation/combo/D*mnase.bw -R peaks/*res.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.mnase_atacPeaks.png --colorMap=Blues --refPointLabel=center"
cmd5="computeMatrix reference-point -S mapping/*bw -R peaks/*res.bed -o peaks/center.gz --referencePoint center -b 1000 -a 1000 --skipZeros; plotHeatmap -m peaks/center.gz -out peaks/plotHeatmap.atac_atacPeaks.png --colorMap=Blues --refPointLabel=center"
cat(cmd1,cmd2,cmd3,cmd4,cmd5, file="runDeeptools.DNS.sh",sep="\n", append=TRUE)
system("cat runDeeptools.DNS.sh")
system("bash runDeeptools.DNS.sh")

--book
### Correlaion of ATAC-seq and DNS coverages
multiBigwigSummary bins -b mapping/*bw ../aggregation/D*0-130.bw ../aggregation/D*130-260.bw ../aggregation/D*full.bw ../isegRes/MOA/D*bw -out mapping.results.npz --outRawCounts mapping.results.tab
head mapping.results.tab.txt
plotCorrelation -in mapping.results.npz -o plotHeatmap_pearson.pdf --corMethod pearson --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix plotHeatmap_pearson.tab.txt


multiBigwigSummary bins -b mapping/*bw ../isegRes/MOA/D*bw ../aggregation_RPGC/D*full.bw -out mapping.results.npz --outRawCounts mapping.results.tab.txt
head mapping.results.tab.txt
plotCorrelation -in mapping.results.npz -o plotHeatmap_pearson.pdf --corMethod pearson --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix plotHeatmap_pearson.tab.txt


################
## CONCLUSION ##
################
# 1. the quality of ATAC-seq datasets is questionable: two replicates contain 21 and 16 miliion raw PE reads (93-95% mapped), corresponding to 46% and 16% of genome region being covered by non zero reads
# 2. Homer and MACS2 peak detection were performed for both replicates, and the intersection between replicates resulted in 2.8-3.0MB peak regions, about 0.4% of the 780MB genome
# 3. iSeg peak detection failed due to overall MAD=0
# 4. Aggregation profile of DNS-seq (H, L, D, L_0-130, mnase) showed tiny bumps over ATAC peak centers; hard to tell resemablance of L_0-130 with ATAC.
#
