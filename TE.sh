# Deeptools analysis Log

# speedy.ent.iastate.edu
# /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingA_new

##### Assess reproducibility
multiBigwigSummary bins -b *q20_unified.bw -out NOC.results.npz --outRawCounts NOC.results.tab
head NOC.results.tab
plotCorrelation -in NOC.results.npz -o plotHeatmap_pearson.pdf --corMethod pearson --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix plotHeatmap_pearson.tab.txt
plotPCA -in NOC.results.npz -o plotPCA.pdf -T "PCA of NOC profiles"

##### Check sequencing coverage and fragment size
module load samtools
samtools sort A6Hn_q20.bam -o A6Hn_q20.sort.bam; samtools index A6Hn_q20.sort.bam
samtools sort A6H_q20.bam -o A6H_q20.sort.bam; samtools index A6H_q20.sort.bam
samtools sort A6Ln_q20.bam -o A6Ln_q20.sort.bam; samtools index A6Ln_q20.sort.bam
samtools sort A6L_q20.bam -o A6L_q20.sort.bam; samtools index A6L_q20.sort.bam
plotCoverage -b *sort.bam -o plotCoverage.pdf --ignoreDuplicates --minMappingQuality 20
bamPEFragmentSize -b *sort.bam -hist plotPEFragmentSize.png

##### Visualize nucleosome occupancy coverage on genomic features"
# gene body

#######
module load bedops/2.4.35-gl7y6z6
module load py-deeptools/2.5.2-py2-lgbtqfe

cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/TEannotatinon
awk '{print $3}' TEannotation/D5_TE.gff  |sort |uniq

grep "Gypsy" D5_TE.gff | sed 's/^D5_//g' |gff2bed >D5.gypsy.bed
grep "Copia" D5_TE.gff | sed 's/^D5_//g' |gff2bed >D5.copia.bed

computeMatrix scale-regions -S /work/LAS/jfw-lab/hugj2006/cottonLeaf/aggregation/combo/D*full.bw -R D5.gypsy.bed D5.copia.bed -o scaleBody.gz --regionBodyLength 3000 -b 1500 -a 1500 --skipZeros
plotHeatmap -m scaleBody.gz -out plotHeatmap_scaleBody.png

# TSS
computeMatrix reference-point -S *q20_unified.bw -R A2.gene.bed -o TSS.gz --referencePoint TSS -b 1500 -a 1500 --skipZeros
plotHeatmap -m TSS.gz -out plotHeatmap_TSS.png
# TTS
computeMatrix reference-point -S *q20_unified.bw -R A2.gene.bed -o TES.gz --referencePoint TES -b 1500 -a 1500 --skipZeros
plotHeatmap -m TES.gz -out plotHeatmap_TES.png
# get gene body tab file for R plot later
plotProfile -m scaleGeneBody.gz --perGroup -out plotProfileGroup.png --outFileNameData plotProfile.tab



mkdir dtOut
mv plotPCA.pdf dtOut/
mv plotHeatmap_pearson.pdf dtOut/
mv plotCoverage.pdf dtOut/
mv plotPEFragmentSize.png dtOut/