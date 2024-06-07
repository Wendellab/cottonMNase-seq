cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/Qregulation_diffACRs

module load deeptools/3.5.1

cd /Fref
# At
multiBigwigSummary bins -b bwRPGC/A6Dn_q20.full.bw F2Da_q20.full.bw  F3Da_q20.full.bw  M1Da_q20.full.bw  M2Da_q20.full.bw bwRPGC/FcDa_q20.full.bw  bwRPGC/McDa_q20.full.bw -l A_r1 F1_r2 F1_r3 M_r1 M_r2 F_c M_c -out results.npz
plotCorrelation -in results.npz -o plotHeatmap_pearson.At.pdf --corMethod pearson --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix plotHeatmap_pearson.At.txt
plotCorrelation -in results.npz -o plotHeatmap_spearman.At.pdf --corMethod spearman --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix plotHeatmap_spearman.At.txt
plotPCA -in results.npz -o plotPCA.At.pdf
#Dt
multiBigwigSummary bins -b D1D_q20.full.bw  D2D_q20.full.bw F2Dd_q20.full.bw  F3Dd_q20.full.bw  M1Dd_q20.full.bw  M2Dd_q20.full.bw bwRPGC/DcD_q20.full.bw   bwRPGC/FcDd_q20.full.bw  bwRPGC/McDd_q20.full.bw -l D_r1 D_r2 F1_r2 F1_r3 M_r1 M_r2 D_c F_c M_c -out results.npz
plotCorrelation -in results.npz -o plotHeatmap_pearson.Dt.pdf --corMethod pearson --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix plotHeatmap_pearson.Dt.txt
plotCorrelation -in results.npz -o plotHeatmap_spearman.Dt.pdf --corMethod spearman --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix plotHeatmap_spearman.Dt.txt
plotPCA -in results.npz -o plotPCA.Dt.pdf


cd /AD1ref
multiBigwigSummary bins -b bwRPGC/ScD_q20.full.bw bwRPGC/FcD_q20.full.bw bwRPGC/McD_q20.full.bw -l A2D5 F1 AD1 -out results.npz
plotCorrelation -in results.npz -o plotHeatmap_pearson.pdf --corMethod pearson --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix plotHeatmap_pearson.txt
plotCorrelation -in results.npz -o plotHeatmap_spearman.pdf --corMethod spearman --skipZeros --removeOutliers --whatToPlot heatmap --colorMap RdYlBu --plotNumbers --outFileCorMatrix plotHeatmap_spearman.txt
plotPCA -in results.npz -o plotPCA.pdf
