# heatmap
# multiple input BW will be plotted side by side horizontally
# multiple input BED will be overlap in profile, and place vertically for heatmap

module load py-deeptools

#### nucleosome TSS ####
# AD1
computeMatrix reference-point -S nucleosome/McH_q20.nucleosome.coverage.bw -R refGenomes/AD1.gene.*bed -o aggregation/nuc_AD1_TSS.gz --referencePoint TSS -b 1000 -a 1000 --skipZeros
plotHeatmap -m aggregation/nuc_AD1_TSS.gz -out aggregation/nuc_AD1_TSS_Heatmap.png
plotProfile -m aggregation/nuc_AD1_TSS.gz --perGroup -out aggregation/nuc_AD1_TSS_ProfileGroup.png --outFileNameData aggregation/nuc_AD1_TSS_ProfileGroup.tab
# F1
computeMatrix reference-point -S nucleosome/FcH_q20.nucleosome.coverage.bw -R refGenomes/F1.gene.*bed -o aggregation/nuc_F1_TSS.gz --referencePoint TSS -b 1000 -a 1000 --skipZeros
plotHeatmap -m aggregation/nuc_F1_TSS.gz -out aggregation/nuc_F1_TSS_Heatmap.png
plotProfile -m aggregation/nuc_F1_TSS.gz --perGroup -out aggregation/nuc_F1_TSS_ProfileGroup.png --outFileNameData aggregation/nuc_F1_TSS_ProfileGroup.tab

# A2
# D5




# A2, D5, F1.At, F1.Dt, AD1.At, AD1.Dt each H/L/D with expression Q0-Q4

# AD1
computeMatrix scale-regions -S mappingMnew/Mc*bw -R expNuc/AD1.A*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out expNuc/plotHeatmap.quantile.AD1.A.png
computeMatrix scale-regions -S mappingMnew/Mc*bw -R expNuc/AD1.D*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out expNuc/plotHeatmap.quantile.AD1.D.png

# A2
computeMatrix scale-regions -S mappingA_new/A*n_q20*bw -R expNuc/A2.Q* -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out expNuc/plotHeatmap.quantile.A2.png

# D5
computeMatrix scale-regions -S mappingD/Dc*bw -R expNuc/D5.Q* -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out expNuc/plotHeatmap.quantile.D5.png

# F1
computeMatrix scale-regions -S mappingF/Fc*bw -R expNuc/F1.A.Q* -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out expNuc/plotHeatmap.quantile.F1.A.png
computeMatrix scale-regions -S mappingF/Fc*bw -R expNuc/F1.D.Q* -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out expNuc/plotHeatmap.quantile.F1.D.png

## A vs D, based on D normalized across genome
ls iseg/*bg
# combine normalized A2 and D5 D bg
cat <(sed 's/Chr/A/' iseg/A6Dn_q20_chr.size_w20_qnorm_qnorm.bg)  <(sed 's/Chr/D/' iseg/DcD_q20_chr.size_w20_qnorm_qnorm.bg) > iseg/ScD_q20_chr.size_w20_qnorm_qnorm.bg
# convert to bw
~/kent/bedGraphToBigWig iseg/ScD_q20_chr.size_w20_qnorm_qnorm.bg mappingF/chr.size.txt iseg/ScD_q20_chr.size_w20_qnorm_qnorm.bw
~/kent/bedGraphToBigWig iseg/FcD_q20_chr.size_w20_qnorm_qnorm.bg mappingF/chr.size.txt iseg/FcD_q20_chr.size_w20_qnorm_qnorm.bw
~/kent/bedGraphToBigWig iseg/McD_q20_chr.size_w20_qnorm_qnorm.bg mappingMnew/chr.size.txt iseg/McD_q20_chr.size_w20_qnorm_qnorm.bw
# S(A2+D5) At vs Dt
computeMatrix scale-regions -S iseg/ScD_q20_chr.size_w20_qnorm_qnorm.bw -R refGenomes/F1*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out expNuc/plotHeatmap.Diploids.AvsD.png
plotProfile -m scaleGeneBody.gz --perGroup -out expNuc/plotProfileGroup.Diploids.AvsD.png --outFileNameData expNuc/plotProfile.Diploids.AvsD.tab
# F1 & S(A2+D5) At vs Dt
computeMatrix scale-regions -S iseg/FcD_q20_chr.size_w20_qnorm_qnorm.bw -R refGenomes/F1*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out expNuc/plotHeatmap.F1.AvsD.png
plotProfile -m scaleGeneBody.gz --perGroup -out expNuc/plotProfileGroup.F1.AvsD.png --outFileNameData expNuc/plotProfile.F1.AvsD.tab
# AD1 At vs Dt, make sure BED file doesn't contain '1.1e+07'
computeMatrix scale-regions -S iseg/McD_q20_chr.size_w20_qnorm_qnorm.bw -R refGenomes/AD1*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out expNuc/plotHeatmap.AD1.AvsD.png
plotProfile -m scaleGeneBody.gz --perGroup -out expNuc/plotProfileGroup.AD1.AvsD.png --outFileNameData expNuc/plotProfile.AD1.AvsD.tab


