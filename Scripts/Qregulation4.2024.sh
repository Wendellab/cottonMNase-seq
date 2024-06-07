cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase

module load deeptools/3.5.1


echo 'AD1 to AD1ref'
bamCoverage -b AD1/Gh_DHS_rep1.nodup.q20.bam -o AD1/Gh_DHS_rep1.nodup.q20.rpkm.bw --binSize 20 --normalizeUsing RPKM
bamCoverage -b AD1/Gh_DHS_rep2.nodup.q20.bam -o AD1/Gh_DHS_rep2.nodup.q20.rpkm.bw --binSize 20 --normalizeUsing RPKM
bamCoverage -b AD1/Gh_DHS_rep3.nodup.q20.bam -o AD1/Gh_DHS_rep3.nodup.q20.rpkm.bw --binSize 20 --normalizeUsing RPKM

echo 'F1 to F1ref'
bamCoverage -b F1/F1_DHS_rep1.nodup.q20.bam -o F1/F1_DHS_rep1.nodup.q20.rpkm.bw --binSize 20 --normalizeUsing RPKM
bamCoverage -b F1/F2_DHS_rep1.nodup.q20.bam -o F1/F2_DHS_rep1.nodup.q20.rpkm.bw --binSize 20 --normalizeUsing RPKM
bamCoverage -b F1/F3_DHS_rep1.nodup.q20.bam -o F1/F3_DHS_rep1.nodup.q20.rpkm.bw --binSize 20 --normalizeUsing RPKM

---

echo 'F1 to AD1ref'
ls F1_AD1ref/*bw
# F1_AD1ref/F1_DHS_rep1.nodup.q20.rpkm.bw  F1_AD1ref/F1_DHS_rep3.nodup.q20.rpkm.bw F1_AD1ref/F1_DHS_rep2.nodup.q20.rpkm.bw
bamCoverage -b AD1/F1_DHS_rep1.nodup.q20.bam -o F1/F1_DHS_rep1.nodup.q20.rpkm.bw --binSize 20 --normalizeUsing RPKM
bamCoverage -b AD1/F1_DHS_rep1.nodup.q20.bam -o F1/F1_DHS_rep1.nodup.q20.rpkm.bw --binSize 20 --normalizeUsing RPKM
bamCoverage -b AD1/F1_DHS_rep1.nodup.q20.bam -o F1/F1_DHS_rep1.nodup.q20.rpkm.bw --binSize 20 --normalizeUsing RPKM

module load bedtools2

echo 'A2'
bamCoverage -b A2/Ga_DHS_rep1.nodup.q20.bam -of bedgraph -o A2/Ga_DHS_rep1.nodup.q20.rpkm.bg --binSize 20 --normalizeUsing RPKM
bamCoverage -b A2/Ga_DHS_rep2.nodup.q20.bam -of bedgraph -o A2/Ga_DHS_rep2.nodup.q20.rpkm.bg --binSize 20 --normalizeUsing RPKM
bamCoverage -b A2/Ga_DHS_rep3.nodup.q20.bam -of bedgraph -o A2/Ga_DHS_rep3.nodup.q20.rpkm.bg --binSize 20 --normalizeUsing RPKM
# Calculate average RPKM
bedtools unionbedg -i A2/Ga_DHS_rep1.nodup.q20.rpkm.bg A2/Ga_DHS_rep2.nodup.q20.rpkm.bg A2/Ga_DHS_rep3.nodup.q20.rpkm.bg >A2/Ga_DHS_union.nodup.q20.rpkm.bg
cat A2/Ga_DHS_union.nodup.q20.rpkm.bg| awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }'|sed -E 's/^Chr/A/g' >A2/Ga_DHS_mean.nodup.q20.rpkm.bg
#
echo 'D5'
bamCoverage -b D5/Gr_DHS_rep1.nodup.q20.bam -of bedgraph -o D5/Gr_DHS_rep1.nodup.q20.rpkm.bg --binSize 20 --normalizeUsing RPKM
bamCoverage -b D5/Gr_DHS_rep2.nodup.q20.bam -of bedgraph -o D5/Gr_DHS_rep2.nodup.q20.rpkm.bg --binSize 20 --normalizeUsing RPKM
bamCoverage -b D5/Gr_DHS_rep3.nodup.q20.bam -of bedgraph -o D5/Gr_DHS_rep3.nodup.q20.rpkm.bg --binSize 20 --normalizeUsing RPKM
# Calculate average RPKM
bedtools unionbedg -i D5/Gr_DHS_rep1.nodup.q20.rpkm.bg D5/Gr_DHS_rep2.nodup.q20.rpkm.bg D5/Gr_DHS_rep3.nodup.q20.rpkm.bg >D5/Gr_DHS_union.nodup.q20.rpkm.bg
cat D5/Gr_DHS_union.nodup.q20.rpkm.bg| awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }'|sed -E 's/^Chr/D/g' >D5/Gr_DHS_mean.nodup.q20.rpkm.bg
#
head A2/Ga_DHS_mean.nodup.q20.rpkm.bg
head D5/Gr_DHS_mean.nodup.q20.rpkm.bg
#
cat D5/Gr_DHS_mean.nodup.q20.rpkm.bg | awk '{sum += $4} END {print sum/NR}' #9.41894
cat A2/Ga_DHS_mean.nodup.q20.rpkm.bg | awk '{sum += $4} END {print sum/NR}' #5.03381
# Quantile to calculate (e.g., 0.5 for median)
# ../../tools/kent/bedGraphToBigWig "$avg_coverage_bedgraph" A2/chr.size.txt "$avg_coverage_bigwig"
cat D5/Gr_DHS_mean.nodup.q20.rpkm.bg | awk -v a=2 '{printf "%s\t%s\t%s\t%.6f\n", $1, $2, $3, $4/a}' | cat A2/Ga_DHS_mean.nodup.q20.rpkm.bg - >F1/A2D5_DHS_mean.nodup.q20.rpkm.bg
../../tools/kent/bedGraphToBigWig F1/A2D5_DHS_mean.nodup.q20.rpkm.bg F1/chr.size.txt F1/A2D5_DHS_mean.nodup.q20.rpkm.bw
