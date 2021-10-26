## call DNS peaks with MACS2

# Parvathaneni GB 2020: To validate results from iSeg, we also performed MACS2 [103] using the light digest as “treatment” and heavy digest as “control”, and parameters “-g 2.3e9 – nolamda –q 0.01”. Using a bc of 2.0 on DNS data, iSeg identified > 90% of the HS regions called by MACS2, in addition to many more putative regulatory regions (Additional file 1: Tables S2 and S3).
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingD

# merge reps as macs2 input
module load samtools
samtools merge DcL_q20.bam D1L_q20.bam D2L_q20.bam
samtools merge DcH_q20.bam D1H_q20.bam D2H_q20.bam

module purge
module load py-macs2
module load bedtools2
macs2 callpeak -t DcL_q20.bam -c DcH_q20.bam -f BAMPE -g 7.5e8 --nolambda -q 0.05 -n Dc --outdir macs2 2> macs2/Dc-macs2.log
macs2 callpeak -t DcL_q20.bam -c DcH_q20.bam -f BAMPE -g 7.5e8 --nolambda -q 0.01 -n Dc1 --outdir macs2 2> macs2/Dc1-macs2.log
macs2 callpeak -t D1L_q20.bam -c D1H_q20.bam -f BAMPE -g 7.5e8 --nolambda -q 0.05 -n D1 --outdir macs2 2> macs2/D1-macs2.log
macs2 callpeak -t D2L_q20.bam -c D2H_q20.bam -f BAMPE -g 7.5e8 --nolambda -q 0.05 -n D2 --outdir macs2 2> macs2/D2-macs2.log
bedtools intersect -a macs2/D1_peaks.narrowPeak -b macs2/D2_peaks.narrowPeak > macs2/reps.intersect.bed
cat macs2/D1_peaks.narrowPeak macs2/D2_peaks.narrowPeak| bedtools sort | bedtools merge > macs2/rep.combine.bed
wc macs2/*narrowPeak
#   2092   20920  145671 macs2/D1_peaks.narrowPeak
#   6884   68840  477709 macs2/D2_peaks.narrowPeak
#  12185  121850  857840 macs2/Dc1_peaks.narrowPeak
#  48502  485020 3377536 macs2/Dc_peaks.narrowPeak
#  57478  574780 4000916 total
wc macs2/*bed
#   2092   10460   92863 macs2/D1_summits.bed
#   6884   34420  307265 macs2/D2_summits.bed
#  48502  242510 2204721 macs2/Dc_summits.bed
#   7365   22095  174428 macs2/rep.combine.bed
#   1611   16110  112834 macs2/reps.intersect.bed
#  66454  325595 2892111 total
