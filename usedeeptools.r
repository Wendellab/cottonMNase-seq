# https://github.com/gthar/nucleServ/blob/master/bin/nucleR.R

library(data.table)
library(GenomicRanges)
library(nucleR)
library(travis)
library(rtracklayer)
thread=8



#############################################################################################
## deeptools plot nuclesome around TSS
cd nucleosome/
# AD1
computeMatrix reference-point -S Mc*bw -R ../expNuc/AD1.D*bed -o tss.gz --referencePoint TSS -b 1000 -a 1000 --skipZeros
plotHeatmap -m tss.gz -out plotHeatmap.quantile.AD1.D.png
plotProfile -m tss.gz -out plotProfileGroup.AD1.D.png
plotProfile -m tss.gz -out plotProfileGroup.AD1.D.group.png --perGroup

computeMatrix reference-point -S Mc*bw -R ../expNuc/AD1.D*bed -o tes.gz --referencePoint TES -b 1000 -a 1000 --skipZeros
plotHeatmap -m tes.gz -out plotHeatmap.quantile.AD1.D.png
plotProfile -m tes.gz -out plotProfileTES.AD1.D.png

computeMatrix reference-point -S ../mappingMnew/McD_q20_chr.size_w20_qnorm.bw -R ../expNuc/AD1.D*bed -o tss.gz --referencePoint TSS -b 1000 -a 1000 --skipZeros
plotProfile -m tss.gz -out plotProfile.dns.AD1.D.png

#############################################################################################
## check size range 0-141 vs 142-200
bedToBam -i mappingF/FcH_q20.bed -g mappingF/chr.size.txt > mappingF/FcH_q20.bam
samtools index mappingF/FcH_q20.bam
bamCoverage -b mappingF/FcH_q20.bam -o nucleosome/FcH_q20_w100k_0-200.bw --binSize 100000 --minFragmentLength 0 --maxFragmentLength 200
bamCoverage -b mappingF/FcH_q20.bam -o nucleosome/FcH_q20_w100k_142-200.bw --binSize 100000 --minFragmentLength 142 --maxFragmentLength 200
bamCoverage -b mappingF/FcH_q20.bam -o nucleosome/FcH_q20_w100k_0-141.bw --binSize 100000 --minFragmentLength 0 --maxFragmentLength 141
bigwigCompare -b1 nucleosome/FcH_q20_w100k_142-200.bw -b2 nucleosome/FcH_q20_w100k_0-200.bw --binSize 100000 -o nucleosome/FcH_q20_w100k_142-200.log2ratio.bw
bigwigCompare -b1 nucleosome/FcH_q20_w100k_0-141.bw -b2 nucleosome/FcH_q20_w100k_0-200.bw --binSize 100000 -o nucleosome/FcH_q20_w100k_0-141.log2ratio.bw
#start r
library(rtracklayer)
grL<- import.bw("~/Downloads/FcH_q20_w100k_142-200.log2ratio.bw")
grS<- import.bw("~/Downloads/FcH_q20_w100k_0-141.log2ratio.bw")
use=which(seqnames(grL)=="A01")
plot(end(grL[use]),grL[use]$score, xlab = "A01", col = 'blue', type="l", ylim=c(-2,0))
lines(end(grS[use]),grS[use]$score,type="l", col = 'red')


bamCoverage -b F2H_q20.bam -o F2H_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 1878118078
computeMatrix scale-regions -S *-*bw -R ../expNuc/F1.A.Q*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out plotSizediff.A.png
computeMatrix scale-regions -S *-*bw -R ../expNuc/F1.D.Q*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out plotSizediff.D.png
# F2L 142-200 77.7%, 0-141 21.2%, L should have lower portion of small fragment than H
# F2H 142-200 68.1%, 0-141 30.6%
~/kent/bigWigToBedGraph F2H_q20_0-141.bw F2H_q20_0-141.bg
head F2H_q20_0-141.bg


cd mappingMnew
bamCoverage -b M1L_q20.bam -o M1L_q20_142-200.bw --binSize 20 --minFragmentLength 142 --maxFragmentLength 200 --normalizeTo1x 1840638571
bamCoverage -b M1L_q20.bam -o M1L_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 1840638571
bamCoverage -b M1H_q20.bam -o M1H_q20_142-200.bw --binSize 20 --minFragmentLength 142 --maxFragmentLength 200 --normalizeTo1x 1840638571
bamCoverage -b M1H_q20.bam -o M1H_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 1840638571
computeMatrix scale-regions -S *-*bw -R ../expNuc/AD1.A.Q*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out plotSizediff.A.png
computeMatrix scale-regions -S *-*bw -R ../expNuc/AD1.D.Q*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out plotSizediff.D.png
# M1L 142-200 72.2%, 0-141 26.8%, L should have lower portion of small fragment than H
# M1H 142-200 50.4%, 0-141 49.6%

cd ../mappingA_new/
bamCoverage -b A6Ln_q20.bam -o A6Ln_q20_142-200.bw --binSize 20 --minFragmentLength 142 --maxFragmentLength 200 --normalizeTo1x 1167736467
bamCoverage -b A6Ln_q20.bam -o A6Ln_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 1167736467
bamCoverage -b A6Hn_q20.bam -o A6Hn_q20_142-200.bw --binSize 20 --minFragmentLength 142 --maxFragmentLength 200 --normalizeTo1x 1167736467
bamCoverage -b A6Hn_q20.bam -o A6Hn_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 1167736467
computeMatrix scale-regions -S *-*bw -R ../expNuc/A2.Q*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out plotSizediff.png
# A6Ln 142-200 66.4%, 0-141 30.4%, L should have lower portion of small fragment than H
# A6Hn 142-200 53.6%, 0-141 44.9%

cd ../mappingD/
bamCoverage -b D1L_q20.bam -o D1L_q20_142-200.bw --binSize 20 --minFragmentLength 142 --maxFragmentLength 200 --normalizeTo1x 716593377
bamCoverage -b D1L_q20.bam -o D1L_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 716593377
bamCoverage -b D1H_q20.bam -o D1H_q20_142-200.bw --binSize 20 --minFragmentLength 142 --maxFragmentLength 200 --normalizeTo1x 716593377
bamCoverage -b D1H_q20.bam -o D1H_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 716593377
computeMatrix scale-regions -S *-*bw -R ../expNuc/D5.Q*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out plotSizediff.png
# D1L 142-200 66.7%, 0-141 32.5%, L should have lower portion of small fragment than H
# D1H 142-200 58.7%, 0-141 40.5%
# IGV shows no apparent difference in read distribution of different size range



#############################################################################################
## use NucTools to caculate NRL
/work/LAS/jfw-lab/hugj2006/NucTools
grep '^A' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingF/FcH_q20.bed >FcH_q20.A.bed
grep '^D' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingF/FcH_q20.bed >FcH_q20.D.bed
grep '^A' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingMnew/McH_q20.bed >McH_q20.A.bed
grep '^D' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingMnew/McH_q20.bed >McH_q20.D.bed
awk -F"\t" '{print > "McH_q20.A" $1 ".bed"}' McH_q20.A.bed
awk -F"\t" '{print > "McH_q20.D" $1 ".bed"}' McH_q20.D.bed
awk -F"\t" '{print > "FcH_q20.A" $1 ".bed"}' FcH_q20.A.bed
awk -F"\t" '{print > "FcH_q20.D" $1 ".bed"}' FcH_q20.D.bed
awk -F"\t" '{print > "DcH_q20." $1 ".bed"}' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingD/DcH_q20.bed
awk -F"\t" '{print > "A6Hn_q20." $1 ".bed"}' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingA_new/A6Hn_q20.bed
perl nucleosome_repeat_length.pl -in /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingD/DcH_q20.bed -out DcH.NRL.txt
./misc/plotNRL.R --input=DcH.NRL.txt --dir=. --sample=DcH
perl nucleosome_repeat_length.pl -in /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingA_new/A6Hn_q20.bed -out A6Hn.NRL.txt
./misc/plotNRL.R --input=A6Hn.NRL.txt --dir=. --sample=A6Hn
perl nucleosome_repeat_length.pl -in McH_q20.A.bed -out McHA.NRL.txt
./misc/plotNRL.R --input=McHA.NRL.txt --dir=. --sample=McHA
perl nucleosome_repeat_length.pl -in McH_q20.D.bed -out McHD.NRL.txt
./misc/plotNRL.R --input=McHD.NRL.txt --dir=. --sample=McHD

# R code
fileL=list.files(pattern="bed$")
for(file in fileL){
    cmd=paste0("perl nucleosome_repeat_length.pl -in ",file," -out NRL.txt; ./misc/plotNRL.R --input=NRL.txt --dir=. --sample=",gsub(".bed","",file))
    message(cmd)
    system(cmd)
}
NRL=data.frame(A2=c(198, 198, 197, 196, 198, 198, 196, 197, 127, 195, 196, 198,196, 198), D5=c(194,198,194,197,196,194,197,198,199,194,194,196,198, 194), F1.A=c(199,198,198,198,199,197, 197,197,150, 198,197,197, 197,199), F1.D=c(198,195,197,197,193,194,200,195,198,198,148,196,196,198), AD1.A=c(201, 199,201,197,201,198,199,199,196,200, 201,197,200,201),AD1.D=c(201,196,202,198,201,152, 200,198,193,202,130,201,201 ,201)); NRL
apply(NRL,2,mean)
sem<-sd(x)/sqrt(length(x))
apply(NRL,2,sem)
apply(NRL,2,sd)
# why some chrs are so packed??
# need to modifly "plotNRL.R" to re do all chrs, get SEs
# histone proteins
annot = read.table("~/Dropbox/MNase-seq/cottonMNase-seq_git/cotton.v2.1.mbh.arabi.txt",sep="\t", header=TRUE)
library(gdata)
histone = read.xls("~/Dropbox/MNase-seq/cottonMNase-seq_git/Gorai.key.histone.xlsx")
core = histone[grep("PF00125", histone$PFAM),]
id=as.character(core$locusName)
t.test(rowMeans(A2$rpkm[id,]), rowMeans(D5$rpkm[id,]) )
