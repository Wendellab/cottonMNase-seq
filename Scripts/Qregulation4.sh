# prepare BW for comparison

cd /work/LAS/jfw-lab/hugj2006/cottonLeaf
ls -lh isegRes/*bg
#isegRes/A6Dn_q20_chr.size_w20_qnorm_qnorm.bg  isegRes/FcD_q20_chr.size_w20_qnorm_qnorm.bg
#isegRes/DcD_q20_chr.size_w20_qnorm_qnorm.bg   isegRes/McDa_q20_chr.size_w20_qnorm_qnorm.bg
#isegRes/FcDa_q20_chr.size_w20_qnorm_qnorm.bg  isegRes/McDd_q20_chr.size_w20_qnorm_qnorm.bg
#isegRes/FcDd_q20_chr.size_w20_qnorm_qnorm.bg  isegRes/McD_q20_chr.size_w20_qnorm_qnorm.bg

# combine normalized A2 and D5 D bg
cat <(sed 's/Chr/A/' isegRes/A6Dn_q20_chr.size_w20_qnorm_qnorm.bg)  <(sed 's/Chr/D/' isegRes/DcD_q20_chr.size_w20_qnorm_qnorm.bg) > isegRes/ScD_q20_chr.size_w20_qnorm_qnorm.bg

# convert to bw
../tools/kent/bedGraphToBigWig isegRes/ScD_q20_chr.size_w20_qnorm_qnorm.bg mappingF2020/chr.size.txt Qregulation4/ScD_q20_chr.size_w20_qnorm_qnorm.bw
../tools/kent/bedGraphToBigWig isegRes/FcD_q20_chr.size_w20_qnorm_qnorm.bg mappingF2020/chr.size.txt Qregulation4/FcD_q20_chr.size_w20_qnorm_qnorm.bw
../tools/kent/bedGraphToBigWig isegRes/McD_q20_chr.size_w20_qnorm_qnorm.bg mappingM_UTX/chr.size.txt Qregulation4/McD_q20_chr.size_w20_qnorm_qnorm.bw

../tools/kent/bedGraphToBigWig isegRes/A6Dn_q20_chr.size_w20_qnorm_qnorm.bg mappingA_WHU/chr.size.txt Qregulation4/A6Dn_q20_chr.size_w20_qnorm_qnorm.bw
../tools/kent/bedGraphToBigWig isegRes/DcD_q20_chr.size_w20_qnorm_qnorm.bg mappingD/chr.size.txt Qregulation4/DcD_q20_chr.size_w20_qnorm_qnorm.bw

../tools/kent/bedGraphToBigWig isegRes/FcDa_q20_chr.size_w20_qnorm_qnorm.bg mappingA_WHU/chr.size.txt Qregulation4/FcDa_q20_chr.size_w20_qnorm_qnorm.bw
../tools/kent/bedGraphToBigWig isegRes/FcDd_q20_chr.size_w20_qnorm_qnorm.bg mappingD/chr.size.txt Qregulation4/FcDd_q20_chr.size_w20_qnorm_qnorm.bw

../tools/kent/bedGraphToBigWig isegRes/McDa_q20_chr.size_w20_qnorm_qnorm.bg mappingA_WHU/chr.size.txt Qregulation4/McDa_q20_chr.size_w20_qnorm_qnorm.bw
../tools/kent/bedGraphToBigWig isegRes/McDd_q20_chr.size_w20_qnorm_qnorm.bg mappingD/chr.size.txt Qregulation4/McDd_q20_chr.size_w20_qnorm_qnorm.bw

# locally installed deeptools 3.5.1
# BAM2BW try BPM normalization
cd /work/LAS/jfw-lab/hugj2006/tools/
git clone https://github.com/deeptools/deepTools.git
cd deepTools
python3 setup.py install --user
# Or
conda install -c bioconda deeptools
Collecting package metadata (current_repodata.json): done


#
module purge
ml py-pybigwig/0.3.4-py3-sgh7lsu
#
ml python/3.6.3-u4oaxsb py-pip/9.0.1-py3-dpds55c py-numpy/1.15.2-py3-i2pxd4u
pip3 install deeptools --user

deeptoolsintervals

ml py-deeptools

## Need A2 and D5 reads both mapped to F1 reference, normalized for direct comparison
ml samtools
# H
samtools view -b F2H_q20.f.sort.bam A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 A11 A12 A13 >F2Ha_q20.f.sort.bam
samtools view -b F2H_q20.f.sort.bam D01 D02 D03 D04 D05 D06 D07 D08 D09 D10 D11 D12 D13 >F2Hd_q20.f.sort.bam
samtools view -b F3H_q20.f.sort.bam A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 A11 A12 A13 >F3Ha_q20.f.sort.bam
samtools view -b F3H_q20.f.sort.bam D01 D02 D03 D04 D05 D06 D07 D08 D09 D10 D11 D12 D13 >F3Ha_q20.f.sort.bam
samtools view -b F2L_q20.f.sort.bam A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 A11 A12 A13 >F2La_q20.f.sort.bam
samtools view -b F2L_q20.f.sort.bam D01 D02 D03 D04 D05 D06 D07 D08 D09 D10 D11 D12 D13 >F2Ld_q20.f.sort.bam
samtools view -b F3L_q20.f.sort.bam A01 A02 A03 A04 A05 A06 A07 A08 A09 A10 A11 A12 A13 >F3La_q20.f.sort.bam
samtools view -b F3L_q20.f.sort.bam D01 D02 D03 D04 D05 D06 D07 D08 D09 D10 D11 D12 D13 >F3La_q20.f.sort.bam
#
module purge
ml py-deeptools
bamCoverage -b F2Ha_q20.f.sort.bam -o F2Ha_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 1509187826
bamCoverage -b F2Hd_q20.f.sort.bam -o F2Hd_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 749228090
bamCoverage -b F3Ha_q20.f.sort.bam -o F3Ha_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 1509187826
bamCoverage -b F3Hd_q20.f.sort.bam -o F3Hd_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 749228090
bamCoverage -b F2La_q20.f.sort.bam -o F2La_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 1509187826
bamCoverage -b F2Ld_q20.f.sort.bam -o F2Ld_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 749228090
bamCoverage -b F3La_q20.f.sort.bam -o F3La_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 1509187826
bamCoverage -b F3Ld_q20.f.sort.bam -o F3Ld_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 749228090

# L-H
bigwigCompare -b1 F2La_q20.full.bw -b2 F2Ha_q20.full.bw -o F2Da_q20.full.bw --ratio subtract
bigwigCompare -b1 F2Ld_q20.full.bw -b2 F2Hd_q20.full.bw -o F2Dd_q20.full.bw --ratio subtract
bigwigCompare -b1 F3La_q20.full.bw -b2 F3Ha_q20.full.bw -o F3Da_q20.full.bw --ratio subtract
bigwigCompare -b1 F3Ld_q20.full.bw -b2 F3Hd_q20.full.bw -o F3Dd_q20.full.bw --ratio subtract
# mean
bigwigCompare -b1 F2Da_q20.full.bw -b2 F3Da_q20.full.bw -o bwRPGC/FcDa_q20.full.bw --ratio mean
bigwigCompare -b1 F2Dd_q20.full.bw -b2 F3Dd_q20.full.bw -o bwRPGC/FcDd_q20.full.bw --ratio mean



# I want to know the difference compared to DNS-qnorm-qnorm
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/Qregulation_diffACRs/Fref
#
bamCoverage -b D1H_q20.f.sort.bam -o D1H_q20.full.bw --binSize 20 --normalizeUsingRPKM
bamCoverage -b D1L_q20.f.sort.bam -o D1L_q20.full.bw --binSize 20 --normalizeUsingRPKM
bamCoverage -b D2H_q20.f.sort.bam -o D2H_q20.full.bw --binSize 20 --normalizeUsingRPKM
bamCoverage -b D2L_q20.f.sort.bam -o D2L_q20.full.bw --binSize 20 --normalizeUsingRPKM
#
bigwigCompare -b1 D1L_q20.full.bw -b2 D1H_q20.full.bw -o D1D_q20.full.bw --ratio subtract
bigwigCompare -b1 D2L_q20.full.bw -b2 D2H_q20.full.bw -o D2D_q20.full.bw --ratio subtract
#
bigwigCompare -b1 D1D_q20.full.bw -b2 D2D_q20.full.bw -o DcD_q20.full.bw --ratio mean
#
bamCoverage -b A6Hn_q20.f.sort.bam -o A6Hn_q20.full.bw --binSize 20 --normalizeUsingRPKM
bamCoverage -b A6Ln_q20.f.sort.bam -o A6Ln_q20.full.bw --binSize 20 --normalizeUsingRPKM
#
bigwigCompare -b1 A6Ln_q20.full.bw -b2 A6Hn_q20.full.bw -o A6Dn_q20.full.bw --ratio subtract


##  Need F1 and AD1 reads both mapped to AD1 reference, Wr
###### RPKM #####
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/Qregulation_diffACRs/AD1ref/mapping
#
bamCoverage -b F2H_q20.sort.bam -o F2H_q20.full.bw --binSize 20 --normalizeUsingRPKM
bamCoverage -b F2L_q20.sort.bam -o F2L_q20.full.bw --binSize 20 --normalizeUsingRPKM
bamCoverage -b F3H_q20.sort.bam -o F3H_q20.full.bw --binSize 20 --normalizeUsingRPKM
bamCoverage -b F3L_q20.sort.bam -o F3L_q20.full.bw --binSize 20 --normalizeUsingRPKM
#
bigwigCompare -b1 F2L_q20.full.bw -b2 F2H_q20.full.bw -o F2D_q20.full.bw --ratio subtract
bigwigCompare -b1 F3L_q20.full.bw -b2 F3H_q20.full.bw -o F3D_q20.full.bw --ratio subtract
#
bigwigCompare -b1 F2D_q20.full.bw -b2 F3D_q20.full.bw -o FcD_q20.full.bw --ratio mean
#
bamCoverage -b M1H_q20.sort.bam -o M1H_q20.full.bw --binSize 20 --normalizeUsingRPKM
bamCoverage -b M1L_q20.sort.bam -o M1L_q20.full.bw --binSize 20 --normalizeUsingRPKM
bamCoverage -b M2H_q20.sort.bam -o M2H_q20.full.bw --binSize 20 --normalizeUsingRPKM
bamCoverage -b M2L_q20.sort.bam -o M2L_q20.full.bw --binSize 20 --normalizeUsingRPKM
#
bigwigCompare -b1 M1L_q20.full.bw -b2 M1H_q20.full.bw -o M1D_q20.full.bw --ratio subtract
bigwigCompare -b1 M2L_q20.full.bw -b2 M2H_q20.full.bw -o M2D_q20.full.bw --ratio subtract
#
bigwigCompare -b1 M1D_q20.full.bw -b2 M2D_q20.full.bw -o McD_q20.full.bw --ratio mean
#
mkdir bwRPKM
mv McD_q20.full.bw bwRPKM/
mv FcD_q20.full.bw bwRPKM/
rm *bw

####### RPGC #####
bamCoverage -b F2H_q20.sort.bam -o F2H_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 2258415916
bamCoverage -b F2L_q20.sort.bam -o F2L_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 2258415916
bamCoverage -b F3H_q20.sort.bam -o F3H_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 2258415916
bamCoverage -b F3L_q20.sort.bam -o F3L_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 2258415916
#
bigwigCompare -b1 F2L_q20.full.bw -b2 F2H_q20.full.bw -o F2D_q20.full.bw --ratio subtract
bigwigCompare -b1 F3L_q20.full.bw -b2 F3H_q20.full.bw -o F3D_q20.full.bw --ratio subtract
#
bigwigCompare -b1 F2D_q20.full.bw -b2 F3D_q20.full.bw -o FcD_q20.full.bw --ratio mean
#
bamCoverage -b M1H_q20.sort.bam -o M1H_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 2281618630
bamCoverage -b M1L_q20.sort.bam -o M1L_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 2281618630
bamCoverage -b M2H_q20.sort.bam -o M2H_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 2281618630
bamCoverage -b M2L_q20.sort.bam -o M2L_q20.full.bw --binSize 20 -p 4 --normalizeTo1x 2281618630
#
bigwigCompare -b1 M1L_q20.full.bw -b2 M1H_q20.full.bw -o M1D_q20.full.bw --ratio subtract
bigwigCompare -b1 M2L_q20.full.bw -b2 M2H_q20.full.bw -o M2D_q20.full.bw --ratio subtract
#
bigwigCompare -b1 M1D_q20.full.bw -b2 M2D_q20.full.bw -o McD_q20.full.bw --ratio mean
#
mkdir bwRPGC
rm *bw
