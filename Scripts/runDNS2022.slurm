#!usr/bin/bash
echo ''
echo ''
echo "Starting Job on "
stampStart=`date`
echo $stampStart 

SUFFIX=_1.fq.gz

cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/
module load fastqc/0.11.3
module load py-multiqc
mkdir rawfastq
mv Z30/FP200007304_L01_125_1.fq.gz rawfastq/Z30_1.fq.gz
mv Z30/FP200007304_L01_125_2.fq.gz rawfastq/Z30_2.fq.gz
mv Q30/FP200007304_L01_124_1.fq.gz rawfastq/Q30_1.fq.gz
mv Q30/FP200007304_L01_124_2.fq.gz rawfastq/Q30_2.fq.gz
mkdir trimmed/
mkdir QCreport
mkdir QCreport/raw
mkdir QCreport/trimmed


ml trimgalore
SUFFIX=_1.fq.gz
for j in $( ls rawfastq/| grep _1.fq.gz$ );    do
echo ''
echo "==========Running Trim_galore for"
echo ${j%$SUFFIX}
echo ""
trim_galore --paired -o trimmed/ rawfastq/$j rawfastq/${j%$SUFFIX}_2.fq.gz
done

module purge
ml fastqc/0.11.7-3flwcvl py-multiqc/1.5-py2-lqqx3ht
fastqc -o QCreport/raw rawfastq/*gz
fastqc -o QCreport/trimmed trimmed/*gz

# results were collected from all samples into a single report for easy comparison with MultiQC (Ewels et al., 2016)
cd QCreport/raw
multiqc .
cd ../trimmed
multiqc .

echo ''
echo ''
echo "Ending  Job on "
stampEnd=`date`
echo $stampEnd
