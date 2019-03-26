#!usr/bin/bash
echo ''
echo ''
echo "Starting Job on "
stampStart=`date`
echo $stampStart 

module load bowtie2
echo "This is mapping against D ref"

for j in $( ls trimmed/ | grep 'D.*val_1' ); do

echo ''
echo ${j%%_*}
echo ${j%%_*}_2_val_2.fq.gz
echo ${j%%_*}.sam
echo ${j%%_*}.log

bowtie2 -q -p 6 -t --no-mixed --no-discordant --no-unal --dovetail -x refGenomes/D5 -1 <(zcat trimmed/$j) -2 <(zcat trimmed/${j%%_*}_2_val_2.fq.gz) -S mapping/${j%%_*}.sam 2>mapping/${j%%_*}.log

samtools view -bS mapping/${j%%_*}.sam | samtools sort - -o mapping/${j%%_*}.sort.bam ; samtools index mapping/${j%%_*}.sort.bam
# rm mapping/$j
# rm mapping/${j/SUFFIX/}.bam


done

echo ''
echo ''
echo "Ending  Job on "
stampEnd=`date`
echo $stampEnd
