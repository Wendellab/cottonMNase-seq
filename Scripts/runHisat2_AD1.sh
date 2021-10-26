#!usr/bin/bash

# core N
t=8

cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/
# mkdir Qregulation_diffACRs
module load hisat2/2.1.0-tnpjajc
module load samtools/1.10-py3-xuj7ylj
#module load py-htseq/0.11.2-py2-4757mqt
module load subread/1.6.0-ak6vxhs

### build ref: hisat2-build [options]* <reference_in> <ht2_index_base>
# cd Qregulation_diffACRs
# hisat2-build ../refGenomes/TM1utx_26.fasta TM1utx_26
# hisat2-build ../refGenomes/F2020_26.fasta F2020

### nmap all to AD1 ref, count for AD1 ACRs
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/
ref=Qregulation_diffACRs/TM1utx_26
gff=Qregulation_diffACRs/AD1acr.gff3
for j in $( ls trimmed/ | grep '_1_val_1.fq.gz' ); do
echo ''
echo $j
echo ${j%%_*}_2_val_2.fq.gz
echo ${j%%_*}.sam
echo ${j%%_*}.log
# mapping
hisat2 -x $ref -p $t -1 trimmed/$j -2 trimmed/${j%%_*}_2_val_2.fq.gz -S Qregulation_diffACRs/${j%%_*}.sam 2>Qregulation_diffACRs/AD1ref/${j%%_*}.log
samtools view -@ $t -b -q 20 -f 0x2 Qregulation_diffACRs/${j%%_*}.sam | samtools sort - -o Qregulation_diffACRs/AD1ref/${j%%_*}_q20.sort.bam
samtools index Qregulation_diffACRs/AD1ref/${j%%_*}_q20.sort.bam
htseq-count -f bam -t ACR --idattr Name --stranded=no Qregulation_diffACRs/AD1ref/${j%%_*}_q20.sort.bam $gff >Qregulation_diffACRs/AD1ref/${j%%_*}_q20.txt
rm Qregulation_diffACRs/${j%%_*}.sam
done

gtf=Qregulation_diffACRs/AD1acr.gtf
featureCounts -p -t acr -g acr_id -a $gtf -o Qregulation_diffACRs/AD1ref/featureCounts.txt Qregulation_diffACRs/AD1ref/*q20.sort.bam 2>Qregulation_diffACRs/AD1ref/featureCounts.txt.log


## check mapping results
ls Qregulation_diffACRs/AD1ref/*bam
head Qregulation_diffACRs/AD1ref/*q20.txt

cd Qregulation_diffACRs/AD1ref/
# total fragments +  aligned
paste <(grep 'reads; of these:' *log|sed 's/ reads; of these://g') <(grep 'aligned concordantly exactly 1 time' *log|sed 's/ aligned concordantly exactly 1 time//g')
# featureCounts: total q20 fragments, assigned
paste <(grep 'Process BAM' featureCounts.txt.log |sed s/'\s'//g) <(grep 'Total fragments' featureCounts.txt.log|sed s/'\s'//g) <(grep 'Successfully assigned fragments' featureCounts.txt.log|sed s/'\s'//g)
# HT-seq
for m in $(ls *q20.txt); do echo $m; grep -v '^_' $m |cut -f2 |awk '{total = total + $1}END{print total}'; done # HT-seq count: total assigned counts
for m in $(ls *q20.txt); do echo $m; cut -f2 $m |awk '{total = total + $1}END{print total}'; done  # HT-seq count: total counts, including ambiguous and no feature

## Conclusion: HT-seq count is very slow and somehow incorrect
