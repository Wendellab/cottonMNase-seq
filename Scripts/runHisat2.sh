#!usr/bin/bash

# core N
t=8

cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/
# mkdir Qregulation_diffACRs
module load hisat2 samtools py-htseq
# hisat2/2.1.0-tnpjajc
# samtools/1.10-py3-xuj7ylj
# py-htseq/0.11.2-py2-4757mqt

### build ref: hisat2-build [options]* <reference_in> <ht2_index_base>
# cd Qregulation_diffACRs
# hisat2-build ../refGenomes/TM1utx_26.fasta TM1utx_26
# hisat2-build ../refGenomes/F2020_26.fasta F2020

### nmap all to AD1 ref, count for AD1 ACRs
ls AD1acr.gff3
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
## check mapping results
ls Qregulation_diffACRs/AD1ref/*q20.txt|wc


### nmap all to F2020 (A2+D5) ref, count for combined ACRs from A2, D5 and F1
ls ADFmsf.gff3
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/
ref=Qregulation_diffACRs/F2020
#gff=Qregulation_diffACRs/ADFmsf.gff3
for j in $( ls trimmed/ | grep '_1_val_1.fq.gz' ); do
echo ''
echo $j
echo ${j%%_*}_2_val_2.fq.gz
echo ${j%%_*}.sam
echo ${j%%_*}.log
# mapping
hisat2 -x $ref -p $t -1 trimmed/$j -2 trimmed/${j%%_*}_2_val_2.fq.gz -S Qregulation_diffACRs/${j%%_*}.f.sam 2>Qregulation_diffACRs/${j%%_*}.f.log
samtools view -@ $t -b -q 20 -f 0x2 Qregulation_diffACRs/${j%%_*}.f.sam | samtools sort - -o Fref/${j%%_*}_q20.f.sort.bam
samtools index Qregulation_diffACRs/Fref/${j%%_*}_q20.f.sort.bam
#htseq-count -f bam -t MSF --idattr Name --stranded=no Qregulation_diffACRs/${j%%_*}_q20.f.bam $gff >Qregulation_diffACRs/${j%%_*}_q20.f.txt
rm Qregulation_diffACRs/${j%%_*}.f.sam
done
## check mapping results
ls *q20.f.txt|wc

### ht-seq counted reads = aligned concordantly eaxactly 1 time
#for m in $(ls *counts.txt); do echo $m; cut $m -f2 |awk '{total = total + $1}END{print total}'; grep 'aligned concordantly exactly 1 time' ${m%%_*}.log; done
