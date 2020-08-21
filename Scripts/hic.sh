#!/bin/bash

###########################
## 0. Set up Working Dir ##
###########################
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/HiC
mkdir D5
mkdir rawFastq/done # once trimmed move insider
mkdir trimmedFastq
mkdir trimmedFastq/done # once mapped move inside
mkdir QCreport
mkdir QCreport/raw
mkdir QCreport/trimmed
mkdir refGenome
mkdir mapping
tree

# install juicer
module load git jdk
git clone https://github.com/theaidenlab/juicer.git
ln -s juicer/SLURM/scripts/ scripts
cd scripts
wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.21.01.jar
ln -s juicer_tools_1.21.01.jar juicer_tools.jar


cd juicer/
ln -s CPU scripts
cd scripts/common
wget https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.21.01.jar
ln -s juicer_tools.1.7.6_jcuda.0.8.jar juicer_tools.jar
cd ..

ncores=13
module load wget

# D5
ref='../../refGenomes/D5'
ref_chr_size='../../mappingD/chr.size.txt'
ref_fasta='../../refGenomes/Dgenome2_13.fasta'
gene_bed='../../refGenomes/D5.gene.bed'
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/HiC/D5
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/008/SRR5885858/SRR5885858_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/008/SRR5885858/SRR5885858_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/007/SRR5885857/SRR5885857_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/007/SRR5885857/SRR5885857_2.fastq.gz

# A2
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/HiC/A2
ref='../../refGenomes/A2WHU'
ref_chr_size='../../mappingA_WHU/chr.size.txt'
ref_fasta='../../refGenomes/A2WHU_13.fasta'
gene_bed='../../refGenomes/A2WHU.gene.bed'
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/005/SRR5885855/SRR5885855_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/005/SRR5885855/SRR5885855_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/006/SRR5885856/SRR5885856_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR588/006/SRR5885856/SRR5885856_2.fastq.gz

# AD1
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/HiC/AD1
ref='../../refGenomes/AD1'
ref_chr_size='../../mappingM_UTX/chr.size.txt'
ref_fasta='../../refGenomes/TM1utx_26.fasta'
gene_bed='../../refGenomes/AD1utx.gene.bed'
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR499/008/SRR4996198/SRR4996198_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR499/008/SRR4996198/SRR4996198_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR499/000/SRR4996200/SRR4996200_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR499/000/SRR4996200/SRR4996200_2.fastq.gz

##########################
## 1. Download and trim ##
##########################

module load trimgalore
SUFFIX=_1.fastq.gz
for j in $( ls | grep _1.fastq.gz$ );    do
echo ''
echo "==========Running Trim_galore for"
echo ${j%$SUFFIX}
echo ""
time trim_galore --cores 4 --paired -o . $j ${j%$SUFFIX}_2.fastq.gz
done

### Quality Control
module load fastqc
module load py-multiqc
fastqc *fq.gz
multiqc .

#####################################################
## 2. General workflow of Hi-C analysis with HOMER ##
#####################################################
module load gcc bowtie2
module load jdk
# use my own Homer
export PATH=$PATH:/work/LAS/jfw-lab/hugj2006/tools/homer/bin/
module load r/3.5.0-py2-ufvuwmm

## From trimmed FASTQ to TagDir and .hic
for j in $( ls | grep _1_val_1.fq.gz$ );    do

echo ''
echo ${j%%_*}

# Trim (HindIII) and map sequences to the genome
echo ''
echo "======Trim (HindIII) and map sequences to the genome"
time homerTools trim -3 AAGCTAGCTT -mis 0 -matchStart 20 -min 20 $j
time homerTools trim -3 AAGCTAGCTT -mis 0 -matchStart 20 -min 20 ${j%%_*}_2_val_2.fq.gz
time bowtie2 -p $ncores -x $ref -U $j.trimmed > ${j%%_*}_1.sam 2>${j%%_*}_1.log
time bowtie2 -p $ncores -x $ref -U ${j%%_*}_2_val_2.fq.gz.trimmed > ${j%%_*}_2.sam 2>${j%%_*}_2.log

# Create Tag Directories and examine general Hi-C experiment characteristics/QC
echo ''
echo '======Create Tag Directories'
tagdir=${j%%_*}_TagDir
filteredTagdir=${j%%_*}_Fil_TagDir
makeTagDirectory $tagdir ${j%%_*}_1.sam,${j%%_*}_2.sam -tbp 1
# Further filter, easier this way to twick filtering parameters
cp -r $tagdir $filteredTagdir
makeTagDirectory $filteredTagdir -update -genome $ref_fasta -removePEbg -restrictionSite AAGCTT -both -removeSelfLigation -removeSpikes 10000 5

# Create hic for juicebox
echo ''
echo '======Create hic for juicebox'
# make short list
tagDir2hicFile.pl $filteredTagdir -p $ncores -genome $ref_chr_size -short short.txt
sort --parallel=$ncores -k2,2d -k6,6d short.txt >$filteredTagdir/short.sorted.txt
java -Xmx2g -jar ../scripts/juicer_tools.jar pre $filteredTagdir/short.sorted.txt $filteredTagdir/${j%%_*}.hic $ref_chr_size

# output contact matrix/list:
echo ''
echo '======Output contact matrix'
analyzeHiC $filteredTagdir -res 1000000 -window 2000000 -raw -balance >${j%%_*}.raw.res1MB.balance.txt  #1MB
analyzeHiC $filteredTagdir -res 10000000 -window 20000000 -balance >${j%%_*}.norm.res10MB.balance.txt  #10MB

# A/B compartment, use gene regions as seed regions to assign proper sign of PC1 result
runHiCpca.pl ${j%%_*}.pca $filteredTagdir -res 50000 -cpu $ncores -pc 1 -active $gene_bed
findHiCCompartments.pl ${j%%_*}.pca.PC1.txt >${j%%_*}.pca.A.txt  # active/permissive A compartment
findHiCCompartments.pl ${j%%_*}.pca.PC1.txt -opp >${j%%_*}.pca.B.txt # inactive/inert B compartment

# Finding TADs and Loops
echo ''
echo '======Call TADs and loops for this rep'
findTADsAndLoops.pl find $filteredTagdir -cpu $ncores

done

module load r/3.5.0-py2-ufvuwmm
echo '======Merge reps and score TADs and Loops'
# merge into a common set of features
merge2Dbed.pl *Fil_TagDir/*loop.2D.bed -loop > merged.loop.2D.bed
merge2Dbed.pl *Fil_TagDir/*tad.2D.bed -tad > merged.tad.2D.bed
# score features for comparison
findTADsAndLoops.pl score -tad merged.tad.2D.bed -loop merged.loop.2D.bed -d *Fil_TagDir -cpu 10 -o repsOut
# a wrapper for DESeq2/edgeR/limma to perform differential calculations
getDiffExpression.pl repsOut.loop.scores.txt rep1 rep2 > repsOut.diff.loop.txt
getDiffExpression.pl repsOut.tad.scores.txt rep1 rep2 > repsOut.diff.tad.txt
# check correlation between reps
getHiCcorrDiff.pl reps *Fil_TagDir

####################
# reps.corrDiff.txt - peak file containing correlation values for each region
# reps.corrDiff.bedGraph - UCSC upload file showing correlation values across the genome
# repsOut.loop.scores.txt -  merged loops with scores
# repsOut.tad.scores.txt -  merged tads with scores
# *Fil_Tagdir/*.hic - for juicebox visulization

