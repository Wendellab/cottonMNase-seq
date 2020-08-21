module load samtools
samtools sort -n mapping/raimondii_R25_ATAC.bam -o mapping/raimondii_R25_ATAC.sortN.bam
samtools sort -n mapping/raimondii_R26_ATAC.bam -o mapping/raimondii_R26_ATAC.sortN.bam

/work/LAS/jfw-lab/hugj2006/tools/Genrich/Genrich -t mapping/raimondii_R25_ATAC.sortN.bam,mapping/raimondii_R26_ATAC.sortN.bam \
  -o peaks/genrich.narrowPeak  \
  -f peaks/genrich.PQvalue.txt -k peaks/genrich.pileup.txt \
  -r  -x  -q 0.05  -a 20.0  -v  \
  -E blacklist.bed\
  >peaks/genrich.log 2>&1

 
# split pileup files to bedgrapgh each rep
sed s/chr/#chr/g peaks/genrich.pileup.txt | cut -f1-4 | csplit -z - /'# experimental file'/ {*}
# bedgrapgh to bigwig  
/work/LAS/jfw-lab/hugj2006/tools/kent/bedGraphToBigWig xx00 ../mappingD/chr.size.txt peaks/genrich.pileup.R25.bw
/work/LAS/jfw-lab/hugj2006/tools/kent/bedGraphToBigWig xx01 ../mappingD/chr.size.txt peaks/genrich.pileup.R26.bw
  
 