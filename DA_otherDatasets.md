# Data analysis of other datasets 

## TE annotation

The [EDTA](https://github.com/oushujun/EDTA) pipeline was used to perform high-quality and consistentl TE annotation for three cotton genome reference genomes: **D5** ~800Mb, the smallest and highest quality assembly [lss](/work/LAS/jfw-lab/genomes/D5/Dgenome2_13.fasta) [ftp](ftp://ftp.bioinfo.wsu.edu/species/Gossypium_raimondii/JGI_221_G.raimondii_Dgenome/assembly/G.raimondii_JGI_221_v2.0.assembly.fasta.gz). **A2** ~1600Mb [lss](/lss/research/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Du2018/1.PacBio-Gar-Assembly-v1.0/G.arboreum.Chr.v1.0.fasta.gz) [ftp](ftp://ftp.bioinfo.wsu.edu/species/Gossypium_arboreum/CRI-Updated_G.arboreum_A2genome/assembly/G.arboreum_CRI-A2_assembly_v1.0.fasta.gz). **AD1** ~2200 Mb [lss](/lss/research/jfw-lab/GenomicResources/archived_resources/AD1Saski/assembly/Ghirsutum_458_v1.0.fa.gz) [ftp](ftp://ftp.bioinfo.wsu.edu/species/Gossypium_hirsutum/Tx-JGI_G.hirsutum_AD1genome/assembly/Tx-JGI_G.hirsutum_v1.1.fa.gz)

Two previous versions (v091919 and v111519) were generated, and the latest one is **v010621** which conducted TE and panTE annotation for a total of 15 available *Gossypium* reference genomes.

## Genomic annotation

Visualization and exploration of sequencing tracks (`BigWig`) around genomic features (e.g. TSS) is the key of downstream analysis. [makeTxDb.r](Scripts/makeTxDb.r) prepared `txdb.<refGenome>.sqlite` for handling gene anonation, prepared `TEannotation.rdata` for TE annotation, and generated corresponding BED files. Alternatively, BEDOPS tool [gff2bed](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gff2bed.html) can be used to convert gff to bed for deeptools, e.g. [makeGenicBed.sh](Scripts/makeGenicBed.sh) .





## ATAC-seq, D5 only

Sep 27, 2018 - Josh transferred two ATAC-seq datasets (*G.raimondii* and *G. longicalyx* with 2 reps each) from BYU. Those were part of an aborted ATAC-seq project aiming to sample D5, F, A, B and E genomes, and libraries were made by Bob Schmitz et al.

    ls /work/LAS/jfw-lab/ATAC/G.raimondii/ATAC

The D5 dataset was analyzed in comparison with dns-MNase-seq profiles, in order to better understand different types of open chromatin profiles by ATAC-seq and lightly digested MNase-seq 0-130bp fragments.

The analytic pipeline [atac_callPeaks.r](scripts/atac.callPeaks.r) conducted quality QC and trimming of raw fastq reads, mapping against reference genome, and peak calling.

Regarding peak calling:

* [Deal lab protocol](http://plant-plasticity.github.io/resources/3_ATAC-seq%20data%20processing.pdf) used HOMER findPeaks in "region" mode: `findpeaks <tag.directory> -o <output> -gsize <effective.mappable.genome.size_7.1e8> minDist 150 -region`
* [Harvard FAS Informatics ATAC-seq Guidelines](https://informatics.fas.harvard.edu/atac-seq-guidelines.html) uses their own program [Genrich](https://github.com/jsh58/Genrich), and also previously [MACS2 in the order version](https://informatics.fas.harvard.edu/atac-seq-guidelines-old-version.html#peak). 



## DNase-seq

 [dnase.sh](scripts/dnase.sh) 



## Hi-C

The young leaves Hi-C datasets were downloaded from the NCBI Sequence Read Archive database. Analyses were conducted to identify valid chromatin conformation interactions, A/B compartment, TADs, and loops using HOMER with [hic.sh](scripts/hic.sh)  [hic.r](scripts/hic.r) 

* G. raimondii - SRX3051289 
* G. arboreum - SRX3051297
* G. hirsutum - SRX2330709

```bash

bash hic.sh

######## Basic Statistics
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/HiC
#total read pairs
grep 'Total reads processed' */*report.txt
#trimmed and filtered read pairs
grep 'Reads written' */*report.txt 
# mapping input
grep 'reads; of these' */*log
# uniquely mapped read pairs for each end 
grep ' aligned exactly 1 time' */*log
# Total Tags before and after filter -  paired ends counted twice
grep '^genome' */*/tagInfo.txt
# interactions from short
wc -l  */*/short.sorted.txt
# count intra contacts
for j in $( ls */*/short.sorted.txt );    do 
echo $j
awk '$2==$6{c++} END{print c+0}' $j
done

Rscript hic.r
```



## Archived

**Conserved Noncoding Sequence (CNSs)**: The COGE CNS Discovery Pipeline was used to find CNSs between A2 vs D5, At vs Dt, A2 vs At, D5 vs Dt genomes.The resulting common CNSs were used to inquiry whether or not chromatin accessible regions are conversed across genomes. [cns.sh](scripts/cns.sh)
