# RNA-seq analysis corresponding to cotton MNase-seq datasets
---

## 1. Preprocessing and mapping of RNA-seq datasets

### Locate FASTQ files
For both MNase-seq and RNA-seq experiments, mature leaf tissue was harvested from flowering branches at 5pm, and immediately flash frozen in liquid nitrogen and stored at −80°C. Four *Gossypium* accessions were used including a natural allopolyploid, *G. hirsutum* (AD1) cultivar Acala Maxxa, and models of its A- and D-genome diploid progenitors, *G. arboreum* (A2) and *G. raimondii* (D5), as well as the corresponding interspecific diploid F1 hybrid (A2 × D5).
    
    cd /home/hugj2006/jfw-lab/Projects/MNase-seq/cottonLeaf/RNAseq/rawfastq
    ln -s /lss/research/jfw-lab/RawData/Duplicated_Networks/transcriptomic/SD5-Maxxa* .
    ln -s /lss/research/jfw-lab/RawData/Duplicated_Networks/transcriptomic/SD5-D5* .
    ln -s /lss/research/jfw-lab/RawData/Duplicated_Networks/transcriptomic/SD5-A2* .
    ls

### Quality trimming

    bash runTrimGalore.sh >runTrimGalore.071117.txt 2>&1
    grep -E "Total reads|Reads written" trimmed/*txt >summaryStat/trimStat.txt

### Cotton reference transcriptomes
[CottonGen](https://www.cottongen.org/data/download/genome#Ass) compiles all published cotton genomes, and I used 4 different individual reference genomes for AD1, A2, D5 and A2xD5.

    cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/
    mkdir refTranscriptomes
    cd refTranscriptomes
    module load bowtie2
    
    ### D5_JGI - Paterson et al. 2012 Nature
    wget ftp://ftp.bioinfo.wsu.edu/species/Gossypium_raimondii/JGI_221_G.raimondii_Dgenome/genes/G.raimondii_JGI_221_v2.1.transcripts.fasta.gz
    gunzip G.raimondii_JGI_221_v2.1.transcripts.fasta.gz
    grep -c '>' G.raimondii_JGI_221_v2.1.transcripts.fasta #77267, need only primary
    cat G.raimondii_JGI_221_v2.1.transcripts.fasta |sed 's/ ID=.*$//g'| awk -v F="[.]1$" '$1 ~ F {printf "%s", RS $0} ' RS='>' FS="\n" >D5.transcripts0.fasta
    grep -c '>' D5.transcripts0.fasta #37505
    awk '/^>/ {P=index($0,"Gorai.N")==0} {if(P) print} ' D5.transcripts0.fasta >D5.transcripts.fasta
    grep -c '>' D5.transcripts.fasta #37223
    # ---------------------------------------I used CD regions before, which seemed to lower mapping rate
    ### cp /lss/research/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/Dgenome2_13.cds.fasta >D5.transcripts.fasta
    
    # --------- 2020 -----------------
    
    ### A2_WHU_v1 - Huang et al, 2020
    # Extract transcript references for RSEM and
    optionally build BOWTIE/BOWTIE2/STAR indices
    rsem-prepare-reference --gff3 ../../refGenomes/Garboreum_Shixiya1_WHUv3.0rc.gene.standard.gff3 --gff3-RNA-patterns mRNA ../../refGenomes/A2WHU_13.fasta  A2_WHU
    # 43278 -> 41739 in 13 Chr
    # -------------------
    # only CDs from geneGen
    # wget ftp://ftp.bioinfo.wsu.edu/species/Gossypium_arboreum/WHU_A2_Updated/genes/Garboreum_Shixiya1_WHUv3.0rc.gene.cds.standard.fa.gz 
    # gunzip Garboreum_Shixiya1_WHUv3.0rc.gene.cds.standard.fa.gz
    # mv Garboreum_Shixiya1_WHUv3.0rc.gene.cds.standard.fa A2_WHU.cds.fasta
    
    ### A2_WHU + D5 as ref for A2xD5
    cat A2_WHU.transcripts.fa D5.transcripts.transcripts.fa >F2020.transcripts.fasta
    
    ### AD1_UTX_v2.1 [2020-04-22]	
    wget ftp://ftp.bioinfo.wsu.edu/species/Gossypium_hirsutum/UTX-TM1_v2.1/genes/Ghirsutum_527_v2.1.transcript_primaryTranscriptOnly.fa.gz
    gunzip Ghirsutum_527_v2.1.transcript_primaryTranscriptOnly.fa.gz
    sed 's/ .*$//g' Ghirsutum_527_v2.1.transcript_primaryTranscriptOnly.fa | awk '/^>/ {P=index($0,"Gohir.1Z")==0} {if(P) print} ' >AD1_UTX.transcripts.fasta
    grep -c '>' AD1_UTX.transcripts.fasta #74902

In addition to mapping and expression analysis based on genome-specific transcript, I used the D5 reference transcriptomes based on SNPs 4.0 and 4.1.

    # ls /lss/research/jfw-lab/GenomicResources/pseudogenomes/
    ls /work/LAS/jfw-lab/hugj2006/refTranscriptomes/D5RefPseudo/
    more /work/LAS/jfw-lab/hugj2006/refTranscriptomes/D5RefPseudo/README.txt 
    # for diploids and F1
    grep -c ">" /work/LAS/jfw-lab/hugj2006/refTranscriptomes/D5RefPseudo/A2D5.transcripts.fa
    # 74446
    grep -c ">" /work/LAS/jfw-lab/hugj2006/refTranscriptomes/D5RefPseudo/AtDt.transcripts.fa
    # 74446

Below is the previous 2018 reference genomes followed by RSEM mapping.

    # --------- 2018 -----------------
    ### A2Du2018 - Du et al, 2018 
    zcat  /lss/research/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Du2018/1.PacBio-Gar-Assembly-v1.0/G.arboreum.Chr.v1.0.cds.v1.0.fasta.gz >A2du.transcripts.fasta
    grep -c '>' A2du.transcripts.fasta
    # 40960
    
    ### A2Du + D5 as ref for A2xD5
    cat A2du.transcripts.fasta D5.transcripts.transcripts.fa >F1.transcripts.fasta
    
    ### AD1_458 - Saski et al, 2017
    zcat /lss/research/jfw-lab/GenomicResources/archived_resources/AD1Saski/annotation/Ghirsutum_458_v1.1.transcript_primaryTranscriptOnly.fa.gz >TM1new.transcripts.fasta
    grep -c '>' TM1new.transcripts.fasta
    # 66577

### Bowtie2-RSEM mapping (2018)

    # build rsem-bowtie2 ref
    rsem-prepare-reference D5.transcripts.fasta --bowtie2 D5.transcripts
    rsem-prepare-reference TM1new.transcripts.fasta --bowtie2 TM1new.transcripts
    rsem-prepare-reference F1.transcripts.fasta --bowtie2 F1.transcripts
    rsem-prepare-reference A2du.transcripts.fasta --bowtie2 A2du.transcripts
    rsem-prepare-reference F2020.transcripts.fasta --bowtie2 F2020.transcripts
    rsem-prepare-reference A2_WHU.transcripts.fasta --bowtie2 A2_WHU.transcripts
    rsem-prepare-reference AD1_UTX.transcripts.fasta --bowtie2 AD1_UTX.transcripts
    
    bash runBowtie2.D.sh >runBowtie2.D.101318.txt 2>&1
    bash runBowtie2.A.sh >runBowtie2.A.101318.txt 2>&1
    bash runBowtie2.F1.sh >runBowtie2.F1.101318.txt 2>&1
    bash runBowtie2.AD1.sh >runBowtie2.AD1.101218.txt 2>&1
    bash runBowtie2.Dref.sh >runBowtie2.Dref.101918.txt 2>&1

Mapping rate using D5 reference transcriptomes based on SNPs 4.0 and 4.1.


    grep 'reads; of these:' mappingDref/*log | paste - <(grep 'concordantly 0 times' mappingDref/*log) <(grep 'overall alignment rate' mappingDref/*log) >summaryStat/mapDref.txt

Mapping rate using individual reference transcriptomes.

    grep 'reads; of these:' mappingIndiv/*log | paste - <(grep 'concordantly 0 times' mappingIndiv/*log) <(grep 'overall alignment rate' mappingIndiv/*log) >summaryStat/mapIndiv2018.txt

### Kallisto mapping (2020)

```bash
# build Kallisto ref
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/RNAseq/refTranscriptomes
kallisto index -i D5.transcripts.kallisto D5.transcripts.fasta
kallisto index -i F2020.transcripts.kallisto F2020.transcripts.fasta
kallisto index -i A2_WHU.transcripts.kallisto A2_WHU.transcripts.fa
kallisto index -i AD1_UTX.transcripts.kallisto AD1_UTX.transcripts.fasta
kallisto index -i A2D5.transcripts.kallisto /work/LAS/jfw-lab/hugj2006/refTranscriptomes/D5RefPseudo/A2D5.transcripts.fa
kallisto index -i AtDt.transcripts.kallisto /work/LAS/jfw-lab/hugj2006/refTranscriptomes/D5RefPseudo/AtDt.transcripts.fa

# Kallisto mapping: 12 samples each map to indivdual and Dref references, e.g.
for j in $( ls trimmed/ | grep 'SD5-Maxxa.*_1_val_1.fq.gz'); do
echo $j
echo ${j%%_1_val*}_2_val_2.fq.gz
kallisto quant -t $ncores -i refTranscriptomes/AD1_UTX.transcripts.kallisto -o MappingIndiv/${j%%_*} trimmed/$j trimmed/${j%%_1_val*}_2_val_2.fq.gz 2>MappingIndiv/${j%%_*}.log.txt
kallisto quant -t $ncores -i refTranscriptomes/AtDt.transcripts.kallisto -o MappingDref/${j%%_*} trimmed/$j trimmed/${j%%_1_val*}_2_val_2.fq.gz 2>MappingDref/${j%%_*}.log.txt
done

# summary
grep "pseudoaligned" MappingIndiv/*log.txt |cut -f3,5 >summaryStat/mapIndiv2020.txt
grep "pseudoaligned" MappingDref/*log.txt |cut -f3,5 >summaryStat/mapDref2020.txt
```

## 2. Analyze gene expression profiles

Ortholog/homoeolog relationships need to be used to link genes from individual genomes

```
# A2_WHU, D5, AD1_UTX
ls /work/LAS/jfw-lab/hugj2006/cottonLeaf/orthohomoeologQuadruplets101218.txt

# A2_Du, D5, TM1 sacki by Justin
ls /work/LAS/jfw-lab/hugj2006/cottonLeaf/orthohomoeologQuadruplets101218.txt
```

The r code [expTests.r](scripts/expTests.r) performs following:

* Prepare raw count and rpkm tables from mapping results
* Differential expression analysis
* Cis-trans analysis
* Impact of genome evolution: Hr, Pr, Wr
* Expression dominance analyisis
