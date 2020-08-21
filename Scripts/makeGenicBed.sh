# 08/07/18 module bedops no longer supported, run gff2bed locally
  cd mappingD
  head ~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/Dgenome2_13.gene.gff
  gff2bed < ~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/Dgenome2_13.gene.gff >D5.gene.bed
  cd ..
  
  cd mappingM
  head ~/jfw-lab/GenomicResources/archived_resources/gmapdb/AD1TM1/Gossypium_hirsutum_v1.1.gene.gff3
  gff2bed < <(grep 'gene' ~/jfw-lab/GenomicResources/archived_resources/gmapdb/AD1TM1/Gossypium_hirsutum_v1.1.gene.gff3) > TM1.gene.bed
  grep '^A' TM1.gene.bed >TM1.gene.A.bed
  grep '^D' TM1.gene.bed >TM1.gene.D.bed
  cd ..
  
  cd mappingMnew
  zcat ~/jfw-lab/GenomicResources/archived_resources/AD1Saski/annotation/Ghirsutum_458_v1.1.gene_exons.gff3.gz |head
  
  cd mappingA
  head ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Li/A2Li.exons.gff
  gff2bed < <(grep "mRNA" ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Li/A2Li.gene.gff) > A2.gene.bed
  cd ..
   
  cd mappingAnew
  head ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Du2018/1.PacBio-Gar-Assembly-v1.0/G.arboreum.Chr.v1.0.gff
  gff2bed < <(grep "mRNA" G.arboreum.Chr.v1.0.gff) > A2.gene.bed
  cd ..
  
  cd mappingF
  cat mappingA/A2.gene.bed <(sed 's/^Chr/D5_chr/g' mappingDref/mappingD/D5.gene.bed) >F.gene.bed
  grep '^A' F.gene.bed >F.gene.A.bed
  grep '^D' F.gene.bed >F.gene.D.bed
  cd ..
