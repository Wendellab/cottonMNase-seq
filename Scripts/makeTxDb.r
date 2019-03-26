# R
# Generating gene annotations for peak calls and evaluating overall distributions across dataset
# Functional enrichment of gene annotations obtained from peak calls


# Load libraries
library(GenomicFeatures)


####### TM1 saski
txdb<- makeTxDbFromGFF("~/jfw-lab/GenomicResources/archived_resources/AD1Saski/annotation/Ghirsutum_458_v1.1.gene_exons.gff3", format="gff3")
txdb
## restrict features to chromomes, excluding scaffolds
seqlevels(txdb)
seqlevels(txdb)<-seqlevels(txdb.TM1.saski)[1:26]
# original levels are seqlevels0(txdb)
saveDb(txdb, file="refGenomes/txdb.TM1saski.sqlite")
txdb.TM1saski <- loadDb("refGenomes/txdb.TM1saski.sqlite")


####### A2
txdb <- makeTxDbFromGFF("~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Du2018/1.PacBio-Gar-Assembly-v1.0/G.arboreum.Chr.v1.0.gff", format="gff3")
txdb
## restrict features to chromomes, excluding scaffolds
seqlevels(txdb)
seqlevels(txdb)<-seqlevels(txdb)[1:13]
# original levels are seqlevels0(txdb)
saveDb(txdb, file="refGenomes/txdb.A2du.sqlite")
txdb.A2 <- loadDb("refGenomes/txdb.A2du.sqlite")


####### D5
txdb <- makeTxDbFromGFF("~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/Graimondii_221_v2.1.gene_exons.gff3", format="gff3")
txdb
## restrict features to chromomes, excluding scaffolds
seqlevels(txdb)
seqlevels(txdb)<-seqlevels(txdb)[1:13]
# original levels are seqlevels0(txdb)
saveDb(txdb, file="refGenomes/txdb.D5.sqlite")
txdb.D5 <- loadDb("refGenomes/txdb.D5.sqlite")


########## F1 = A2 x D5
system("sed 's/^Chr/A/g' ~/jfw-lab/GenomicResources/archived_resources/gmapdb/A2Du2018/1.PacBio-Gar-Assembly-v1.0/G.arboreum.Chr.v1.0.rn.gff >refGenomes/At.gff")
system("sed 's/^Chr/D/g' ~/jfw-lab/GenomicResources/archived_resources/gmapdb/D5/Graimondii_221_v2.1.gene_exons.gff3 >refGenomes/Dt.gff")
system("cat refGenomes/At.gff refGenomes/Dt.gff >refGenomes/F1.gff ")
txdb <- makeTxDbFromGFF("refGenomes/F1.gff", format="gff3")
txdb
## restrict features to chromomes, excluding scaffolds
seqlevels(txdb)
seqlevels(txdb)<-seqlevels(txdb)[1:26]
# original levels are seqlevels0(txdb)
saveDb(txdb, file="refGenomes/txdb.F1.sqlite")
txdb.F1 <- loadDb("refGenomes/txdb.F1.sqlite")



########## Features
promoters(txdb, upstream=2000, downstream=400)
genes(txdb)
transcripts(txdb)
exons(txdb)
cds(txdb)