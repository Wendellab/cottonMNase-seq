## find conserved noncoding regions
# https://bioconductor.org/packages/3.10/bioc/vignettes/CNEr/inst/doc/PairwiseWholeGenomeAlignment.html#last-aligner
# BiocManager::install("CNEr")

# module load r/3.5.0-py2-ufvuwmm bedtools2 last gcc parallel py-pillow/3.2.0-py2-ods24od
# export PATH=$PATH:/work/LAS/jfw-lab/hugj2006/tools/kent/

## soft mask genome FASTA files
# cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/TEannotation
# time RepeatMasker -pa 36 -div 40 -xsmall -lib phase2/G.arboreum_CRI-A2_assembly_v1.0.fasta.EDTA.TElib.fa -cutoff 225 ../A2Du_13.fasta
# time RepeatMasker -pa 36 -div 40 -xsmall -lib phase2/G.raimondii_JGI_221_v2.0.assembly.fasta.EDTA.TElib.fa -cutoff 225 ../Dgenome2_13.fasta
# time RepeatMasker -pa 36 -div 40 -xsmall -lib phase2/phase2/Tx-JGI_G.hirsutum_v1.1.fa.EDTA.TElib.fa -cutoff 225 ../TM1new_26.fasta
# cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/refGenomes/
# mv A2Du_13.fasta.masked LAST/A2.fasta
# mv Dgenome2_13.fasta.masked LAST/D5.fasta

library(CNEr)
## Locate soft masked FASTA files for LAST alignment
assemblyDir <- "refGenomes/LAST"
axtDir <-"refGenomes/LAST"
# "A2Du_13.fasta"
# "Dgenome2_13.fasta"
# "TM1new_26.fasta"
A2 = "A2"
D5 = "D5"
At = "At"
Dt = "Dt"
df.contrast = data.frame(g1 = c(A2,At,A2,D5),g2=c(D5,Dt,At,Dt ))
df.contrast
for(i in 1:4){

G1 = df.contrast$g1[i]
G2 = df.contrast$g2[i]

[CODE]
}

.................
###################################
## 1 Pairwise Geonomic Alignment ##
###################################

## Generate 2bit files from fasta
# faToTwoBit fasta 2bit
system2(command="faToTwoBit", args=c(file.path(assemblyDir, paste0(G1,".fasta")), file.path(assemblyDir, paste0(G1,".2bit"))))
system2(command="faToTwoBit",  args=c(file.path(assemblyDir, paste0(G2,".fasta")), file.path(assemblyDir, paste0(G2,".2bit"))))

## Build the lastdb index for target genome sequence
# lastdb -c refDb ref.fasta
# -c soft-mask lowercase in seq
system2(command="lastdb", args=c("-c", file.path(assemblyDir, G1), file.path(assemblyDir, paste0(G1,".fasta"))))
system2(command="lastdb", args=c("-c", file.path(assemblyDir, G2), file.path(assemblyDir, paste0(G2,".fasta"))))

## Run last aligner
# lastal refDb query.fasta > myalns.maf
lastal(db=file.path(assemblyDir, G1),
       queryFn=file.path(assemblyDir, paste0(G2,".fasta")),
       outputFn=file.path(axtDir, paste0(G1,".",G2,".maf")),
       distance="far", binary="lastal", mc.cores=4L)
lastal(db=file.path(assemblyDir, G2),
       queryFn=file.path(assemblyDir, paste0(G1,".fasta")),
       outputFn=file.path(axtDir, paste0(G2,".",G1,".maf")),
       distance="far", binary="lastal", mc.cores=4L)
# parallel-fasta -j 4 --compress "lastal -a 400 -b 30 -e 6000 -p /tmp/RtmpsCcic5/file71e242d8c54c.lastzMatrix -s 2 -f 1 refGenomes/A2Du_13 " < refGenomes/Dgenome2_13.fasta > refGenomes/LAST/A2Du_13.Dgenome2_13.maf
# lastal refGenomes/Dgenome2_13 refGenomes/A2Du_13.fasta > refGenomes/LAST/Dgenome2_13.A2Du_13.maf
# near 48hrs for unmasked run, 3hr for TE masked run
# cat A2.D5.maf |grep -c '^a score' # 30,894,598
# cat D5.A2.maf |grep -c '^a score' # 36,779,989

## maf to psl
# maf-convert psl my-alignments.maf > my-alignments.psl
system2(command="maf-convert", args=c("psl",
         file.path(axtDir, paste0(G2,".",G1,".maf")),
         ">", file.path(assemblyDir, paste0(G2,".",G1,".psl"))))
system2(command="maf-convert", args=c("psl",
         file.path(axtDir, paste0(G1,".",G2,".maf")),
         ">", file.path(axtDir, paste0(G1,".",G2,".psl"))))
         
### G1 as DB, G2 as query
## Join close alignments
psls =file.path(axtDir, paste0(G1,".",G2,".psl"))
assemblyTarget = file.path(assemblyDir, paste0(G1,".2bit"))
assemblyQuery  = file.path(assemblyDir, paste0(G2,".2bit"))
chains <- axtChain(psls, assemblyTarget=assemblyTarget,
                   assemblyQuery=assemblyQuery, distance="far",
                   removePsl=FALSE, binary="axtChain")
## Sort and combine
allChain <- chainMergeSort(chains, assemblyTarget, assemblyQuery,
            allChain=file.path(axtDir, paste0(G1, ".", G2, ".all.chain")),
            removeChains=FALSE, binary="chainMergeSort")
## Filtering out chains
allPreChain <- chainPreNet(allChain, assemblyTarget, assemblyQuery,
                           allPreChain=file.path(axtDir,paste0(G1,".",G2,".all.pre.chain")),
                           removeAllChain=FALSE, binary="chainPreNet")
## Keep the best chain and add synteny information
netSyntenicFile <- chainNetSyntenic(allPreChain, assemblyTarget, assemblyQuery,
                     netSyntenicFile=file.path(axtDir, paste0(G1, ".", G2, ".noClass.net")),
                     binaryChainNet="chainNet", binaryNetSyntenic="netSyntenic")
## create the .net.axt file from the previous net and chain files.
netToAxt(netSyntenicFile, allPreChain, assemblyTarget, assemblyQuery,
         axtFile=file.path(axtDir, paste0(G1, ".", G2, ".net.axt")),
             removeFiles=FALSE,
             binaryNetToAxt="netToAxt", binaryAxtSort="axtSort")

### G2 as DB, G1 as query
## Join close alignments
psls =file.path(axtDir, paste0(G2,".",G1,".psl"))
assemblyTarget = file.path(assemblyDir, paste0(G2,".2bit"))
assemblyQuery  = file.path(assemblyDir, paste0(G1,".2bit"))
chains <- axtChain(psls, assemblyTarget=assemblyTarget,
                   assemblyQuery=assemblyQuery, distance="far",
                   removePsl=FALSE, binary="axtChain")
## Sort and combine
allChain <- chainMergeSort(chains, assemblyTarget, assemblyQuery,
            allChain=file.path(axtDir, paste0(G2, ".", G1, ".all.chain")),
            removeChains=FALSE, binary="chainMergeSort")
## Filtering out chains
allPreChain <- chainPreNet(allChain, assemblyTarget, assemblyQuery,
                           allPreChain=file.path(axtDir,paste0(G2,".",G1,".all.pre.chain")),
                           removeAllChain=FALSE, binary="chainPreNet")
## Keep the best chain and add synteny information
netSyntenicFile <- chainNetSyntenic(allPreChain, assemblyTarget, assemblyQuery,
                     netSyntenicFile=file.path(axtDir,
                                               paste0(G2, ".", G1, ".noClass.net")),
                     binaryChainNet="chainNet", binaryNetSyntenic="netSyntenic")
## create the .net.axt file from the previous net and chain files.
netToAxt(netSyntenicFile, allPreChain, assemblyTarget, assemblyQuery,
         axtFile=file.path(axtDir, paste0(G2, ".", G1, ".net.axt")),
             removeFiles=FALSE,
             binaryNetToAxt="netToAxt", binaryAxtSort="axtSort")

#######################
## 2 Input alignment ##
#######################
             
## read these axt files and inspect
axtG1G2 = readAxt(file.path(axtDir, paste0(G1, ".", G2, ".net.axt")))
axtG2G1 = readAxt(file.path(axtDir, paste0(G2, ".", G1, ".net.axt")))
axtG1G2
axtG2G1

## need genome size
G1.size = system(paste0("faSize -detailed -tab refGenomes/LAST/",G1,".fasta"),intern=T)
nn = gsub("\t.*","",G1.size)
G1.size = as.numeric(gsub(".*\t","",G1.size))
names(G1.size)=nn
G2.size = system(paste0("faSize -detailed -tab refGenomes/LAST/",G2,".fasta"),intern=T)
nn = gsub("\t.*","",G2.size)
G2.size = as.numeric(gsub(".*\t","",G2.size))
names(G2.size)=nn

pdf(paste0(axtDir,"/plotSynteny.",G1,".",G2,".pdf"))
## Distribution of matched alignments; Given an Axt alignment, plot a heatmap with percentage of each matched alignment
matchDistribution(axtG1G2)
matchDistribution(axtG2G1)

## plot syntenic dot plot
chrs =c(paste0("Chr0",1:9),paste0("Chr1",0:3))
syntenicDotplot(axtG1G2, firstSeqlengths=G1.size, secondSeqlengths=G2.size, firstChrs=chrs, secondChrs=chrs,col=c("blue", "red"), type= "dot")
syntenicDotplot(axtG2G1, firstSeqlengths=G2.size, secondSeqlengths=G1.size, firstChrs=chrs, secondChrs=chrs,col=c("blue", "red"), type= "dot")
dev.off()

##########################
## 3 Filter information ##
##########################

## 3.2 load genomic intervals to be filtered out, including genes and repeats
G1Filter <- readBed(...) 
G2Filter <- readBed(...) 

## 3.3 Create CNE class to store all metadata of running the pipeline for identifying a set of CNEs between two species, including the intermediate and final results. CNE can be created by providing the paths of the twoBit files of two assemblies, and the paths of axt files, with each assembly as reference.
cneG1G2 <- CNE(
  assembly1Fn=file.path(assemblyDir, paste0(G1,".2bit")),
  assembly2Fn=file.path(assemblyDir, paste0(G2,".2bit"))
  axt12Fn=file.path(axtDir, paste0(G1, ".", G2, ".net.axt")),
  axt21Fn=file.path(axtDir, paste0(G2, ".", G1, ".net.axt")),
  cutoffs1=4L, cutoffs2=4L)
cneG1G2

##########################
## 4 CNE identification ##
##########################

## 4.1 Scan alignment to identify CNEs, at two diffrent window sizes 30 and 50 with several similarity criterias (I/C) range from 70% to 100%. 
identities <- c(45L, 48L, 49L)
windows <- c(50L, 50L, 50L)
# since G1 is the reference, use G1Filter as target tFilter 
cneListG1G2 <- ceScan(x=cneG1G2, tFilter=G1Filter, qFilter=G2Filter, window=windows, identity=identities)
# inspect alignments with 70% 
CNE12(cneListG1G2[["45_50"]]) # G1 as reference
CNE21(cneListG1G2[["45_50"]]) # G2 as reference

## 4.2 Merge CNEs. As we perform two rounds of CNE detection with each genome as reference, some conserved elements overlap on both genomes and should be removed. Elements, however, that overlap only on one of the genomes are kept, so that duplicated elements remain distinct.
cneMergedListG1G2 <- lapply(cneListG1G2, cneMerge)

## 4.3 Realignment of CNEs. Some CNEs might be unannotated repeats. To remove them, currently we use blat (W James Kent 2002) to realign each sequence of CNEs against the respective genomes. When the number of matches exceeds a certain threshold, for instance 8, that CNE will be discarded.
# This step can be very time-consuming when the number of CNEs is large. Other alignment methods could also be considered, for example Bowtie2 or BWA (provided that they are installed on the user’s machine)
cneFinalListG1G2 <- lapply(cneMergedListG1G2, blatCNE)

## 4.4 CNE storage and query, using a local SQLite database .
# on individual tables
dbName <- file.path(axtDir, paste0(G1, G2, "CNE.sqlite")),
tableNames <- paste(G1, G2, names(cneFinalListG1G2),
                    sep="_")
for(i in 1:length(cneFinalListG1G2)){
    saveCNEToSQLite(cneFinalListG1G2[[i]], dbName, tableNames[i],
                    overwrite=TRUE)
}
# When querying results from the local SQLite database based on the chr, coordinates and other criterias, a GRanges object is returned.
chr <- "chr6"
start <- 24000000L
end <-  27000000L
minLength <- 50L
tableName <- "A2_D5_45_50"
fetchedCNERanges <- readCNERangesFromSQLite(dbName, tableName, chr, 
                                            start, end, whichAssembly="first",
                                            minLength=minLength)
fetchedCNERanges

## 4.5 CNE length distribution, expected to exhibit power-law distributions, we implemented a function that might be useful in showing interesting patterns in the distribution of the identified elements.
cneGRangePairs <- readCNERangesFromSQLite(dbName=dbName, tableName="A2_D5_45_50",
                      tAssemblyFn=file.path(assemblyDir, paste0(G2,".2bit")),
                      qAssemblyFn=file.path(assemblyDir, paste0(G1,".2bit")))
plotCNEWidth(cneGRangePairs)

## 4.6 Genomic distribution of CNEs along the chromosome
plotCNEDistribution(first(cneGRangePairs)) 

##########################
## 5 CNEs visualization ##
##########################
# https://bioconductor.org/packages/release/bioc/vignettes/CNEr/inst/doc/CNEr.html#creating-a-cne-class
