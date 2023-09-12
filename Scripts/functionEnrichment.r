# Prep functional annotation for AD1-UTX and OG
# Modified based on https://github.com/Wendellab/CisTransRegulation/blob/master/functionEnrichment.r
# ----AND "GOandKEGG.0510.R" by XXP
# Guanjing Hu
# 05-12-2022
setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/4.RNAseqAnalysis")


#########################
## clusterProfile prep ##
#########################
# plots: http://bioconductor.org/packages/devel/bioc/vignettes/enrichplot/inst/doc/enrichplot.html

# prep reference database
# AD1 TM1 2019 saski: "Genes66610.transcript.description.txt"
library(clusterProfiler)
library("goProfiles")
library(tidyr)
library(GO.db)
library(topGO)

## prep GO enrichment function, plot top 10 enriched each for MF, CC and BP and save table
# require custom geneID2GO
topGOenrich = function(geneL, refL, filename="topGOenrich"){
    if(!is.null(filename)){pdf(paste0(filename,".pdf"))}
    geneList <- factor(as.integer(refL %in% geneL))
    names(geneList) <- refL
    for(on in c("MF","BP","CC"))
    {
        print(on)
        # Make topGO object
        GOdata <- new("topGOdata", ontology = on, allGenes = geneList, nodeSize = 5, annot = annFUN.gene2GO, gene2GO = geneID2GO)
        # fisher test
        result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
        results.table <- GenTable(GOdata, result, topNodes = length(result@score))
        # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, FDR <= 5%
        results.table$qval.bh<-p.adjust(results.table[,"result1"],method="BH")
        # label ontology type
        results.table$ont<-on
        # reduce results to GO terms passing Benjamini-Hochberg multiple hypothesis corrected pval <= 0.05, consider FDR <= 5% in future
        keep <- results.table[as.numeric(results.table[,"result1"])<0.01,]
        if(exists("dfL")) dfL<- rbind(dfL, keep)
        if(!exists("dfL")) dfL<-keep
        # draw figure for GO terms pval<=0.05 before FDR correction
        if(is.na(sigNo<-length(keep$ontology))){next}
        showSigOfNodes(GOdata, score(result), firstSigNodes = sigNo, useInfo = "all")
        mtext(on, line=-1)
    }
    if(!is.null(filename)){
        dev.off()
        write.csv(dfL,file=paste0(filename,".csv"))}
    print(table(dfL$ont))
    return(dfL)
}


## prep GO enrichment function, plot top 10 enriched each for MF, CC and BP and save table
# require custom godb and goterm2name
GOenrich = function(geneL, refL, filename="GOenrich"){
    if(!is.null(filename)){pdf(paste0(filename,".pdf"))}
    for(o in c("MF","BP","CC")){
        en <- enricher(geneL, pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=refL, minGSSize = 5, maxGSSize = 500, qvalueCutoff = 0.05, TERM2GENE = godb[godb$ont==o,1:2], TERM2NAME = goterm2name)
        # print(en)
        print( dotplot(en, showCategory=10,font.size = 6, title=o) )
        if(nrow(en)>0){
            print( emapplot(en, showCategory=30,font.size = 6, title=o) )
            df=as.data.frame(en)
            df$ont=o
            df$level2 = df$ID %in% getGOLevel(o, 2)
            df$level3 = df$ID %in% getGOLevel(o, 3)
            df$level4 = df$ID %in% getGOLevel(o, 4)
            if(o=="MF"){dfL=df}
            else if(!exists("dfL")){dfL=df}
            else{dfL=rbind(dfL,df)}
        }
    }
    if(!is.null(filename)){
        dev.off()
        write.csv(dfL,file=paste0(filename,".csv"))}
    print(table(dfL$ont))
    return(dfL)
}
# test
# GOenrich(geneL,refL, filename="GOenrich")

GOcompare = function(geneL, refL, filename="GOcompare"){
    if(!is.null(filename)){pdf(paste0(filename,".pdf"))}
    for(o in c("MF","BP","CC")){
        print(o);
        en=tryCatch(compareCluster(geneL, fun = "enricher", pvalueCutoff = 0.05, pAdjustMethod = "BH", universe=refL, minGSSize = 5, maxGSSize = 500, qvalueCutoff = 0.05, TERM2GENE = godb[godb$ont==o,1:2], TERM2NAME = goterm2name), error=function(e){cat("ERROR :",conditionMessage(e), "\n")}, finally={en=NULL})
        if(!is.null(en)){
            print( dotplot(en, showCategory=10,font.size = 6, title=o) )
            df=as.data.frame(en)
            df$ont=o
            df$level2 = df$ID %in% getGOLevel(o, 2)
            df$level3 = df$ID %in% getGOLevel(o, 3)
            df$level4 = df$ID %in% getGOLevel(o, 4)
            if(!exists("dfL")){dfL=df}else{dfL=rbind(dfL,df)}
        }
    }
    if(!exists("dfL")){dfL=NULL}
    if(!is.null(filename)){
        dev.off()
        write.csv(dfL,file=paste0(filename,".csv"))}
    print(table(dfL[,c("ont","Cluster")]))
    return(dfL)
}

GOcompare.formula = function(formula, data, filename="GOcompare"){
    if(!is.null(filename)){pdf(paste0(filename,".pdf"))}
    for(o in c("MF","BP","CC")){
        en <- compareCluster(formula, data, fun = "enricher", pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500, qvalueCutoff = 0.05, TERM2GENE = godb[godb$ont==o,1:2], TERM2NAME = goterm2name)
        print( dotplot(en, showCategory=10,font.size = 6, title=o) )
        df=as.data.frame(en)
        df$ont=o
        df$level2 = df$ID %in% getGOLevel(o, 2)
        df$level3 = df$ID %in% getGOLevel(o, 3)
        df$level4 = df$ID %in% getGOLevel(o, 4)
        if(o=="MF"){dfL=df}else{dfL=rbind(dfL,df)}
    }
    print(dotplot.df(dfL,title=filename))
    if(!is.null(filename)){
        dev.off()
        write.csv(dfL,file=paste0(filename,".csv"))}
    print(table(dfL[,c("ont","Cluster")]))
    return(dfL)
}
# result=GOcompare.formula(ID~rd,data=rd,filename="GOcompare.RD")

# dependent functions
getGOLevel <- function(ont, level) {
    switch(ont,
    MF = {
        topNode <- "GO:0003674"
        Children <- GOMFCHILDREN
    },
    BP = {
        topNode <- "GO:0008150"
        Children <- GOBPCHILDREN
    },
    CC = {
        topNode <- "GO:0005575"
        Children <- GOCCCHILDREN
    }
    )
    
    max_level <- max(level)
    if (any(level == 1)) {
        all_nodes <- topNode
    } else {
        all_nodes <- c()
    }
    
    Node <- topNode
    for (i in seq_len(max_level-1)) {
        Node <- mget(Node, Children, ifnotfound=NA)
        Node <- unique(unlist(Node))
        Node <- as.vector(Node)
        Node <- Node[!is.na(Node)]
        if ((i+1) %in% level) {
            all_nodes <- c(all_nodes, Node)
        }
    }
    return(all_nodes)
}
dotplot.df=function(result, title="")
{
    require(ggplot2)
    require(DOSE) # theme_dose
    parse_ratio <- function(ratio) {
        gsize <- as.numeric(sub("/\\d+$", "", as.character(ratio)))
        gcsize <- as.numeric(sub("^\\d+/", "", as.character(ratio)))
        return(gsize/gcsize)
    }
    ratio=parse_ratio(result$GeneRatio)
    # give color to ont
    a <- ifelse(result$ont == "MF", "gold", ifelse(result$ont == "BP","steelblue","salmon"))[order(result$Description)]
    ggplot(result, aes_string(x="Cluster", y="Description", size=ratio, color="p.adjust")) + geom_point() + scale_color_continuous(low="red", high="blue", name = "p.adjust", guide=guide_colorbar(reverse=TRUE)) + theme_dose(font.size=6) + scale_size(range=c(3, 8))+ theme(axis.text.y = element_text(colour = a)) + ylab(NULL) + ggtitle(title)
    #scale_color_gradientn(name = "p.adjust", colors=sig_palette, guide=guide_colorbar(reverse=TRUE)) +
    #+ facet_grid(ont~.)
}
# prep GO levels
library("goProfiles")
level1<-c(getGOLevel("MF",1),getGOLevel("BP",1),getGOLevel("CC",1))
level2<-c(getGOLevel("MF",2),getGOLevel("BP",2),getGOLevel("CC",2))
level3<-c(getGOLevel("MF",3),getGOLevel("BP",3),getGOLevel("CC",3))
level4<-c(getGOLevel("MF",4),getGOLevel("BP",4),getGOLevel("CC",4))
level5<-c(getGOLevel("MF",5),getGOLevel("BP",5),getGOLevel("CC",5))
level6<-c(getGOLevel("MF",6),getGOLevel("BP",6),getGOLevel("CC",6))
level7<-c(getGOLevel("MF",7),getGOLevel("BP",7),getGOLevel("CC",7))

ogQ<-read.table("orthohomoeolog052421.txt", sep="\t", header=FALSE,stringsAsFactor = FALSE)
ogQ$V3<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",ogQ$V3)))
ogQ$V4<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",ogQ$V4)))
# prep GOs in "term2og"
og<-ogQ[,c(1,3:4,2,5)]
rownames(og)<-ogQ$V1
names(og)<-c("ogName","At","Dt","A2","D5")


#################################################################################
# 1. AD1 TM1 UTX, annotated by XXP "Ghirsutum_527_v2.1.protein_xxp.annotations" #
#################################################################################
egg1 <- read.csv("Ghirsutum_527_v2.1.protein_xxp.annotations.txt", sep = "\t")
dim(egg1)  #[1] 101038     22
length(unique(unlist(strsplit(egg1$GOs,","))) ) #11233
# extract only primary longest transcripts
gene75736<-read.table("75736.txt", sep = "\t",header = T)
head(gene75736)
dim(egg<-merge(gene75736,egg1,by="query_name")) #69,957
# get gene ID
length(unique(egg$geneName<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",egg$query_name)))))#69,957
# how many has GOs
table(egg$GOs!="")
# FALSE  TRUE
# 34301 35656

# prep GOs in "term2gene"
gomap <- egg[,c("GOs","geneName")]
names(gomap)<-c("GO","geneName")
gomap<-separate_rows(gomap, GO,sep=",")
term2gene=data.frame(gomap)

# in order to seperate ontology first
# godb<-buildGOmap(gomap)
godb<-term2gene
ont=go2ont(godb$GO)
godb$ont = ont$Ontology[match(godb$GO,ont$go_id)]
length(unique(godb$geneName)) # 35,656 genes with GOs
length(unique(godb$GO))  # 11,232 unique GO

# prep ALL GO term2name
goterms <- Term(GOTERM)
dim(goterm2name <- data.frame(term=names(goterms),name=goterms)) #  44,509

# topGO
library(topGO)
temp <- aggregate(GO~geneName, term2gene, I)
geneID2GO <- temp$GO
names(geneID2GO) <-temp$geneName
head(geneID2GO)

# save into rdata
save(egg, godb, goterm2name,geneID2GO, topGOenrich, GOenrich, GOcompare, getGOLevel, dotplot.df,geneID2GO, file="funcEnrich.utx_xxp.rdata")

# TEST
geneL<-term2gene$geneName[1:1000]
refL <-term2gene$geneName
## Test clusterProfiler enrichment, noting that en0 didn't seperate ontology first
GOenrich(geneL, refL, filename="testCP")
en0 <- enricher(gene=geneL, universe=refL,TERM2GENE = term2gene, TERM2NAME = goterm2name, pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05 )
print( dotplot(en0, showCategory=50,font.size = 6, title="test") )
en0  # collect more information
df<-as.data.frame(en)
terms <- go2term(df$ID)
onts <- go2ont(df$ID)
df$ont <- onts[match(df$ID,onts$go_id),"Ontology"]
df$term <- terms[match(df$ID,terms$go_id),"Term"]
df$level2 <- df$ID %in% level2
df$level3 <- df$ID %in% level3
df$level4 <- df$ID %in% level4
## Test clusterProfiler enrichment, now seperate ontology first
en1 <- GOenrich(geneL, refL, filename="testCP") # BP-6
## Test topGO
library(topGO)
en2 <- topGOenrich(geneL, refL, filename="testTG") #BP-15, CC-1, MF-7


## Prep for OG
head(og)
dim(egg[egg$GOs!="",])
eggAt<-egg[match(og$At,egg$geneName),c("query_name","GOs","seed_eggNOG_ortholog","eggNOG.free.text.desc.")]
eggDt<-egg[match(og$Dt,egg$geneName),c("query_name","GOs","seed_eggNOG_ortholog","eggNOG.free.text.desc.")]
table(eggAt$GOs==eggDt$GOs&eggAt$GOs!="")
# FALSE  TRUE
# 10247 12383
table(eggAt$seed_eggNOG_ortholog==eggDt$seed_eggNOG_ortholog)
# FALSE  TRUE
#   522 22094
table(eggAt$eggNOG.free.text.desc.==eggDt$eggNOG.free.text.desc. & eggDt$eggNOG.free.text.desc.!="")
# FALSE  TRUE
#   60 20677

##  How many OGs have egg Annotation
table(og$At %in% egg$geneName | og$Dt %in% egg$geneName)
# FALSE  TRUE
# 231 22658

og$GO<-paste0(eggAt$GOs,",",eggDt$GOs)
gomap <- egg[,c("GOs","geneName")]
names(gomap)<-c("GO","geneName")
gomap<-separate_rows(og[,c("GO","ogName")], GO,sep=",")
dim(godb<-unique(data.frame(gomap))) #731,070
dim(godb<-godb[grepl("GO",godb$GO),]) # after removing empty, 720583
ont=go2ont(godb$GO)
godb$ont = ont$Ontology[match(godb$GO,ont$go_id)]
length(unique(godb$ogName)) # 12,440 OGs with GOs
length(unique(godb$GO))  # 10,910 unique GO
temp <- aggregate(GO~ogName, godb, I)
geneID2GO <- temp$GO
names(geneID2GO) <-temp$ogName
head(geneID2GO)
length(geneID2GO)# 12,440 OGs with GOs

save(og,eggAt,eggDt, godb, geneID2GO, goterm2name,GOenrich,topGOenrich, GOcompare, getGOLevel, dotplot.df, file="funcEnrich.OG_xxp.rdata")

library(simplifyEnrichment)
go_id = godb$GO[godb$ont=="MF"]
# the GO terms should only belong to one same ontology (i.e., BP, CC or MF).
mat = GO_similarity(go_id[1:300])
simplifyGO(mat)

simplifyGOFromMultipleLists()

#######################################
# 2.AD1 TM1 UTX, annotated by NOVOgen #
#######################################

dim(go1 <- read.csv("go.txt", sep = "\t",header=T)) #[1] 804,411     22
length(unique(go1$gene_id)) #40,500
length(unique(go1$go_id))# 2662,WTH so few???
go1$gene_id<-gsub(".v2.1.*","",go1$gene_id)
dim(go1<-unique(go1)) #802,069      4
# extract only primary longest transcripts
gene75736<-read.table("75736.txt", sep = "\t",header = T)
ids<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",gene75736$query_name)))
head(ids)
dim(go2<-go1[go1$gene_id %in% ids,] )# 787,600      4
length(unique(go2$gene_id)) #39,507
#godb
godb<-data.frame(GO=go2$go_id,geneName=go2$gene_id,ont=go2$go_ontology)
length(unique(godb$geneName)) # 339507 genes with GOs
length(unique(godb$GO))  # 2645 unique GO
en3 <- GOenrich(geneL, refL, filename="testNovo")  ## CC-2, NOT GOOD

# og
ogL<-data.frame(ogName=c(og$ogName, og$ogName),gene=c(og$At,og$Dt))
go2$ogName<-ogL$ogName[match(go2$gene_id,ogL$gene)]
#  prep
head( godb_gene <- data.frame(GO=go2$go_id,geneName=go2$gene_id,ont=go2$go_ontology) )
godb_og <- data.frame(GO=go2$go_id,ogName=go2$ogName,ont=go2$go_ontology)
dim(godb_og<-godb_og[!is.na(godb_og$ogName),]) # 536871      3
length(unique(godb_og$ogName)) #13480
length(unique(godb_og$GO)) # 2602
dim( goterm2name <- unique(data.frame(term=go2$go_id, name=go2$go_term)) ) #2645
save(godb,godb_gene,godb_og, goterm2name,GOenrich, GOcompare, getGOLevel, dotplot.df, file="funcEnrich.novogen.rdata")

###################
# 3.D5 annotation #
###################
# Annotated moduls with topGO
# gene-to-GO mappings

library("topGO")
#geneID2GO <- readMappings(file = "~/Research/Sequence Database/G. raimondii Genome/Function_v2.1/cotton.v2.1.stable.txt")
#names(geneID2GO)<-gsub("[.]1$","",names(geneID2GO) )
#save(gsc, geneID2GO,  goraiUniverse, file="~/Dropbox/Scripts/cottonGOenrich.RData")
load('~/Nutstore Files/GHuLab-Me/1.Research/Scripts/cottonGOenrich.RData')->l;l
# "geneID2GO"       "go4profile"  "goraiUniverse"   "gsc"             "level1_BP_terms" "level1_MF_terms"  "level2_BP_terms" "level2_MF_terms" "level3_BP_terms" "level3_MF_terms"   "level4_BP_terms" "level4_MF_terms" "level5_BP_terms" "level5_MF_terms"
length(unique(as.character(go4profile$GeneID))) # 52821
length(unique(as.character(go4profile$GOID))) # 1187
length(geneID2GO) # 20231
length(unique(unlist(geneID2GO))) # 1195

## Prep for OG
ogQ<-read.table("orthohomoeolog052421.txt", sep="\t", header=FALSE,stringsAsFactor = FALSE)
ogQ$V3<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",ogQ$V3)))
ogQ$V4<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",ogQ$V4)))
# prep GOs in "term2og"
og<-ogQ[,c(1,3:4,2,5)]
rownames(og)<-ogQ$V1
names(og)<-c("ogName","At","Dt","A2","D5")
# check representation
table(gsub("[.]1$","",og$D5) %in% go$GeneID)
# FALSE  TRUE
# 19009  3880 few og has GO based on D5, really not good

#######################
## UTX + A2Du +D5JGI ##
#######################

# https://www.cottongen.org/data/download/genome_diploid_A_nd_D5: GO from CottonGen InterProScan
a<-read.table("refGenomeGO/A2_WHU_v1_a1_genes2Go.txt",sep="\t")
dim(a) # 27000     4
length(unique(a[,1])) # 11561
length(unique(a[,2])) # 1760

d<-read.table("refGenomeGO/G.raimondii_JGI_221_v2.1_genes2GO.txt",sep="\t")
dim(d) # 79789     4
length(unique(d[,1])) # 24734
length(unique(d[,2])) # 1451

ad<-read.table("refGenomeGO/Gh_TM1_UTX_v2.1_genes2Go.txt",sep="\t")
dim(ad) # 65099     4
length(unique(ad[,1])) # 27685
length(unique(ad[,2])) # 2068

## Prep for OG
ogQ<-read.table("orthohomoeolog052421.txt", sep="\t", header=FALSE,stringsAsFactor = FALSE)
ogQ$V3<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",ogQ$V3)))
ogQ$V4<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",ogQ$V4)))
# prep GOs in "term2og"
og<-ogQ[,c(1,3:4,2,5)]
rownames(og)<-ogQ$V1
names(og)<-c("ogName","At","Dt","A2","D5")

ad$id<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",ad$V1)))
a$id<-gsub("_WHU-A2_v1","",a$V1)
d$id<-d$V1

table(og$At %in% ad$id) #9109
table(og$Dt %in% ad$id) #3721
table(og$A2 %in% a$id) #7186
table(og$D5 %in% d$id)  #7472

## FUCK, annotation still doesn't look too good.
# AD1 TM1 2019 saski: "Genes66610.transcript.description.txt"
# grep "PF00125" /Users/guanjinghu/Desktop/2.Projects/Project_MNase-seq/4.RNAseqAnalysis/Genes66610.transcript.description.txt > /Users/guanjinghu/Desktop/2.Projects/Project_MNase-seq/4.RNAseqAnalysis/histones.Genes66610.transcript.description.txt
# manual check to keep only histone H2A/H2B/H3

histone<- read.table("histones.Genes66610.transcript.description.txt", sep="\t")
length(hid<-gsub(".1$","",histone$V1)) #112 histones
table(og$At %in% hid) #31
table(og$Dt %in% hid) #29
rownames(histone)<-hid
og$descAt<-histone[og$At,"V12"]
og$descDt<-histone[og$Dt,"V12"]

# D5: JGI v2 download from Phytozome "Graimondii_221_v2.1.annotation_info.txt"
# grep 'PF00125' Graimondii_221_v2.1.annotation_info.txt |cut -f2|sort|uniq |wc -l    #53
# grep 'PF00125' Graimondii_221_v2.1.annotation_info.txt >histone.D5.txt
histoneD<-read.table("histone.D5.txt", sep="\t")
table(og$D5 %in% histoneD$V3)  #32
rownames(histoneD)<-histoneD$V3
og$descD5<-histoneD[og$D5,"V12"]
og$descD5m<-histoneD[og$D5,"V13"]

# histoneDB
yy<-read.table("AD1utx.histones.txt", sep="\t")
zz<-read.table("D5.histones.txt", sep="\t")
table(og$At %in% yy$V2)  #56
table(og$Dt %in% yy$V2)  #57
table(og$D5 %in% zz$V2)  #57
rownames(yy)<-yy$V2
rownames(zz)<-zz$V2
og$hD5<-zz[og$D5,"V1"]
og$hAt<-yy[og$At,"V1"]
og$hDt<-yy[og$Dt,"V1"]

# histonProbst
hh<-read.table("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/4.RNAseqAnalysis/histoneProbst/Ath.histones.Probst.txt",sep="\t",head=T)
mm<-read.table("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/4.RNAseqAnalysis/histoneProbst/cottonAth.txt",head=T,sep="\t")
mm$Ath<-gsub("[.].*","",mm$Query.ID)
dim(mh<-merge(mm,hh,by="Ath"))
mh$loci <- gsub("_",".",gsub("[.].*","",gsub("Gohir.","Gohir_",gsub("Gorai.","Gorai_",mh$Protein))))
dim(mhh<-unique(mh[,c("loci","Histone")])) #260
dim(mhhh<-unique(mh[,c("loci","Histone","Variant")])) #1937
table(mhh$loci %in% og$At) # 56
table(mhh$loci %in% og$Dt) #57
write.table(mhh, "cottonID.class.txt",sep="\t",row.names=F)
cbind(mhh[ogH$At,"Histone"],mhh[ogH$Dt,"Histone"],ogH$hD5) # consistent classification
# no difference from above histoneDB test, think I am ok for now; Why 260 genes awaits future gene famility analysis!!!!
hea
ogh<-rownames(og)[og$At %in% hid | og$Dt %in% hid | og$D5 %in% histoneD$V3 | og$At %in% yy$V2 |og$Dt %in% yy$V2 | og$D5 %in% zz$V2]
ogH<-og[ogh,]
ogH$class<-gsub("[.].*","",gsub(".*_","",ogH$hD5))
save(ogH,"histoneOG.rdata")

load("AD1GNT2021rRNA/diffExprPatternTPM.rdata")
tpm<-totalTPM[ogH[,"D5"],]
boxplot(tpm, las=2)

##############
## Analysis ##
##############

library(clusterProfiler)
library("goProfiles")
library(tidyr)
library(GO.db)

setwd("~/Nutstore Files/GHuLab-Me/2.Projects/Project_MNase-seq/4.RNAseqAnalysis")
con<- read.table("conserved.regPattern.txt",header=T, sep="\t")

load("funcEnrich.OG_xxp.rdata")
result=GOcompare.formula(V1~A,data=con,filename="GOxxp/A")
result=GOcompare.formula(V1~B,data=con,filename="GOxxp/B")
result=GOcompare.formula(V1~Bp,data=con,filename="GOxxp/Bp")
result=GOcompare.formula(V1~Hr,data=con,filename="GOxxp/Hr")
result=GOcompare.formula(V1~Wr,data=con,filename="GOxxp/Wr")
result=GOcompare.formula(V1~Pr,data=con,filename="GOxxp/Pr")
result=GOcompare.formula(V1~cistrans,data=con,filename="GOxxp/cistrans")
result=GOcompare.formula(V1~F1.dominance,data=con,filename="GOxxp/F1.dominance")
result=GOcompare.formula(V1~AD1.dominance,data=con,filename="GOxxp/dominance")
# # # ont   - A<0 A=0 A>0
#   BP 33   0   0  11
#   CC 18   0   0  18
#   MF 10   0   1   0
# # # # # # ont   - B<0 B=0 B>0
#   BP 43   0   2   0
#   CC 40   0   2   0
#   MF 12   0   5   0
# # # # # # # # # ont  - Hr<0 Hr=0 Hr>0
#   BP 8    0    0    0
#   CC 7    0    0    0
#   MF 2    2    0    0
# # # # # # ont   - Wr<0 Wr=0 Wr>0
#   BP 23    0    0    0
#   CC 16    0    0    0
#   MF  3    0    1    0
# # # # # # ont   - Pr<0 Pr=0 Pr>0
#   BP 10    1    0    0
#   CC 10    2    0    0
#   MF  4    0    1    0
# # # # # # ont   - 1.Cis only 2.Trans only 3.Cis+Trans: enhancing 4.Cis+Trans: compensating 5.Compensatory
#   BP 49          0            0                     16                        24              0
#   CC 29          0            0                      5                         1              0
#   MF 17          0            0                      9                         5             12
#     Cluster
# ont  6.Conserved 7.Ambiguous
#   BP           4          57
#   CC           3          15
#   MF           5          15
# # # # # # ont    - A-dominant Additivity D-dominant Other non-additivity Transgressive Down Transgressive Up
#   BP  29          0         69         28                    0                 38              136
#   CC  46          0         16         25                    0                 15                8
#   MF   3          2          3          6                    0                 11                6
# # # # # # ont   - A-dominant Additivity D-dominant Other non-additivity Transgressive Down Transgressive Up
#   BP 94          2         16          3                    0                 17                0
#   CC 25          0          5         15                    0                  0                0
#   MF 14          2          0          0                    2                  0                0

## ressulting GOs can be summarized
library(simplifyEnrichment)
# the GO terms should only belong to one same ontology (i.e., BP, CC or MF).
res=result[result$ont=="MF",]
lt=split(res,res$A)
simplifyGOFromMultipleLists(lt,padj_cutoff = 0.05)

# use p value
lt2 = lapply(lt, function(x) structure(x$p.adjust, names = x$ID))


refL<-og$ogName
##topGO analysis + simplifyEnrich plot
for(tag in c("A","B","Bp","Hr","Pr","Wr","cistrans","F1.dominance","AD1.dominance")){
    # tag="A"
    print(tag)
    gL<-  split(con$V1,con[,tag])
    nn<-names(gL)
    resL<-list()
    for(r in 1:length(gL)){
        resL[[r]]<- topGOenrich(geneL=gL[[r]], refL, filename=paste0("GOxxp/",tag,".",nn[r]))
    }
    names(resL)<-nn
    resL[nn=="-"]<-NULL
    pdf(paste0("GOxxp/simplify.",tag,".pdf"))
    for(o in c("MF","BP","CC")){
        resL1<-lapply(resL, function(x)x[x$ont==o,])
        lt<-lapply(resL1, function(x) structure(as.numeric(x$result1), names = x$GO.ID))
        lt<-lapply(lt,function(x)x[!is.na(x)])
        simplifyGOFromMultipleLists(lt,padj_cutoff = 0.01)
    }
    dev.off()
}



----

load("funcEnrich.novogen.rdata")
godb<-godb_og
result=GOcompare.formula(V1~A,data=con,filename="GOnovo/A")
result=GOcompare.formula(V1~B,data=con,filename="GOnovo/B")
result=GOcompare.formula(V1~Bp,data=con,filename="GOnovo/Bp")
result=GOcompare.formula(V1~Hr,data=con,filename="GOnovo/Hr")
result=GOcompare.formula(V1~Wr,data=con,filename="GOnovo/Wr")
result=GOcompare.formula(V1~Pr,data=con,filename="GOnovo/Pr")
result=GOcompare.formula(V1~cistrans,data=con,filename="GOnovo/cistrans")
result=GOcompare.formula(V1~F1.dominance,data=con,filename="GOnovo/F1.dominance")
result=GOcompare.formula(V1~AD1.dominance,data=con,filename="GOnovo/dominance")
# ont   - A<0 A=0 A>0
#   BP 10   0   0   1
#   CC  4   0   0   9
#   MF  2   2   0   0
# ont   - B<0 B=0 B>0
#   BP  8   0   5   0
#   CC  7   0   3   0
#   MF 14   0   6   0
# ont   - Bp<0 Bp=0 Bp>0
#   BP 11    0    0    0
#   CC  4    0    2    0
#   MF 11    0    0    0
# ont  - Hr<0 Hr=0 Hr>0
#   BP 5    0    0    0
#   CC 4    0    0    0
#   MF 2    0    0    0
# ont  - Wr<0 Wr=0 Wr>0
#   BP 5    0    0    0
#   CC 5    0    3    0
#   MF 3    0    0    0
# ont  - Pr<0 Pr=0 Pr>0
#   BP 5    0    0    0
#   CC 4    0    1    0
#   MF 2    0    0    0
# ont   - 1.Cis only 2.Trans only 3.Cis+Trans: enhancing 4.Cis+Trans: compensating 5.Compensatory
#   BP  7          0            0                      2                        12              0
#   CC  4          0            0                      0                         0              0
#   MF  4          6            0                      8                         3              0
#     Cluster
# ont  6.Conserved 7.Ambiguous
#   BP           4           5
#   CC           0           2
#   MF           5          11
# ont   - A-dominant Additivity D-dominant Other non-additivity Transgressive Down Transgressive Up
#   BP  5          0         15         10                    0                  0                0
#   CC 18          0          5         14                    0                  0                0
#   MF  7          0          5          2                    0                  5               12
# ont   - A-dominant Additivity D-dominant Other non-additivity Transgressive Down Transgressive Up
#   BP 15          0          8          1                    0                  0                0
#   CC  7          0          3         11                    0                  0                0
#   MF 11          0         12          0                    0                  0                0
#


####################
## Flowering gene ##
####################
library(gdata)
load("funcEnrich.OG_xxp.rdata")
con<- read.table("conserved.regPattern.txt",header=T, sep="\t")
table(con$V1==og$ogName)

# append XXP annotation
con$eggAt.og<-eggAt$seed_eggNOG_ortholog
con$eggAt.og.desc<-eggAt$eggNOG.free.text.desc.
con$eggDt.og<-eggDt$seed_eggNOG_ortholog
con$eggDt.og.desc<-eggDt$eggNOG.free.text.desc.

# append UTX At blastp result from cottongen
at<-read.table("flowering/utxAt.txt",head=T,sep="\t")
at$gene<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",at$Query)))
con$atg<-paste(at[match(con$V3,at$gene),"Match"],at[match(con$V4,at$gene),"Match"])

# append ZD's flowering reesults
x<-read.xls("flowering/UTX_flowering.xlsx")
x$gene<-gsub("Gohir_","Gohir.",gsub("[.].*","",gsub("Gohir.","Gohir_",x$id)))
table(x$gene %in% c(og$At,og$Dt))
# FALSE  TRUE
#   180   723 genes correspond to OG
table(og$At %in% x$gene | og$Dt %in% x$gene)
# FALSE  TRUE
# 22506   383 OGs are flowering related
con$AtFT <- x$short[match(con$V3,x$gene)]
con$DtFT <- x$short[match(con$V4,x$gene)]
dim(conFt<-con[!is.na(con$At)|!is.na(con$Dt),]) #383
write.table(conFt,"flowering/flow070522.txt",sep="\t")

select=which(conFt$cistrans%in%c("1.Cis only","2.Trans only","3.Cis+Trans: enhancing", "4.Cis+Trans: compensating","5.Compensatory"))
write.table(conFt[select,],"flowering/flow070522.sl.txt",sep="\t")
