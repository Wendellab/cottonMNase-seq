setwd("/work/LAS/jfw-lab/hugj2006/cottonLeaf/DNase/")
fqL <- c(grep("_1.fq.gz$",list.files("A2", full=T),value=T),grep("_1.fq.gz$",list.files("D5", full=T),value=T),grep("_1.fq.gz$",list.files("AD1", full=T),value=T))
###------Generate sequence and mapping statistics:
sumStat<-data.frame(sample=fqL, Type="DNase-seq", raw.reads=NA, filtered.reads=NA, mapping.input=NA, mapped=NA, mapped.Perc=NA, mapped.unique=NA, mapped.unique.Perc=NA, corr_estFragLen=NA, NSC=NA, RSC=NA,Qtag=NA)
rownames(sumStat) = gsub("_1.fq.gz","",fqL)
for(fq in fqL){
    each = gsub("_1.fq.gz","",fq)
    # from trimming report
    p1 <- readLines(paste0(each,"_1.fastq.gz_trimming_report.txt"))
    p1.raw.reads = gsub(".* ","",grep("Total reads",p1,value=TRUE))
    p1.filtered.reads = gsub(".*  | \\(.*","",grep("Reads written",p1,value=TRUE))
    p2 <- readLines(paste0(each,"_2.fastq.gz_trimming_report.txt"))
    p2.raw.reads = gsub(".* ","",grep("Total reads",p2,value=TRUE))
    p2.filtered.reads= gsub(".*  | \\(.*","",grep("Reads written",p2,value=TRUE))
    sumStat[each,"raw.reads"]=ifelse(p1.raw.reads==p2.raw.reads,p1.raw.reads,"Error: p1!=p2")
    sumStat[each,"filtered.reads"]=ifelse(p1.filtered.reads==p2.filtered.reads,p1.filtered.reads,"error: p1!=p2")
    # from mapping
    m <- readLines(paste0(each,".log"))
    input.reads <- as.numeric(gsub(" reads; of these:","",grep("reads",m,value=TRUE)))
    unmapped <- as.numeric(gsub(" |\\(.*","",grep("0 times",m,value=TRUE)))
    mapped <-input.reads-unmapped
    mappedUni <- as.numeric(gsub(" |\\(.*","",grep("exactly 1 time",m,value=TRUE)))
    sumStat[each,"mapping.input"]=input.reads
    sumStat[each,"mapped"]=mapped
    sumStat[each,"mapped.unique"]=mappedUni
    # from PCARD output
    p <- readLines(paste0(each,".dups.txt"))[8]
    ps<-strsplit(p,split="\t")
    sumStat[each,"duplicates"] =ps[[1]][7]
    sumStat[each,"optical.duplicates"] =ps[[1]][8]
    sumStat[each,"percent.duplication"] =ps[[1]][9]
    sumStat[each,"after.dup.removal"] =ps[[1]][10]
    # from phantompeakqual
    p <- read.table(paste0(each,".nodup.q20.phantomPeak.txt"),sep="\t")
    sumStat[each,"corr_estFragLen"] = p$V3
    sumStat[each,"NSC"] = p$V9
    sumStat[each,"RSC"] =p$V10
    sumStat[each,"Qtag"] =p$V11
}
sumStat$mapped.Perc=round(sumStat$mapped/sumStat$mapping.input*100,2)
sumStat$mapped.unique.Perc=round(sumStat$mapped.unique/sumStat$mapping.input*100,2)
sumStat
write.table(sumStat, "seqSumStats.txt", quote=FALSE, sep="\t")

###-------Report peak calling summary
library(data.table)
bedL <- list.files("D5",pattern="narrowPeak$|intersect.bed$|combine.bed$", full=T)
peakStat = data.frame(Bed=bedL)
rownames(peakStat)=bedL
for(each in bedL){
    df<-fread(each)
    peakStat[each,"peakN"] = nrow(df)
    peakStat[each,"peakBp"] = sum(df$V3-df$V2)
    w<-(df$V3-df$V2)
    peakStat[each,"minWidth"] = min(w)
    peakStat[each,"meanWidth"] = mean(w)
    peakStat[each,"medianWidth"] = median(w)
    peakStat[each,"maxWidth"] = max(w)
}
peakStat$bpPerGenome = peakStat$peakBp/749228090
peakStat
t=peakStat
#
bedL <- list.files("A2",pattern="narrowPeak$|intersect.bed$|combine.bed$", full=T)
peakStat = data.frame(Bed=bedL)
rownames(peakStat)=bedL
for(each in bedL){
    df<-fread(each)
    peakStat[each,"peakN"] = nrow(df)
    peakStat[each,"peakBp"] = sum(df$V3-df$V2)
    w<-(df$V3-df$V2)
    peakStat[each,"minWidth"] = min(w)
    peakStat[each,"meanWidth"] = mean(w)
    peakStat[each,"medianWidth"] = median(w)
    peakStat[each,"maxWidth"] = max(w)
}
peakStat$bpPerGenome = peakStat$peakBp/1509187826
peakStat
t= rbind(t,peakStat)
#
bedL <- list.files("AD1",pattern="narrowPeak$|intersect.bed$|combine.bed$", full=T)
peakStat = data.frame(Bed=bedL)
rownames(peakStat)=bedL
for(each in bedL){
    df<-fread(each)
    peakStat[each,"peakN"] = nrow(df)
    peakStat[each,"peakBp"] = sum(df$V3-df$V2)
    w<-(df$V3-df$V2)
    peakStat[each,"peakBp"] = sum(w)
    peakStat[each,"minWidth"] = min(w)
    peakStat[each,"meanWidth"] = mean(w)
    peakStat[each,"medianWidth"] = median(w)
    peakStat[each,"maxWidth"] = max(w)
}
peakStat$bpPerGenome = peakStat$peakBp/2281618630
t= rbind(t,peakStat)
write.table(t, "peakSumStats.txt", quote=FALSE, sep="\t",row.names=F)

q("no")
