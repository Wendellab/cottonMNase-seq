# https://github.com/gthar/nucleServ/blob/master/bin/nucleR.R

library(data.table)
library(GenomicRanges)
library(nucleR)
library(travis)
library(rtracklayer)
thread=8

# what portions of genomes are covered by numcleosomes? Well-positioned or fuzzy?
res=c("sample","genome", "N.Fuzzy","N.uncertain","N.well","N.well150", "bp.Fuzzy","bp.uncertain","bp.well","bp.well150")
for(sample in c("McH","FcH","A6Hn","DcH")){
    print(sample)
    df=fread(paste0(sample,"_q20.nucleosome.gff"),select=c("V1","V9"))
    width =as.numeric(gsub("width=|;.*","",df$V9))
    class =gsub(".*=","",df$V9)
    w150=which(width<=150 & class=="W")
    stbl = c(sample, substr(sample,1,1),table(class), table(class[w150]))
    # covered genomic regions
    f= sum(width[class=="F"])
    u= sum(width[class=="uncertain"])
    w= sum(width[class=="W"])
    ws= sum(width[w150])
    stbl=c(stbl, f,u,w,ws)
    res=rbind(res,stbl)
    if(sample =="McH"|sample=="FcH"){
        genome =gsub("[0-9]","",df$V1)
        # nuc numbers
        tbl = cbind(rep(sample,2),c("A","D"),t(xtabs(~class+genome)), table(genome[w150]))
        # covered genomic regions
        f= aggregate(width[class=="F"],list(genome[class=="F"]),sum)
        u= aggregate(width[class=="uncertain"],list(genome[class=="uncertain"]),sum)
        w= aggregate(width[class=="W"],list(genome[class=="W"]),sum)
        ws= aggregate(width[w150],list(genome[w150]),sum)
        tbl=cbind(tbl, f$x,u$x,w$x,ws$x)
        res = rbind(res,tbl)
    }
}
res
write.table(res,file="sumTbl.nuclosome.txt",row.names=F,col.names=F)
