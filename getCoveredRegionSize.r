library(data.table)
fileL = grep("_q20_unified.bg$",list.files(),value=TRUE)
res = data.frame(chr = fileL,size=NA)
for(i in 1:length(fileL))
{
    rpm <- fread(fileL[i], select = 4)$V4;
    start <- fread(fileL[i], select = 2)$V2;
    end <- fread(fileL[i], select = 3)$V3;
    ts <- end-start
    cs <- sum(ts[rpm!=0])
    res$size[i] = cs
}
res
write.table(res,file="summaryCoveredGenomicSize.txt",row.names=FALSE, quote=FALSE,sep="\t")
