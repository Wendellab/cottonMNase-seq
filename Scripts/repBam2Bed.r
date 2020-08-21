####### Convert BAM to BED and pool replicates ########

## load libraries and set threads
library(travis)
# options(threads=detectCores())
options(threads=8)
options(verbose=T)

## Start analysis from mapping results with quality >=20
bams=list.files(pattern="q20.bam$")
bams

## convert bam to fragment bed files
for(i in bams){
    cmd=paste0("samtools depth ",i," | awk '{if($3>0) total+=1}END{print total}'")
    message(cmd)
    system(cmd) # get genomic bp covered by reads
    cmd <- paste0("samtools sort -n ",i," | bedtools bamtobed -i - -bedpe  | cut -f 1,2,6,7,8,9 | sort -T . -k1,1 -k2,2n -S 1G  > ",gsub("bam","bed",i))
    message(cmd)
    system(cmd)
}
beds=list.files(pattern="q20.bed$")
beds
# below code works for name sorted bams
# beds=bamToBed(bams,paired=T)

## examine fragment size distribution, error with bedHist
# bedHist(beds,xlims=c(0,200),dens=F,brks=100)
chrom=NULL
ss<-bedSizes(beds)
#rageHist(ss) error
pdf("plotSizeDistribution.pdf")
plot(density(ss[[1]]), col=2, xlim=c(0,500), ylim=c(0,0.025),main = "Distribution of fragment size" )
abline(v=147, col="grey")
text(150,0.025,"147 bp",col="grey",cex=0.7)
abline(v=130, col="grey")
text(120,0.025,"130 bp",col="grey",cex=0.7)
for(set in 2:length(ss))
{
    lines(density(ss[[set]]),col=set+1 )
}
legend("topright", legend = paste0(names(ss),": ",lapply(ss,length)),col=2:5,lwd=3)
dev.off()

## pool replicates
len=unlist(lapply(ss,length))
con =ifelse(grepl("L",names(ss)),"L","H")
size = tapply(len,factor(con),min)
df = data.frame(len,con,size[con])
df$subsample = df$len>df$size.con.
df$cmd =''
for(i in 1:nrow(df))
{
    if(df$subsample[i]){
        df$cmd[i] = paste0("<(shuf -n ",df$size.con.[i]," ",rownames(df)[i],")")
    }else{
        df$cmd[i] = rownames(df)[i]
    }
}
df$out = gsub("\n","c",rownames(df))
# inspect df to make sure
df
write.table(df,file="poolReps.txt",sep="\t")
# run for conditions
genome=substring(beds[1],1,1)
cat("# pool rep",file="poolReps.sh",sep="\n")
for(each in c("H","L"))
{
    cmd=paste0("cat ", paste(df$cmd[df$con==each], collapse=" ")," | sort -T . -k1,1 -k2,2n -S 1G  >",genome,"c",each,"_q20.bed")
    cat(cmd,file="poolReps.sh",append=T,sep="\n")
    message(cmd)
}
system("bash poolReps.sh")
# will have to paste and run command in bash
# cat <(shuf -n 22277062 D1H_q20.bed) D2H_q20.bed | sort -T . -k1,1 -k2,2n -S 1G  >DcH_q20.bed
# cat D1L_q20.bed <(shuf -n 18024273 D2L_q20.bed) | sort -T . -k1,1 -k2,2n -S 1G  >DcL_q20.bed
# cat <(shuf -n 113082817 M1H_q20.bed) M2H_q20.bed | sort -T . -k1,1 -k2,2n -S 1G  >McH_q20.bed
# cat M1L_q20.bed <(shuf -n 77827484 M2L_q20.bed) | sort -T . -k1,1 -k2,2n -S 1G  >McL_q20.bed
# cat <(shuf -n 76347017 F2H_q20.bed) F3H_q20.bed | sort -T . -k1,1 -k2,2n -S 1G  >FcH_q20.bed
# cat F2L_q20.bed <(shuf -n 72395403 F3L_q20.bed) | sort -T . -k1,1 -k2,2n -S 1G  >FcL_q20.bed

q("no")
