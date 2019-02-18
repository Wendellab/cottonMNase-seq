####### DNS-seq pipline ########
## Modified from Daniel Vera's doc https://docs.google.com/document/d/18zrTm_VeCrHEKeWCtG7W0SGTKkYPtCkwuYJ9uNK7LZE/edit?ts=58e7a8c9

## Log in server
# ssh hugj2006@iastate.edu

## Required programs
# * R
# * cutadapt
# * samtools
# * bedtools
# * kent source utilities - "bedGraphToBigWig"

## load libraries and set threads
library(Rsamtools)
library(travis)
# options(threads=detectCores())
options(threads=8)
options(verbose=T)

## Start analysis from mapping results with quality >=20
bams=files("*.q20.bam")

## convert bam to fragment bed files
beds=bamToBed(bams,paired=T)

## examine fragment size distribution, error with bedHist
# bedHist(beds,xlims=c(0,200),dens=F,brks=100)
chrom=NULL
ss<-bedSizes(beds)
#rageHist(ss) error
pdf("sizeDistribution.pdf")
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


## Parse bed files by fragment size
bpl=bedParseLengths(beds,c(0,130,260))
bpl=unlist(bpl)

## calculate fragment densities.
allbeds=files("*.bed")
bgs=bedtoolsGenomeCov(allbeds,chromsizes)

## unify the bedgraphs so they have the same coordinates
ubgs=bgUnify(bgs,filler=0,discardUnshared=F)

## define light and heavy unified bedGraphs. Assumes light-digest files have an L_ in their names, and heavy-digest have H_
l=files("*L*unified.bg")
h=files("*H*unified.bg")

## check to see if light and heavy are paired properly
data.frame(l,h)

## calculate difference between light and heavy
dbgs=bgOps(l,"difference",h,pattern="L",replacement="D")


## convert bedGraph files to bigWig
bws=bedGraphToBigWig(c(l,h,dbgs),chromsizes)
# The resulting bigWig files can be used on most genome browsers

#quit(save="no")
q()
n