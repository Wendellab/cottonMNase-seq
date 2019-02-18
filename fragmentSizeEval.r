# set which species working on
genomes = c("AD1","F1","A2","D5")
dirs = c("mappingMnew","mappingF","mappingA_new","mappingD")
for(i in 1:4)
{
    setwd(paste0("/work/LAS/jfw-lab/hugj2006/cottonLeaf/",dirs[i]))
    genome = genomes[i]
    getwd()
    
    # use combo beds for simplisity
    if(genome=="A2"){beds = list.files(pattern="n_q20.bed")
    }else{ beds = list.files(pattern="c.*q20.bed") }
    
    # check range
    chrom=NULL
    library(travis)
    ss<-bedSizes(beds)
    pdf("plotFragmentSizeEval.pdf")
    plot(density(ss[[1]]), col=2, xlim=c(0,500), ylim=c(0,0.025),main = "Distribution of fragment size" )
    abline(v=100, col="grey")
    text(150,0.025,"100 bp",col="grey",cex=0.7)
    abline(v=142, col="grey")
    text(120,0.025,"142 bp",col="grey",cex=0.7)
    for(set in 2:length(ss))
    {
        lines(density(ss[[set]]),col=set+1 )
    }
    legend("topright", legend = paste0(names(ss),": ",lapply(ss,length)),col=2:5,lwd=3)
    dev.off()
}