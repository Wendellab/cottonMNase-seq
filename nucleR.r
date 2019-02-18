## ------------ NucleR analysis ---------------
# http://bioconductor.org/packages/release/bioc/vignettes/nucleR/inst/doc/nucleR.pdf
# 1. data preprocessing: SE needs fragment length to correct strand-specific mapping that shifts reads, PE do not shift reads but filter out longer fragments; trim narrows fragments to nucleosome core, improving peak calling.
# 2. filtering noise for coverage profile using Fast Fourier Transform (FFT)
# 3. peak detection for nucleosome positioning: 25% threshold checked first quantile data to detect peaks
# 4. define well arranged and loosely positioned nucleosomes
# 5. export peaks and coverage
## --------------------------------------------

# https://github.com/gthar/nucleServ/blob/master/bin/nucleR.R

library(data.table)
library(GenomicRanges)
library(nucleR)
library(travis)
library(rtracklayer)
thread=8

#FUN
df2gff <- function (df, ...)
{   # Convert a data.frame into a form that is easily saved in disk as a gff.
    # It tries to match gff fields to the row names of the data.frame.
    # gff fields not found in the data.frame will be set to `.` and fields in
    # the data.frame not found in the gff fields will be stored accordingly in
    # the `attribute` field, separated by `;`.
    # It also accepts optional keyword arguments to specify gff fields that
    # are not present in the data.frame.
    kwargs <- list(...)
    fields <- c("seqnames", "source", "feature", "start", "end", "score",
    "strand", "frame")
    out.df <- data.frame(matrix(nrow=nrow(df),
    ncol=length(fields) + 1,
    dimnames=list(c(),
    c(fields, "attribute"))))
    for (f in fields) {
        if (f %in% colnames(df)) {
            out.df[[f]] <- as.vector(df[[f]])
        } else if (f %in% names(kwargs)) {
            if (length(kwargs[[f]]) == 1) {
                out.df[[f]] <- rep(kwargs[[f]], nrow(out.df))
            } else {
                out.df[[f]] <- kwargs[[f]]
            }
        } else {
            out.df[[f]] <- rep(".", nrow(out.df))
        }
    }
    nonfield.columns <- colnames(df)[!colnames(df) %in% fields]
    attrVal <- function(i, df) sprintf("%s=%s", i, df[[i]])
    attrs <- do.call(paste,
    c(lapply(nonfield.columns,
    attrVal,
    df),
    sep=";"))
    if (length(attrs)) {
        out.df[["attribute"]] <- attrs
    } else {
        out.df[["attribute"]] <- rep(".", nrow(out.df))
    }
    out.df
}

writeGff <- function (df, outpath)
{
    # Use this to write the output of df2gff to disk.
    write.table(df,
    file=outpath,
    quote=FALSE,
    sep="\t",
    row.names=FALSE,
    col.names=FALSE)
}


#
#source("nucleR.source.r")

## load bed files, loop over files
# prepare loop
genomes =c("AD1","F1","A2","D5")
tags=c("M","F","A","D")
dirs = c("mappingMnew","mappingF","mappingA_new","mappingD")
txdbs=c("refGenomes/txdb.TM1saski.sqlite","refGenomes/txdb.F1.sqlite","refGenomes/txdb.A2du.sqlite","refGenomes/txdb.D5.sqlite")
outdir="nucleosome/"
for(i in 1:4)
{
    genome = genomes[i]
    dir=dirs[i]
    # work on light digestion
    bedfile=list.files(dir,pattern="(cH|6Hn)_q20.bed",full.name=T)
    print(bedfile)
    
    message("...read bed files into data frame")
    df<- fread(bedfile, data.table=FALSE, header=FALSE, select=c(1,2,3,6))
    names(df)<-c("seqnames","start","end","strand")
    # import BED file is 0 based
    df$start=df$start+1
    
    message("...convert data frame to genome ranges")
    # rd<- as(df, "RangedData")
    gr <- makeGRangesFromDataFrame(df, keep.extra.columns=TRUE, start.field="start", end.field="end", starts.in.df.are.0based=FALSE)
    # nuckeR requires "RangedData" from "IRanges" package, or "AlignedRead" object containing the position, length and strand of the sequence reads.
    sl<-read.table(paste0(dir,"/chr.size.txt"))
    seqlengths(gr) = sl$V2
    
    message("...processing reads")
    # -- Dyad Length = 50, trim each read to 50bp around the dyad, use reads <200bp long
    prep <- processReads(gr, type="paired", fragmentLen=200, trim=50)
    
    message("...calculate the coverage, directly in reads per million (r.p.m)")
    cover <- coverage.rpm(prep)
    
    # message("...plot coverage summary by chromosome")
    # check raw depth, are there enough coverage from sequencing?
    # depth<-coverage(rd)
    # pdf(paste0(gsub(".bed","",bedfile),".checkRawDepth.pdf"))
    # par(mfrow=c(1,2))
    # for(i in 1:length(depth)){
    #    chr<-names(depth)[i]
    #    x<-as.data.frame(table(depth[[i]]))
    #    x$cumsum<- cumsum(x$Freq)
    #    x$cumfrac<- x$cumsum/sl$V2[sl$V1==chr]
    #    plot(x$Var1[-1], 1-x$cumfrac[-1], type='n', xlab="coverage", ylab="fraction of bases sampled >= coverage", ylim=c(0,1.0), xlim=c(1,10), main=paste0(chr,": Target Region Coverage"))
    #    plot(x$Var1[-1], x$Freq[-1]/sl$V2[sl$V1==chr], type='n', xlab="coverage", ylab="frequence of bases sampled >= coverage", xlim=c(1,10), main=paste0(chr,": Target Region Coverage"))
    #}
    #dev.off()
    
    message("...write coverage into WIG")
    export.bw(cover,paste0(outdir,gsub(".*/|.bed","",bedfile),".nucleosome.coverage.bw"))
    # export.wig(cover, paste0(gsub(".bed","",bedfile),".coverage"))
    
    message("...filtering noise for coverage profile")
    fft<- filterFFT(cover,pcKeepComp = 0.02, mc.preschedule = FALSE, mc.cores = thread)
    # emptyHandler <- function (f)
    # function (x, ...)
    # if (length(x) > 0) {
    #    f(x, ...)
    # } else {
    #    numeric()
    #}
    # 2% component is suggested for NGS
    #fft <- mclapply(cover, emptyHandler(filterFFT), pcKeepComp = 0.02, mc.preschedule = FALSE, mc.cores = thread)
    
    message("...peak detection for nucleosome positioning")
    peaks <- peakDetection(lapply(fft,function(x)x+2), width  =147, threshold = "35%", score = TRUE, mc.cores  = thread)
    # min.cov=2, +2 for fft
    # width =147, Size of each nucleosome, in bp, to be considered
    # threshold = "35%",: Minimum number of reads (Coverage) to call a nucleosome. Can be given as a percentage (i.e., "25%" means that the peaks with coverage in the 1st quartile of data won't be considered); or as an absolute coverage value (i.e., "20" means that the peaks with less than 20 reads per million of mapped reads won't be considered). Default = 35%.
    # score = TRUE, same as
    # scores <- peakScoring(peaks, fft, threshold = "35%", dyad.length = 50, mc.cores = thread)
    # peak width score (score_w) indicates positioning, higher means sharper peak
    # peak height score (score_h) indicates coverage, higher means talled peak
    
    message("...merging peaks")
    merged <- mergeCalls(peaks, min.overlap = 50, mc.cores  = thread)
    # min.overlap - Minimum overlap between two reads for merge them
    # discard.low=0.2 - Discard low covered calls (i.e. calls with score_h < discard.low will be discarded)
    
    message("...peak classifications")
    # Well-positioned (W): score_w>0.4 & score_h>0.6; Otherwise fuzzy (F)
    getType <- function(score_w, score_h, thresh_w, thresh_h, nmerge=1)
    {
        ifelse(nmerge > 2, "uncertain", ifelse(`&`(score_w > thresh_w, score_h > thresh_h),"W", "F"))
    }
    merged$class <- getType(merged$score_w, merged$score_h, 0.4, 0.6, merged$nmerge)
    
    message("...summarize nucleosome classes")
    print(table(merged$class ))
    
    message("...write nucleosome positions into GFF")
    merged <- as.data.frame(merged)
    names(merged)[names(merged) == "space"] <- "seqname"
    names(merged)[names(merged) == "score_h"] <- "score_height"
    names(merged)[names(merged) == "score_w"] <- "score_width"
    merged$nmerge <- NULL
    writeGff(df2gff(merged, source  = "nucleR", feature = "Nucleosome"), paste0(outdir,gsub(".*/|.bed","",bedfile),".nucleosome.gff") )
    
}

# look at nucleosomes
library(data.table)
# F1
df=fread("FcH_q20.nucleosome.gff",select=c("V1","V9"))
width =as.numeric(gsub("width=|;.*","",df$V9))
class =gsub(".*=","",df$V9)
genome =gsub("[0-9]","",df$V1)
xtabs(~class+genome)
# class             A       D
# F         1474311 1061655
# uncertain 1193338  637100
# W         2352929 1397439
aggregate(width,by=list(genome,class),sum)
# Group.1   Group.2         x
# 1       A         F 252765022
# 2       D         F 180291475
# 3       A uncertain 467251022
# 4       D uncertain 243969318
# 5       A         W 404037780
# 6       D         W 237655721
sum(width)
# 1785970338
sum(width[width<=150& class=="W"])
# 376447890
sum(width[width<=150& class=="W"& genome=="A"])
# 233352357
sum(width[width<=150& class=="W"& genome=="D"])
# 143095533

# AD1
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

q()
n
#############################################################################################
## deeptools plot nuclesome around TSS
cd nucleosome/
# AD1
computeMatrix reference-point -S Mc*bw -R ../expNuc/AD1.D*bed -o tss.gz --referencePoint TSS -b 1000 -a 1000 --skipZeros
plotHeatmap -m tss.gz -out plotHeatmap.quantile.AD1.D.png
plotProfile -m tss.gz -out plotProfileGroup.AD1.D.png
plotProfile -m tss.gz -out plotProfileGroup.AD1.D.group.png --perGroup

computeMatrix reference-point -S Mc*bw -R ../expNuc/AD1.D*bed -o tes.gz --referencePoint TES -b 1000 -a 1000 --skipZeros
plotHeatmap -m tes.gz -out plotHeatmap.quantile.AD1.D.png
plotProfile -m tes.gz -out plotProfileTES.AD1.D.png

computeMatrix reference-point -S ../mappingMnew/McD_q20_chr.size_w20_qnorm.bw -R ../expNuc/AD1.D*bed -o tss.gz --referencePoint TSS -b 1000 -a 1000 --skipZeros
plotProfile -m tss.gz -out plotProfile.dns.AD1.D.png

#############################################################################################
## check size range 0-141 vs 142-200
bedToBam -i mappingF/FcH_q20.bed -g mappingF/chr.size.txt > mappingF/FcH_q20.bam
samtools index mappingF/FcH_q20.bam
bamCoverage -b mappingF/FcH_q20.bam -o nucleosome/FcH_q20_w100k_0-200.bw --binSize 100000 --minFragmentLength 0 --maxFragmentLength 200
bamCoverage -b mappingF/FcH_q20.bam -o nucleosome/FcH_q20_w100k_142-200.bw --binSize 100000 --minFragmentLength 142 --maxFragmentLength 200
bamCoverage -b mappingF/FcH_q20.bam -o nucleosome/FcH_q20_w100k_0-141.bw --binSize 100000 --minFragmentLength 0 --maxFragmentLength 141
bigwigCompare -b1 nucleosome/FcH_q20_w100k_142-200.bw -b2 nucleosome/FcH_q20_w100k_0-200.bw --binSize 100000 -o nucleosome/FcH_q20_w100k_142-200.log2ratio.bw
bigwigCompare -b1 nucleosome/FcH_q20_w100k_0-141.bw -b2 nucleosome/FcH_q20_w100k_0-200.bw --binSize 100000 -o nucleosome/FcH_q20_w100k_0-141.log2ratio.bw
#start r
library(rtracklayer)
grL<- import.bw("~/Downloads/FcH_q20_w100k_142-200.log2ratio.bw")
grS<- import.bw("~/Downloads/FcH_q20_w100k_0-141.log2ratio.bw")
use=which(seqnames(grL)=="A01")
plot(end(grL[use]),grL[use]$score, xlab = "A01", col = 'blue', type="l", ylim=c(-2,0))
lines(end(grS[use]),grS[use]$score,type="l", col = 'red')


bamCoverage -b F2H_q20.bam -o F2H_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 1878118078
computeMatrix scale-regions -S *-*bw -R ../expNuc/F1.A.Q*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out plotSizediff.A.png
computeMatrix scale-regions -S *-*bw -R ../expNuc/F1.D.Q*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out plotSizediff.D.png
# F2L 142-200 77.7%, 0-141 21.2%, L should have lower portion of small fragment than H
# F2H 142-200 68.1%, 0-141 30.6%
~/kent/bigWigToBedGraph F2H_q20_0-141.bw F2H_q20_0-141.bg
head F2H_q20_0-141.bg


cd mappingMnew
bamCoverage -b M1L_q20.bam -o M1L_q20_142-200.bw --binSize 20 --minFragmentLength 142 --maxFragmentLength 200 --normalizeTo1x 1840638571
bamCoverage -b M1L_q20.bam -o M1L_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 1840638571
bamCoverage -b M1H_q20.bam -o M1H_q20_142-200.bw --binSize 20 --minFragmentLength 142 --maxFragmentLength 200 --normalizeTo1x 1840638571
bamCoverage -b M1H_q20.bam -o M1H_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 1840638571
computeMatrix scale-regions -S *-*bw -R ../expNuc/AD1.A.Q*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out plotSizediff.A.png
computeMatrix scale-regions -S *-*bw -R ../expNuc/AD1.D.Q*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out plotSizediff.D.png
# M1L 142-200 72.2%, 0-141 26.8%, L should have lower portion of small fragment than H
# M1H 142-200 50.4%, 0-141 49.6%

cd ../mappingA_new/
bamCoverage -b A6Ln_q20.bam -o A6Ln_q20_142-200.bw --binSize 20 --minFragmentLength 142 --maxFragmentLength 200 --normalizeTo1x 1167736467
bamCoverage -b A6Ln_q20.bam -o A6Ln_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 1167736467
bamCoverage -b A6Hn_q20.bam -o A6Hn_q20_142-200.bw --binSize 20 --minFragmentLength 142 --maxFragmentLength 200 --normalizeTo1x 1167736467
bamCoverage -b A6Hn_q20.bam -o A6Hn_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 1167736467
computeMatrix scale-regions -S *-*bw -R ../expNuc/A2.Q*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out plotSizediff.png
# A6Ln 142-200 66.4%, 0-141 30.4%, L should have lower portion of small fragment than H
# A6Hn 142-200 53.6%, 0-141 44.9%

cd ../mappingD/
bamCoverage -b D1L_q20.bam -o D1L_q20_142-200.bw --binSize 20 --minFragmentLength 142 --maxFragmentLength 200 --normalizeTo1x 716593377
bamCoverage -b D1L_q20.bam -o D1L_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 716593377
bamCoverage -b D1H_q20.bam -o D1H_q20_142-200.bw --binSize 20 --minFragmentLength 142 --maxFragmentLength 200 --normalizeTo1x 716593377
bamCoverage -b D1H_q20.bam -o D1H_q20_0-141.bw --binSize 20 --minFragmentLength 0 --maxFragmentLength 141 --normalizeTo1x 716593377
computeMatrix scale-regions -S *-*bw -R ../expNuc/D5.Q*bed -o scaleGeneBody.gz --regionBodyLength 3000 -b 1000 -a 1000 --skipZeros
plotHeatmap -m scaleGeneBody.gz -out plotSizediff.png
# D1L 142-200 66.7%, 0-141 32.5%, L should have lower portion of small fragment than H
# D1H 142-200 58.7%, 0-141 40.5%
# IGV shows no apparent difference in read distribution of different size range



#############################################################################################
## use NucTools to caculate NRL
/work/LAS/jfw-lab/hugj2006/NucTools
grep '^A' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingF/FcH_q20.bed >FcH_q20.A.bed
grep '^D' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingF/FcH_q20.bed >FcH_q20.D.bed
grep '^A' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingMnew/McH_q20.bed >McH_q20.A.bed
grep '^D' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingMnew/McH_q20.bed >McH_q20.D.bed
awk -F"\t" '{print > "McH_q20.A" $1 ".bed"}' McH_q20.A.bed
awk -F"\t" '{print > "McH_q20.D" $1 ".bed"}' McH_q20.D.bed
awk -F"\t" '{print > "FcH_q20.A" $1 ".bed"}' FcH_q20.A.bed
awk -F"\t" '{print > "FcH_q20.D" $1 ".bed"}' FcH_q20.D.bed
awk -F"\t" '{print > "DcH_q20." $1 ".bed"}' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingD/DcH_q20.bed
awk -F"\t" '{print > "A6Hn_q20." $1 ".bed"}' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingA_new/A6Hn_q20.bed
perl nucleosome_repeat_length.pl -in /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingD/DcH_q20.bed -out DcH.NRL.txt
./misc/plotNRL.R --input=DcH.NRL.txt --dir=. --sample=DcH
perl nucleosome_repeat_length.pl -in /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingA_new/A6Hn_q20.bed -out A6Hn.NRL.txt
./misc/plotNRL.R --input=A6Hn.NRL.txt --dir=. --sample=A6Hn
perl nucleosome_repeat_length.pl -in McH_q20.A.bed -out McHA.NRL.txt
./misc/plotNRL.R --input=McHA.NRL.txt --dir=. --sample=McHA
perl nucleosome_repeat_length.pl -in McH_q20.D.bed -out McHD.NRL.txt
./misc/plotNRL.R --input=McHD.NRL.txt --dir=. --sample=McHD

# R code
fileL=list.files(pattern="bed$")
for(file in fileL){
    cmd=paste0("perl nucleosome_repeat_length.pl -in ",file," -out NRL.txt; ./misc/plotNRL.R --input=NRL.txt --dir=. --sample=",gsub(".bed","",file))
    message(cmd)
    system(cmd)
}
NRL=data.frame(A2=c(198, 198, 197, 196, 198, 198, 196, 197, 127, 195, 196, 198,196, 198), D5=c(194,198,194,197,196,194,197,198,199,194,194,196,198, 194), F1.A=c(199,198,198,198,199,197, 197,197,150, 198,197,197, 197,199), F1.D=c(198,195,197,197,193,194,200,195,198,198,148,196,196,198), AD1.A=c(201, 199,201,197,201,198,199,199,196,200, 201,197,200,201),AD1.D=c(201,196,202,198,201,152, 200,198,193,202,130,201,201 ,201)); NRL
apply(NRL,2,mean)
sem<-sd(x)/sqrt(length(x))
apply(NRL,2,sem)
apply(NRL,2,sd)
# why some chrs are so packed??
# need to modifly "plotNRL.R" to re do all chrs, get SEs
# histone proteins
annot = read.table("~/Dropbox/MNase-seq/cottonMNase-seq_git/cotton.v2.1.mbh.arabi.txt",sep="\t", header=TRUE)
library(gdata)
histone = read.xls("~/Dropbox/MNase-seq/cottonMNase-seq_git/Gorai.key.histone.xlsx")
core = histone[grep("PF00125", histone$PFAM),]
id=as.character(core$locusName)
t.test(rowMeans(A2$rpkm[id,]), rowMeans(D5$rpkm[id,]) )
