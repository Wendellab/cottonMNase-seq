## FUN

# extract genomic regions from txDB
getPoint=function(gr,type=c("start","end","center"),be=1000, af=1000)
{
    grP =gr
    if(type=="start")
    {
        start(grP[strand(grP)=="+",]) = start(gr[strand(gr)=="+",])-be
        # end of the + strand genes must be equalized to start pos
        end(grP[strand(grP)=="+",])  = start(gr[strand(gr)=="+",])+af
        # start of the - strand genes must be equalized to end pos
        start(grP[strand(grP)=="-",]) = end(gr[strand(gr)=="-",])-af
        end(grP[strand(grP)=="-",]) = end(gr[strand(gr)=="-",])+be
    }
    if(type=="end")
    {
        # start of the + strand genes must be equalized to end pos
        start(grP[strand(grP)=="+",])  =end(gr[strand(gr)=="+",]) -be
        end(grP[strand(grP)=="+",]) = end(gr[strand(gr)=="+",]) + af
        # end of the - strand genes must be equalized to start pos
        start(grP[strand(grP)=="-",]) = start(gr[strand(gr)=="-",])-af
        end(grP[strand(grP)=="-",]) = start(gr[strand(gr)=="-",])+be
     }
    if(type=="center")
    {
        mid=round((start(grP)+end(grP))/2)
        # start of the + strand genes must be equalized to end pos
        start(grP[strand(grP)=="+",]) = mid[strand(grP)=="+"] - be
        end(  grP[strand(grP)=="+",]) = mid[strand(grP)=="+"] + af
        # end of the - strand genes must be equalized to start pos
        start(grP[strand(grP)=="-",]) = mid[strand(grP)=="-"] - af
        end(  grP[strand(grP)=="-",]) = mid[strand(grP)=="-"] + be
     }
    # remove duplicated
    grP=grP[!duplicated(grP),]
    seqlevels(grP)=unique(as.character(seqnames(grP)))
    return(grP)
}

getTSS=function(txdb, include=NULL, exclude.pattern=NULL, be=1000, af=1000){
    tss =  genes(txdb)
    if(!is.null(exclude.pattern)) {tss = tss[grep(exclude.pattern,tss$gene_id,invert=T)]}
    if(!is.null(include)) {tss = tss[tss$gene_id %in% include]}
    # +
    end(tss[strand(tss)=="+",])  =start(tss[strand(tss)=="+",])+af
    start(tss[strand(tss)=="+",])  =start(tss[strand(tss)=="+",])-be
    # -
    start(tss[strand(tss)=="-",])=end(tss[strand(tss)=="-",])-af
    end(tss[strand(tss)=="-",])=end(tss[strand(tss)=="-",])+be
    # remove duplicated TSSes ie alternative transcripts
    # this keeps the first instance and removes duplicates
    tss=tss[!duplicated(tss),]
    seqlevels(tss)=unique(as.character(seqnames(tss)))
    return(tss)
}

getTTS=function(txdb, include=NULL, exclude.pattern=NULL, be=1000, af=1000){
    tts =  genes(txdb)
    if(!is.null(exclude.pattern)) {tts = tts[grep(exclude.pattern,tts$gene_id,invert=T)]}
    if(!is.null(include)) {tts = tts[tts$gene_id %in% include]}
    # +
    end(tss[strand(tss)=="+",])   = end(tss[strand(tss)=="+",])+af
    start(tss[strand(tss)=="+",]) = end(tss[strand(tss)=="+",])-be
    # -
    end(tss[strand(tss)=="-",])   = start(tss[strand(tss)=="-",])+be
    start(tss[strand(tss)=="-",]) = start(tss[strand(tss)=="-",])-af
    # remove duplicated TSSes ie alternative transcripts
    # this keeps the first instance and removes duplicates
    tts=tts[!duplicated(tts),]
    seqlevels(tts)=unique(as.character(seqnames(tts)))
    return(tts)
}

plotTSSnTTS=function(bw.files,tss.gr,tts.gr,tag,title)
{
    require(genomation)
    #TSS
    sml = ScoreMatrixList(target = bw.files, windows = tss.gr, bin.num=50, type="bigWig",strand.aware=TRUE)
    names(sml) = tag
    plotMeta(sml,  xcoords = c(-3000, 1000), profile.names=names(sml), ylab="RPKM",dispersion = "se", xlab="TSS", main=title)
    abline(v=0,lty=2)
    #TTS
    sml = ScoreMatrixList(target = bw.files, windows = tts.gr, bin.num=50, type="bigWig",strand.aware=TRUE)
    names(sml) = tag
    plotMeta(sml,  xcoords = c(-1000, 3000), profile.names=names(sml), ylab="RPKM",dispersion = "se", xlab="TTS", main=title)
    abline(v=0,lty=2)
}

