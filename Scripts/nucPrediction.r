## ------------ Nucleosome positioning prediction ---------------
# A summary of tools
# https://generegulation.org/nucleosome_positioning_prediction/
## --------------------------------------------

## NuPop, Xi etal 2010 and Wang et al 2008
# BiocManager::install("NuPoP")
library(NuPoP)
library(data.table)

# Species 10 = A. thaliana; 11 = Maize; 0 = other, will identify a most similiar species from 1-11 according to sequence base composition
# Model = 4 or 1. 4 costs 2.5 times of 1, slower.
fileL = list.files(pattern="fasta$")
sumT = data.frame(fileL)
rownames(sumT)=fileL
for (each in fileL){
    predNuPoP(each,species=0,model=1)
    # summarize
    # res=readNuPop(paste0(each,"_Prediction1.txt"))
    res=fread(paste0(each,"_Prediction1.txt"))
    sumT[each,"chrLen"] = nrow(res)
    sumT[each,"NucLen"] = sum(res$"N/L")
    sumT[each,"NucCov"] = sumT[each,"NucLen"]/sumT[each,"chrLen"]
    #nrl
    res$incre = c(0,res$"N/L"[2:nrow(res)] - res$"N/L"[1:(nrow(res)-1)] )
    # fixed nucleosome length, start and end the same
    useS= res$Position[res$incr==1]
    #useE= res$Position[res$incr==-1]
    nrlS = useS[2:length(useS)] - useS[1:(length(useS)-1)]
    #nrlE = useE[2:length(useE)] - useE[1:(length(useE)-1)]
    nrl=nrlS[nrlS<quantile(nrlS,0.90)]
    sumT[each,"NRL.mean"] = mean(nrl)
    sumT[each,"NRL.median"] = median(nrl)
    sumT[each,"NRL.sd"] = sd(nrl)
    print(sumT[each,])
    system("rm *Prediction1.txt")
}
write.table(sumT,file="predict.txt",sep="\t",row.names=F)
q("no")

----------Bash control

cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/Nucleosome/

awk -F "\t" '/^>/ {F = $1".fasta"} {print > F}' ../refGenomes/Dgenome2_13.fasta
R CMD BATCH nucPrediction.r nucPrediction.D.Rout.txt
mv predict.txt predict.D.txt
rm *fasta

awk -F "\t" '/^>/ {F = $1".fasta"} {print > F}' ../refGenomes/A2WHU_13.fasta
R CMD BATCH nucPrediction.r nucPrediction.A.Rout.txt
mv predict.txt predict.A.txt
rm *fasta

awk -F "\t" '/^>/ {F = $1".fasta"} {print > F}' ../refGenomes/TM1utx_26.fasta
R CMD BATCH nucPrediction.r nucPrediction.M.Rout.txt
mv predict.txt predict.M.txt
rm *fasta
