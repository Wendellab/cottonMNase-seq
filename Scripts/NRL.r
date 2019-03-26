# Nucleosome repeat length was calculated using NucTools
# https://homeveg.github.io/nuctools/

## installation
cd /work/LAS/jfw-lab/hugj2006/
git clone https://github.com/homeveg/nuctools.git NucTools
ls NucTools/

## prep input bed files: parsing into subgenomes, then parsing chromosomes
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/nucleosome
mkdir NRL; cd NRL
grep '^A' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingF/FcH_q20.bed >FcH_q20.A.bed
grep '^D' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingF/FcH_q20.bed >FcH_q20.D.bed
grep '^A' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingMnew/McH_q20.bed >McH_q20.A.bed
grep '^D' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingMnew/McH_q20.bed >McH_q20.D.bed
awk -F"\t" '{print > "McH_q20.A" $1 ".bed"}' McH_q20.A.bed
awk -F"\t" '{print > "McH_q20.D" $1 ".bed"}' McH_q20.D.bed
awk -F"\t" '{print > "FcH_q20.A" $1 ".bed"}' FcH_q20.A.bed
awk -F"\t" '{print > "FcH_q20.D" $1 ".bed"}' FcH_q20.D.bed
ln -s /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingD/DcH_q20.bed
awk -F"\t" '{print > "DcH_q20." $1 ".bed"}' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingD/DcH_q20.bed
ln -s /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingA_new/A6Hn_q20.bed
awk -F"\t" '{print > "A6Hn_q20." $1 ".bed"}' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingA_new/A6Hn_q20.bed

## use NucTools to caculate NRL
perl /work/LAS/jfw-lab/hugj2006/NucTools/nucleosome_repeat_length.pl -in DcH_q20.bed -out NRL.txt
/work/LAS/jfw-lab/hugj2006/NucTools/misc/plotNRL.R --input=NRL.txt --dir=. --sample=DcH_q20


# ------------------ R code below ------------------------
fileL=list.files(pattern="bed$")
for(file in fileL){
    cmd=paste0("perl /work/LAS/jfw-lab/hugj2006/NucTools/nucleosome_repeat_length.pl -in ",file," -out ", gsub(".bed","",file),".NRL.txt; /work/LAS/jfw-lab/hugj2006/NucTools/misc/plotNRL.R --input=NRL.txt --dir=. --sample=",gsub(".bed","",file))
    message(cmd)
    system(cmd)
}
# about 16 hours


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
