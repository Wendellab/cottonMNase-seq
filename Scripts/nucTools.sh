# Nucleosome repeat length was calculated using NucTools
# https://homeveg.github.io/nuctools/

## installation
cd /work/LAS/jfw-lab/hugj2006/
git clone https://github.com/homeveg/nuctools.git NucTools
ls NucTools/

cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/Nucleosome
mkdir NRL

## prep input bed files: parsing into subgenomes, then parsing chromosomes
cd /work/LAS/jfw-lab/hugj2006/cottonLeaf/nucleosomeNRL
grep '^A' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingF2020/FcH_q20.bed >FcH_q20.A.bed
grep '^D' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingF2020/FcH_q20.bed >FcH_q20.D.bed
grep '^A' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingM_UTX/McH_q20.bed >McH_q20.A.bed
grep '^D' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingM_UTX/McH_q20.bed >McH_q20.D.bed
awk -F"\t" '{print > "McH_q20.A" $1 ".bed"}' McH_q20.A.bed
awk -F"\t" '{print > "McH_q20.D" $1 ".bed"}' McH_q20.D.bed
awk -F"\t" '{print > "FcH_q20.A" $1 ".bed"}' FcH_q20.A.bed
awk -F"\t" '{print > "FcH_q20.D" $1 ".bed"}' FcH_q20.D.bed
ln -s /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingD/DcH_q20.bed
awk -F"\t" '{print > "DcH_q20." $1 ".bed"}' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingD/DcH_q20.bed
ln -s /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingA_new/A6Hn_q20.bed
awk -F"\t" '{print > "A6Hn_q20." $1 ".bed"}' /work/LAS/jfw-lab/hugj2006/cottonLeaf/mappingA_WHU/A6Hn_q20.bed

## use NucTools to caculate NRL
for bed in $(ls *bed); do
   echo ${bed%%.bed}
   perl /work/LAS/jfw-lab/hugj2006/NucTools/nucleosome_repeat_length.pl -in $bed -out NRL.txt
   /work/LAS/jfw-lab/hugj2006/NucTools/misc/plotNRL.R --input=NRL.txt --dir=. --sample=${bed%%.bed}
done
# about 16 hours
