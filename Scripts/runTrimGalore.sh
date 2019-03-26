#!usr/bin/bash
echo ''
echo ''
echo "Starting Job on "
stampStart=`date`
echo $stampStart 

SUFFIX=_1.fq.gz

for j in $( ls rawfastq/| grep _1.fq.gz$ );	do
   
echo ''
echo "==========Running Trim_galore for"
echo ${j%$SUFFIX}
echo ""
~/trim_galore_v0/trim_galore --paired -o trimmed/ rawfastq/$j rawfastq/${j%$SUFFIX}_2.fq.gz

done

echo ''
echo ''
echo "Ending  Job on "
stampEnd=`date`
echo $stampEnd
