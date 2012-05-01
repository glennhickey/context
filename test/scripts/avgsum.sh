#!/bin/bash

outPathBase=${1}

len=${#outPathBase}
len=$(($len - 1))
lchar=${outPathBase:${len}}    
if [ $lchar == "/" ]
then
    outPathBase=`echo ${outPathBase} |sed "s/.$//"`
fi

pfile=${outPathBase}/params.txt

name=${outPathBase##*/}

mu=`tail -1 ${pfile} | awk '$1 == "<*t=" {sum += $4; r+=1} END {print sum / r}'`

lenAvg=`cat ${outPathBase}/diffsum.txt | awk 'NR > 1 {sum += $2} END {print sum / (NR - 1)}'`

lenSD=`cat ${outPathBase}/diffsum.txt | awk -v x=${lenAvg} 'NR > 1 {delta = $2 - x; sum += (delta >= 0 ? delta : -delta)} END {print sum / (2 * (NR - 1))}'`

diffAvg=`cat ${outPathBase}/diffsum.txt | awk 'NR > 1 {sum += $3} END {print sum / (NR - 1)}'`

diffSD=`cat ${outPathBase}/diffsum.txt | awk -v x=${diffAvg} 'NR > 1 {delta = $3 - x; sum += (delta >= 0 ? delta : -delta)} END {print sum / (2 * (NR - 1))}'`

diffAvgNC=`cat ${outPathBase}/diffsum.txt | awk 'NR > 1 {sum += $5} END {print sum / (NR - 1)}'`

diffSDNC=`cat ${outPathBase}/diffsum.txt | awk -v x=${diffAvgNC} 'NR > 1 {delta = $5 - x; sum += (delta >= 0 ? delta : -delta)} END {print sum / (2 * (NR - 1))}'`

diffDelta=`cat ${outPathBase}/diffsum.txt | awk 'NR > 1 {sum += ($5 - $3)} END {print sum / (NR - 1)}'`

diffDeltaSD=`cat ${outPathBase}/diffsum.txt | awk -v x=${diffDelta} 'NR > 1 {delta = ($5 - $3) - x; sum += (delta >= 0 ? delta : -delta)} END {print sum / (2 * (NR - 1))}'`

bicAvg=`cat ${outPathBase}/emtabStats.txt | awk 'NR > 1 {sum += $5} END {print sum / (NR - 1)}'`

bicSD=`cat ${outPathBase}/emtabStats.txt | awk -v x=${bicAvg} 'NR > 1 {delta = ($5 - x); sum += (delta >= 0 ? delta : -delta)} END {print sum / (2 * (NR - 1))}'`

bicAvgNC=`cat ${outPathBase}/emtabStats.txt | awk 'NR > 1 {sum += $4} END {print sum / (NR - 1)}'`

bicSDNC=`cat ${outPathBase}/emtabStats.txt | awk -v x=${bicAvgNC} 'NR > 1 {delta = ($4 - x); sum += (delta >= 0 ? delta : -delta)} END {print sum / (2 * (NR - 1))}'`

bicDelta=`cat ${outPathBase}/emtabStats.txt | awk 'NR > 1 {sum += ($4 - $5)} END {print sum / (NR - 1)}'`

bicDeltaSD=`cat ${outPathBase}/emtabStats.txt | awk -v x=${bicDelta} 'NR > 1 {delta = ($4 - $5) -x; sum += (delta >= 0 ? delta : -delta)} END {print sum / (2 * (NR - 1))}'`

echo "#name mu length lengthSD  diff diffSD  diffNC diffNCSD diffDelta diffDeltaSD  BIC BICSD BICNC BICNCSD BICDelta BICDeltaSD" > ${outPathBase}/row.txt

echo "$name  $mu $lenAvg $lenSD $diffAvg $diffSD $diffAvgNC $diffSDNC $diffDelta $diffDeltaSD  $bicAvg $bicSD $bicAvgNC $bicSDNC $bicDelta $bicDeltaSD" >> ${outPathBase}/row.txt