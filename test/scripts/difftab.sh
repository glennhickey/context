#!/bin/bash

outPathBase=${1}
numRuns=${2}

echo "#trial  length  diff  lengthNC  diffNC" > ${outPathBase}/diffsum.txt

for (( i=1; i<=${numRuns}; i++ ))
do
    # get the total length
    len=`tail -1 ${outPathBase}/emtab${i}.txt | awk '{print $2}'`
    # sanity check
    lenNC=`tail -1 ${outPathBase}/emtab${i}_nc.txt | awk '{print $2}'`

    # get the total difference
    tdiff=`grep -a tdiff ${outPathBase}/emout${i}.txt | awk '{print $NF; exit 0}'`
    tdiffNC=`grep -a tdiff ${outPathBase}/emout${i}_nc.txt | awk '{print $NF; exit 0}'`

    echo $i  ${len}  ${tdiff}  ${lenNC}  ${tdiffNC} >> ${outPathBase}/diffsum.txt
    
done

