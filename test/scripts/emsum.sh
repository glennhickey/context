#!/bin/bash

outNameBase=${1}

exeFile="../../src/traintest"
outDir="../results/em/${outNameBase}"

pushd . > /dev/null

cd ${outDir}

echo "#mu  indel  context  pvMin  pvMax  pvAvg  dbMin  dbMax  dbAvg" > emall.txt
echo "#mu  indel  context lMin  lMax  lAvg  dncMin  dncMax  dncAvg  dcMin  dcMax  dcAvg" > diffall.txt

for subname in out*
do
    echo $subname
    pfile=${subname}/params.txt
    statfile=${subname}/emtabStats.txt
    difffile=${subname}/diffsum.txt

    mu=`cat ${pfile} | awk '{print $4}'`
    indel=`cat ${pfile} | awk '{print $8}'`
    context=`cat ${pfile} | awk '{print $12}'`

     dbAvg=`sed '1d' ${statfile} | awk '{sum+=($5-$4)} END { print sum/NR}'` 
     dbMin=`sed '1d' ${statfile} | awk 'NR==1 {min=$5-$4} $5-$4 < min {min=($5-$4)} END {print min}'` 
     dbMax=`sed '1d' ${statfile} | awk 'NR==1 {max=$5-$4} $5-$4 > max {max=($5-$4)} END {print max}'`
     pvAvg=`sed '1d' ${statfile} | awk '{sum+=$6} END { print sum/NR}'` 
     pvMin=`sed '1d' ${statfile} | awk 'NR==1 {min=$6} $6 < min {min=$6} END { print min}'`
     pvMax=`sed '1d' ${statfile} | awk 'NR==1 {max=$6} $6 > max {max=$6} END { print max}'`
     
     echo ${mu}  ${indel}  ${context}  ${pvMin}  ${pvMax}  ${pvAvg}  ${dbMin}  ${dbMax}  ${dbAvg} >> emall.txt

     lAvg=`sed '1d' ${difffile} | awk '{sum+=$2} END { print sum/NR}'` 
     lMin=`sed '1d' ${difffile} | awk 'NR==1 {min=$2} $2 < min {min=$2} END { print min}'`
     lMax=`sed '1d' ${difffile} | awk 'NR==1 {max=$2} $2 > max {max=$2} END { print max}'`

     dncAvg=`sed '1d' ${difffile} | awk '{sum+=$5} END { print sum/NR}'` 
     dncMin=`sed '1d' ${difffile} | awk 'NR==1 {min=$5} $5 < min {min=$5} END { print min}'`
     dncMax=`sed '1d' ${difffile} | awk 'NR==1 {max=$5} $5 > max {max=$5} END { print max}'`

     dcAvg=`sed '1d' ${difffile} | awk '{sum+=$3} END { print sum/NR}'` 
     dcMin=`sed '1d' ${difffile} | awk 'NR==1 {min=$3} $3 < min {min=$3} END { print min}'`
     dcMax=`sed '1d' ${difffile} | awk 'NR==1 {max=$3} $3 > max {max=$3} END { print max}'`

     echo ${mu}  ${indel}  ${context}  ${lMin}  ${lMax}  ${lAvg}  ${dncMin}  ${dncMax}  ${dncAvg}  ${dcMin}  ${dcMax}  ${dcAvg} >> diffall.txt

done

popd > /dev/null