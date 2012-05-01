#!/bin/bash

outNameBase=${1}
numSeq=${2}
minLen=${3}
maxLen=${4}
numEMIts=${5}
offset=${6}
win=${7}
fromscratch=${8}
optits=250

exeSam="../../src/axtsample"
exeTrain="../../src/train"
exeAlign="../../src/align"

# toggle for cluster / desktop
subCmd=""
#subCmd="qsub"

outDir="../results/hky/${outNameBase}"
mkdir "${outDir}" 2> /dev/null

printf "species" "sf84_djc" "sf84_djc_nc" "sf84_dk" "sk_djc" "sk_djc_nc" "sk_dk" "sjc_djc" "sjc_djc_nc" "sjc_dk"

printf "\n"
for i in mouse dog tarsier marmoset macaque gorilla orang chimp 
do
    printf "$i  "
    for sufname in "sf84_djc" "sf84_djc_nc" "sf84_dk" "sk_djc" "sk_djc_nc" "sk_dk" "sjc_djc" "sjc_djc_nc" "sjc_dk"
    do
        fname="${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}_cross_out.txt"
        val=`tail -1 $fname | awk '{print $NF}'`
        printf "%s" "$val  "

	if [ $fromscratch -ne "0" ]
	    then
	    fnamefs="${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}_fs_cross_out.txt"
	    val=`tail -1 $fnamefs | awk '{print $NF}'`
	    printf "%s" "$val  "
	fi

    done
        
    printf "\n"
done
