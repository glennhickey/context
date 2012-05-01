#!/bin/bash

outNameBase=${1}
numSeq=${2}
minLen=${3}
maxLen=${4}
numEMIts=${5}
offset=${6}
win=${7}
resample=${8}
retrain=${9}
optits=250

exeSam="../../src/axtsample"
exeTrain="../../src/train"
exeAlign="../../src/align"

outDir="../results/hky/${outNameBase}"
mkdir "${outDir}" 2> /dev/null

args="minlen=${minLen} maxlen=${maxLen} n=${numSeq} win=${win} emits=${numEMIts} offset=${offset} optruns=1 optits=250"

if [ $resample -ne "0" ]
    then
    
    ${exeSam} infile=../../data/chr22.hg19.mm9.net.axt outfile=${outDir}/mouse_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.canFam2.net.axt outfile=${outDir}/dog_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.tarSyr1.net.axt outfile=${outDir}/tarsier_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.calJac3.net.axt outfile=${outDir}/marmoset_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.rheMac2.net.axt outfile=${outDir}/macaque_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.gorGor1.net.axt outfile=${outDir}/gorilla_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.panTro2.net.axt outfile=${outDir}/chimp_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

fi

if [ $retrain -ne "0" ]
then
#    for i in mouse dog tarsier marmoset macaque gorilla chimp
    for i in mouse dog macaque gorilla chimp
    do
#        $exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_fhky.txt &

#        $exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_fjc.txt sf84=false df84=false

#        $exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_fk.txt pflat=true fixa=0 qflat=true fixc=0 &

#        $exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_shkydk.txt qflat=true

        $exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_shkydjc.txt qflat=true fixc=0 &

        $exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_shkydjc0.txt qflat=true fixc=0  fixd=0.000000001

        $exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_shkydjc1.txt qflat=true fixc=0 fixd=0.1 &

        $exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_shkydjc2.txt qflat=true fixc=0 fixd=0.2

        $exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_shkydjc3.txt qflat=true fixc=0 fixd=0.3

    done
fi

for i in mouse dog tarsier marmoset macaque gorilla chimp
do
    printf "\n\n------- $i ---------\n"
    printf "Single: HKY  Double: HKY\n"
    $exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_fhky.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_fhky.txt

    printf "Single: JC  Double: JC\n"
    $exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_fjc.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_fjc.txt sf84=false df84=false

    printf "Single: KIM  Double: KIM\n"
    $exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_fk.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_fk.txt

    printf "Single: HKY  Double: KIM\n"
    $exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_shkydk.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_shkydk.txt

    printf "Single: HKY  Double: JC\n"
    $exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_shkydjc.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_shkydjc.txt

    printf "Single: HKY  Double: JC0\n"
    $exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_shkydjc0.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_shkydjc0.txt

    printf "Single: HKY  Double: JC1\n"
    $exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_shkydjc1.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_shkydjc1.txt

    printf "Single: HKY  Double: JC2\n"
    $exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_shkydjc2.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_shkydjc2.txt

    printf "Single: HKY  Double: JC3\n"
    $exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_shkydjc3.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_shkydjc3.txt

done
