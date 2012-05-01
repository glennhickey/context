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
optits=1000

exeSam="../../src/axtsample"
exeTrain="../../src/train"
exeAlign="../../src/align"

# toggle for cluster / desktop
subCmd=""
#subCmd="qsub"

outDir="../results/hky/${outNameBase}"
mkdir "${outDir}" 2> /dev/null

args="minlen=${minLen} maxlen=${maxLen} n=${numSeq} win=${win} emits=${numEMIts} offset=${offset} optruns=1 optits=250"

jobname="${outNameBase}_${numSeq}_${minLen}_${maxLen}_${win}"

head="f84test.sh.head"
echo '#!/bin/bash' > ${head}
chmod u+x ${head}
echo "pushd . > /dev/null" >> ${head}
echo "cd ~/Documents/devel/context/test/scripts" >> ${head}
tail="popd > /dev/null"

if [ $resample -ne "0" ]
    then

    ${exeSam} infile=../../data/chr22.hg19.galGal3.net.axt outfile=${outDir}/chicken_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}
    
    ${exeSam} infile=../../data/chr22.hg19.mm9.net.axt outfile=${outDir}/mouse_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.canFam2.net.axt outfile=${outDir}/dog_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.tarSyr1.net.axt outfile=${outDir}/tarsier_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.calJac3.net.axt outfile=${outDir}/marmoset_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.rheMac2.net.axt outfile=${outDir}/macaque_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.gorGor1.net.axt outfile=${outDir}/gorilla_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.ponAbe2.net.axt outfile=${outDir}/orang_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.panTro2.net.axt outfile=${outDir}/chimp_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr19.mm9.rn4.net.axt outfile=${outDir}/mouse_rat_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr19.mm9.oryCun2.net.axt outfile=${outDir}/mouse_rabbit_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr38.canFam2.felCat3.net.axt outfile=${outDir}/dog_cat_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

    ${exeSam} infile=../../data/chr4.dm3.droSec1.net.axt outfile=${outDir}/drom_dros_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}

fi

if [ $retrain -ne "0" ]
then
    for i in chicken mouse dog tarsier marmoset macaque gorilla orang chimp mouse_rat mouse_rabbit dog_cat drom_dros
    do
        
        sufname="sjc_djc"
        sname="${jobname}_${i}_${sufname}.sh"
       cp ${head} ${sname}
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt sf84=false df84=false" >> ${sname}
        echo ${tail} >> ${sname}
#        ${subCmd} ./${sname}

        sufname="sk_dk"
        sname="${jobname}_${i}_${sufname}.sh"
       cp ${head} ${sname}
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt pflat=true qflat=true" >> ${sname}
        echo ${tail} >> ${sname}
#        ${subCmd} ./${sname}

        sufname="sf84_dk"
        sname="${jobname}_${i}_${sufname}.sh"
       cp ${head} ${sname}
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt qflat=true" >> ${sname}
        echo ${tail} >> ${sname}
#        ${subCmd} ./${sname}

        sufname="sf84_djc"
        sname="${jobname}_${i}_${sufname}.sh"
       cp ${head} ${sname}
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt qflat=true fixc=0" >> ${sname}
        echo ${tail} >> ${sname}
        ${subCmd} ./${sname}

        for d in 0.000000001 0.05 0.1 0.15 0.2 0.25 0.3
        do
            sufname="sf84_djc_${d}"
            sname="${jobname}_${i}_${sufname}.sh"
           cp ${head} ${sname}
            echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt qflat=true fixc=0  fixd=${d}" >> ${sname}
            echo ${tail} >> ${sname}
            ${subCmd} ./${sname}
        done
    done
fi

sufname="align"
aname="${jobname}_${sufname}.sh"
cp ${head} ${aname}

for i in chicken mouse dog tarsier marmoset macaque gorilla orang chimp mouse_rat mouse_rabbit dog_cat drom_dros
do

    sufname="sjc_djc"
    sname="${jobname}_${i}_${sufname}_al.sh"
    cp ${head} ${sname}
    echo "$exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}.txt sf84=false df84=false > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}_out.txt" >> ${sname}
    echo ${tail} >> ${sname}
#    echo "${subCmd} ./${sname}" >> ${aname}
    
#    for sufname in "sk_dk" "sf84_dk" "sk_dk" "sf84_djc"
    for sufname in "sf84_djc"
    do
        sname="${jobname}_${i}_${sufname}_al.sh"
        cp ${head} ${sname}
        echo "$exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}.txt > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}_out.txt" >> ${sname}
        echo ${tail} >> ${sname}
        echo "${subCmd} ./${sname}" >> ${aname}
    done

    for d in 0.000000001 0.05 0.1 0.15 0.2 0.25 0.3
    do
        sufname="sf84_djc_${d}"
        sname="${jobname}_${i}_${sufname}_al.sh"
        cp ${head} ${sname}
        echo "$exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}.txt > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}_out.txt" >> ${sname}
        echo ${tail} >> ${sname}
        echo "${subCmd} ./${sname}" >> ${aname}
    done
done
