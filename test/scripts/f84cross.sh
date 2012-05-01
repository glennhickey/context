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
fromscratch=${10}
optits=1000
optruns=2
rorder=true

exeSam="../../src/axtsample"
exeTrain="../../src/train"
exeAlign="../../src/align"

# toggle for cluster / desktop
subCmd=""
#subCmd="qsub"

outDir="../results/hky/${outNameBase}"
mkdir "${outDir}" 2> /dev/null

args="minlen=${minLen} maxlen=${maxLen} n=${numSeq} win=${win} emits=${numEMIts} offset=${offset} optruns=${optruns} optits=${optits} rorder=${rorder}"

jobname="${outNameBase}_cross_${numSeq}_${minLen}_${maxLen}_${win}"

head="f84test.sh.head"
echo '#!/bin/bash' > ${head}
chmod u+x ${head}
echo "pushd . > /dev/null" >> ${head}
echo "cd ~/Documents/devel/context/test/scripts" >> ${head}
tail="popd > /dev/null"

if [ $resample -ne "0" ]
    then
    
    ${exeSam} infile=../../data/chr22.hg19.mm9.net.axt outfile=${outDir}/mouse_${numSeq}_${minLen}_${maxLen}_${win}_cross.txt ${args}  

    ${exeSam} infile=../../data/chr22.hg19.canFam2.net.axt outfile=${outDir}/dog_${numSeq}_${minLen}_${maxLen}_${win}_cross.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.tarSyr1.net.axt outfile=${outDir}/tarsier_${numSeq}_${minLen}_${maxLen}_${win}_cross.txt ${args}

    ${exeSam} infile=../../data/chr22.hg19.calJac3.net.axt outfile=${outDir}/marmoset_${numSeq}_${minLen}_${maxLen}_${win}_cross.txt ${args}  

    ${exeSam} infile=../../data/chr22.hg19.rheMac2.net.axt outfile=${outDir}/macaque_${numSeq}_${minLen}_${maxLen}_${win}_cross.txt ${args}  

    ${exeSam} infile=../../data/chr22.hg19.gorGor1.net.axt outfile=${outDir}/gorilla_${numSeq}_${minLen}_${maxLen}_${win}_cross.txt ${args} 

    ${exeSam} infile=../../data/chr22.hg19.ponAbe2.net.axt outfile=${outDir}/orang_${numSeq}_${minLen}_${maxLen}_${win}_cross.txt ${args} 

    ${exeSam} infile=../../data/chr22.hg19.panTro2.net.axt outfile=${outDir}/chimp_${numSeq}_${minLen}_${maxLen}_${win}_cross.txt ${args} 

fi


sufname="align"
aname="${jobname}_${sufname}.sh"
cp ${head} ${aname}

for i in mouse dog tarsier marmoset macaque gorilla orang chimp 
do

    
#    for sufname in "sf84_djc sf84_dk sf84_df84"
#    for sufname in "sf84_djc" "sf84_dk"
    for sufname in "sf84_djc" "sf84_djc_nc" "sf84_dk" "sk_djc" "sk_djc_nc" "sk_dk" "sjc_djc" "sjc_djc_nc" "sjc_dk"
    do
        sname="${jobname}_${i}_${sufname}_al.sh"
        cp ${head} ${sname}
        echo "$exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_cross.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}_cross.txt > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}_cross_out.txt" >> ${sname}
         if [ $fromscratch -ne "0" ]
         then
             echo "$exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_cross.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_fs.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}_cross_fs.txt > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}_fs_cross_out.txt" >> ${sname}
         fi

         echo ${tail} >> ${sname}
         echo "${subCmd} ./${sname}" >> ${aname}

    done
done
