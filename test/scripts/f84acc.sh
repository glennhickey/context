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

exeSim="../../src/simulate"
exeTrain="../../src/train"
exeAlign="../../src/align"
exeComp="../../src/compare"

# toggle for cluster / desktop
subCmd=""
#subCmd="qsub"

outDir="../results/hky/${outNameBase}"
mkdir "${outDir}" 2> /dev/null

args="minlen=${minLen} maxlen=${maxLen} n=${numSeq} win=${win} emits=${numEMIts} offset=${offset} optruns=${optruns} optits=${optits} rorder=${rorder} sym=true"

jobname="${outNameBase}_${numSeq}_${minLen}_${maxLen}_${win}_acc"

head="f84acc.sh.head"
echo '#!/bin/bash' > ${head}
chmod u+x ${head}
echo "pushd . > /dev/null" >> ${head}
echo "cd ~/Documents/devel/context/test/scripts" >> ${head}
tail="popd > /dev/null"

ognsuf="sf84_dk"

if [ $resample -ne "0" ]
    then
    for i in mouse dog tarsier marmoset macaque gorilla orang chimp 
    do
        ${exeSim} outfile=${outDir}/${i}_sim_${ognsuf}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}  fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${ognsuf}.txt fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}.txt
    done
fi

if [ $retrain -ne "0" ]
then
    for i in mouse dog tarsier marmoset macaque gorilla orang chimp 
#for i in mouse
    do


        sufname="sjc_djc"
        sname="${jobname}_${i}_${sufname}.sh"
        snamefs="${jobname}_${i}_${sufname}_fs.sh"
       cp ${head} ${sname}
       cp ${head} ${snamefs}

        echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}.txt 0 1 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_sf84_djci.txt" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_sf84_djci.txt qflat=true pflat=true fixa=0 fixc=0" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}_fs.txt qflat=true pflat=true fixa=0 fixc=0" >> ${snamefs}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}_nc.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_sf84_djci.txt qflat=true pflat=true fixa=0 context=false" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}_nc_fs.txt qflat=true pflat=true fixa=0 context=false" >> ${snamefs}

       echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt 0 0 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}i.txt" >> ${sname}

        sufname="sjc_dk"
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_sjc_djci.txt qflat=true pflat=true fixa=0" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}_fs.txt qflat=true pflat=true fixa=0" >> ${snamefs}

       echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt 0 0 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}i.txt" >> ${sname}

        sufname="sk_djc"
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_sjc_djci.txt qflat=true pflat=true fixc=0" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}_fs.txt qflat=true pflat=true fixc=0" >> ${snamefs}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}_nc.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_sjc_djci.txt qflat=true pflat=true context=false" >> ${sname}

 echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}_nc_fs.txt qflat=true pflat=true context=false" >> ${snamefs}

       echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt 0 0 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}i.txt" >> ${sname}

        sufname="sk_dk"
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_sk_djci.txt qflat=true pflat=true" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}_fs.txt qflat=true pflat=true" >> ${snamefs}

       echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt 0 0 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}i.txt" >> ${sname}

        sufname="sf84_djc"
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_sk_djci.txt qflat=true fixc=0" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}_fs.txt qflat=true fixc=0" >> ${snamefs}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}_nc.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_sk_djci.txt qflat=true context=false" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}_nc_fs.txt qflat=true context=false" >> ${snamefs}

       echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt 0 0 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}i.txt" >> ${sname}

        sufname="sf84_dk"
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_sf84_djci.txt qflat=true" >> ${sname}

 echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}_fs.txt qflat=true" >> ${snamefs}

       echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt 0 0 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}i.txt" >> ${sname}

        sufname="sf84_df84"
#        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_sf84_dki.txt" >> ${sname}

        echo ${tail} >> ${sname}
        echo ${tail} >> ${snamefs}
s
        ${subCmd} ./${sname}

        if [ $fromscratch -ne "0" ]
        then
            ${subCmd} ./${snamefs}
        fi
    done
fi

sufname="align"
aname="${jobname}_${sufname}.sh"
cp ${head} ${aname}

for i in mouse dog tarsier marmoset macaque gorilla orang chimp 
do
    for sufname in "sf84_djc" "sf84_djc_nc" "sf84_dk" "sk_djc" "sk_djc_nc" "sk_dk" "sjc_djc" "sjc_djc_nc" "sjc_dk"
    do
        sname="${jobname}_${i}_${sufname}_al.sh"
        cp ${head} ${sname}
        echo "$exeAlign infile=${outDir}/${i}_sim_${ognsuf}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_sim_${ognsuf}_${sufname}.txt > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_sim_${ognsuf}_${sufname}_out.txt" >> ${sname}

        echo "$exeComp infile=${outDir}/${i}_sim_${ognsuf}_${numSeq}_${minLen}_${maxLen}_${win}.txt infile2=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_sim_${ognsuf}_${sufname}.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_sim_${ognsuf}_${sufname}_comp.txt" >> ${sname}

         if [ $fromscratch -ne "0" ]
         then

        echo "$exeAlign infile=${outDir}/${i}_sim_${ognsuf}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sim_${ognsuf}_${sufname}_fs.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_sim_${ognsuf}_${sufname}_fs.txt > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_sim_${ognsuf}_${sufname}_fs_out.txt" >> ${sname}

        echo "$exeComp infile=${outDir}/${i}_sim_${ognsuf}_${numSeq}_${minLen}_${maxLen}_${win}.txt infile2=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_sim_${ognsuf}_${sufname}_fs.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_sim_${ognsuf}_${sufname}_fs_comp.txt" >> ${sname}
         fi

         echo ${tail} >> ${sname}
         echo "${subCmd} ./${sname}" >> ${aname}

    done

done
