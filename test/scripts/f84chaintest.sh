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

jobname="${outNameBase}_${numSeq}_${minLen}_${maxLen}_${win}"

head="f84test.sh.head"
echo '#!/bin/bash' > ${head}
chmod u+x ${head}
echo "pushd . > /dev/null" >> ${head}
echo "cd ~/Documents/devel/context/test/scripts" >> ${head}
tail="popd > /dev/null"

if [ $resample -ne "0" ]
    then
    
    ${exeSam} infile=../../data/chr22.hg19.mm9.net.axt outfile=${outDir}/mouse_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}  fpout=${outDir}/mouse_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_axt.txt

    ${exeSam} infile=../../data/chr22.hg19.canFam2.net.axt outfile=${outDir}/dog_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}  fpout=${outDir}/dog_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_axt.txt

    ${exeSam} infile=../../data/chr22.hg19.tarSyr1.net.axt outfile=${outDir}/tarsier_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}  fpout=${outDir}/tarsier_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_axt.txt

    ${exeSam} infile=../../data/chr22.hg19.calJac3.net.axt outfile=${outDir}/marmoset_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}  fpout=${outDir}/marmoset_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_axt.txt

    ${exeSam} infile=../../data/chr22.hg19.rheMac2.net.axt outfile=${outDir}/macaque_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args}  fpout=${outDir}/macaque_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_axt.txt

    ${exeSam} infile=../../data/chr22.hg19.gorGor1.net.axt outfile=${outDir}/gorilla_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpout=${outDir}/gorilla_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_axt.txt

    ${exeSam} infile=../../data/chr22.hg19.ponAbe2.net.axt outfile=${outDir}/orang_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpout=${outDir}/orang_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_axt.txt

    ${exeSam} infile=../../data/chr22.hg19.panTro2.net.axt outfile=${outDir}/chimp_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpout=${outDir}/chimp_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_axt.txt

echo     ${exeSam} infile=../../data/chr22.hg19.panTro2.net.axt outfile=${outDir}/chimp_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpout=${outDir}/chimp_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_axt.txt

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

        echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_axt.txt 0 1 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sf84_djci.txt" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sf84_djci.txt qflat=true pflat=true fixa=0 fixc=0" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_fs.txt qflat=true pflat=true fixa=0 fixc=0" >> ${snamefs}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_nc.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sf84_djci.txt qflat=true pflat=true fixa=0 context=false" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_nc_fs.txt qflat=true pflat=true fixa=0 context=false" >> ${snamefs}

       echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt 0 0 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}i.txt" >> ${sname}

        sufname="sjc_dk"
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sjc_djci.txt qflat=true pflat=true fixa=0" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_fs.txt qflat=true pflat=true fixa=0" >> ${snamefs}

       echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt 0 0 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}i.txt" >> ${sname}

        sufname="sk_djc"
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sjc_djci.txt qflat=true pflat=true fixc=0" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_fs.txt qflat=true pflat=true fixc=0" >> ${snamefs}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_nc.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sjc_djci.txt qflat=true pflat=true context=false" >> ${sname}

 echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_nc_fs.txt qflat=true pflat=true context=false" >> ${snamefs}

       echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt 0 0 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}i.txt" >> ${sname}

        sufname="sk_dk"
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sk_djci.txt qflat=true pflat=true" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_fs.txt qflat=true pflat=true" >> ${snamefs}

       echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt 0 0 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}i.txt" >> ${sname}

        sufname="sf84_djc"
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sk_djci.txt qflat=true fixc=0" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_fs.txt qflat=true fixc=0" >> ${snamefs}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_nc.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sk_djci.txt qflat=true context=false" >> ${sname}

        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_nc_fs.txt qflat=true context=false" >> ${snamefs}

       echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt 0 0 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}i.txt" >> ${sname}

        sufname="sf84_dk"
        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sf84_djci.txt qflat=true" >> ${sname}

 echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_fs.txt qflat=true" >> ${snamefs}

       echo "./paramjc2k.sh ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt 0 0 > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}i.txt" >> ${sname}

        sufname="sf84_df84"
#        echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_sf84_dki.txt" >> ${sname}

        echo ${tail} >> ${sname}
        echo ${tail} >> ${snamefs}

        ${subCmd} ./${sname}

        if [ $fromscratch -ne "0" ]
        then
            ${subCmd} ./${snamefs}
        fi

        for d in 0.000000001 0.05 0.1 0.15 0.2 0.25 0.3
        do
            sufname="sf84_djc_${d}"
#            sname="${jobname}_${i}_${sufname}.sh"
#           cp ${head} ${sname}
#            echo "$exeTrain infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt win=${win} ${args} fpout=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt qflat=true fixc=0  fixd=${d}" >> ${sname}
#            echo ${tail} >> ${sname}
#            ${subCmd} ./${sname}
        done
    done
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
        echo "$exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}.txt > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}_out.txt" >> ${sname}
         if [ $fromscratch -ne "0" ]
         then
             echo "$exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}_fs.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}_fs.txt > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}_fs_out.txt" >> ${sname}
         fi

         echo ${tail} >> ${sname}
         echo "${subCmd} ./${sname}" >> ${aname}

    done

#    for d in 0.000000001 0.05 0.1 0.15 0.2 0.25 0.3
#    do
#        sufname="sf84_djc_${d}"
#        sname="${jobname}_${i}_${sufname}_al.sh"
#        cp ${head} ${sname}
#        echo "$exeAlign infile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}.txt ${args} fpfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_par_${sufname}.txt outfile=${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}.txt > ${outDir}/${i}_${numSeq}_${minLen}_${maxLen}_${win}_${offset}_${numEMIts}_al_${sufname}_out.txt" >> ${sname}
#        echo ${tail} >> ${sname}
#        echo "${subCmd} ./${sname}" >> ${aname}
#    done
done
