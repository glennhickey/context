#!/bin/bash

srcNameBase=${1}
outNameBase=${2}
numRuns=${3}
numSeq=${4}
minLen=${5}
maxLen=${6}
numEMIts=${7}
fprob=${8}
win=${9}
ecr=3
ecrth=1e-5
optits=500

# toggle for cluster / desktop
subCmd=""
#subCmd="qsub"

exeFile="../../src/traintest"
outDir="../results/sim/${outNameBase}"
srcDir="../results/axt/${srcNameBase}"
mkdir "${outDir}" 2> /dev/null

doneFile="${outNameBase}_sim_finish.sh"
echo '#!/bin/bash' > ${doneFile}
chmod u+x ${doneFile}

for (( i=1; i<=${numRuns}; i++ ))
do
    workDir="${outDir}"
    mkdir ${workDir} 2> /dev/null
    parFile="${workDir}/params.txt"
    cp ${srcDir}/paramsEst.txt ${parFile}
    
    rnum=$RANDOM
    echo "${i}) ${workDir} ${inFile}"

    sfile=${outNameBase}_sim_${i}.sh
    echo '#!/bin/bash' > ${sfile}
    chmod u+x ${sfile}
    echo "pushd . > /dev/null" >> ${sfile}
    echo "cd ~/Documents/devel/context/test/scripts" >> ${sfile}

    echo ${exeFile} fpfile=${parFile} n=${numSeq} minlen=${minLen} maxlen=${maxLen} win=${win} optits=${optits} emthreshold=-100 emits=${numEMIts} ecr=${ecr} ecrthreshold=${ecrth} seed=${rnum} fprob=${fprob} sym=true outfile="${workDir}/emtab${i}.txt" ">" "${workDir}/emout${i}.txt" >> ${sfile}

    echo ${exeFile} fpfile=${parFile} n=${numSeq} minlen=${minLen} maxlen=${maxLen} win=2 optits=${optits} emthreshold=-100 emits=${numEMIts} ecr=${ecr} ecrthreshold=${ecrth} seed=${rnum} context=false fprob=${fprob} sym=true outfile="${workDir}/emtab${i}_nc.txt" ">" "${workDir}/emout${i}_nc.txt" >> ${sfile}

    echo "popd > /dev/null" >> ${sfile}
    echo "rm ${sfile}" >> ${doneFile}
    ${subCmd} ./${sfile}
done

echo "R --vanilla --slave --args ${workDir}/emtab ${numRuns} < lrtemtab.R > lrtemtab.Rout" >> ${doneFile}

echo "R --vanilla --slave --args ${workDir}/emtab ${numRuns} < avgemtab.R > avgemtab.Rout" >> ${doneFile}

echo "./difftab.sh ${workDir} ${numRuns}" >> ${doneFile}

echo "./avgpar.sh ${workDir}" >> ${doneFile}

echo "./avgsum.sh ${workDir}" >> ${doneFile}

echo "rm ${doneFile}" >> ${doneFile}
