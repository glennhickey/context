#!/bin/bash

inFile=${1}
outNameBase=${2}
numRuns=${3}
numSeq=${4}
maxLen=${5}
numEMIts=${6}
win=${7}
ecr=3
ecrth=1e-5
optits=500

# toggle for cluster / desktop
subCmd=""
#subCmd="qsub"

exeFile="../../src/axtscan"
outDir="../results/axt/${outNameBase}"
mkdir "${outDir}" 2> /dev/null

doneFile="${outNameBase}_axt_finish.sh"
echo '#!/bin/bash' > ${doneFile}
chmod u+x ${doneFile}

for (( i=1; i<=${numRuns}; i++ ))
do
    workDir="${outDir}"
    mkdir ${workDir} 2> /dev/null
    parFile="${workDir}/params.txt"
    echo "<*t= 1 mu= 0 ga= 0 RD= 0 RI= 0 RMD= 0 RMI= 0 PE= 0 PCD= 0 PCI= 0 KD= 0 KI= 0 >" > ${parFile}
    
    rnum=$RANDOM
    echo "${i}) ${workDir} ${inFile}"

    sfile=${outNameBase}_axt_${i}.sh
    echo '#!/bin/bash' > ${sfile}
    chmod u+x ${sfile}
    echo "pushd . > /dev/null" >> ${sfile}
    echo "cd ~/Documents/devel/context/test/scripts" >> ${sfile}

    echo ${exeFile} infile=${inFile} fpfile=${parFile} n=${numSeq} maxlen=${maxLen} win=${win} optits=${optits} emthreshold=-100 emits=${numEMIts} ecr=${ecr} ecrthreshold=${ecrth} seed=${rnum} outfile="${workDir}/emtab${i}.txt" distfile="${workDir}/gldist${i}.txt" sym=true ">" "${workDir}/emout${i}.txt" >> ${sfile}
                              
    echo ${exeFile} infile=${inFile} fpfile=${parFile} n=${numSeq} maxlen=${maxLen} win=2 optits=${optits} emthreshold=-100 emits=${numEMIts} ecr=${ecr} ecrthreshold=${ecrth} seed=${rnum} context=false outfile="${workDir}/emtab${i}_nc.txt" distfile="${workDir}/gldist${i}_nc.txt" sym=true ">" "${workDir}/emout${i}_nc.txt" >> ${sfile}

    echo "popd > /dev/null" >> ${sfile}
    echo "rm ${sfile}" >> ${doneFile}
    ${subCmd} ./${sfile}
done
                        
echo "R --vanilla --slave --args ${workDir}/emtab ${numRuns} < lrtemtab.R > lrtemtab.Rout" >> ${doneFile}

echo "R --vanilla --slave --args ${workDir}/emtab ${numRuns} < avgemtab.R > avgemtab.Rout" >> ${doneFile}

echo "./difftab.sh ${workDir} ${numRuns}" >> ${doneFile}

echo "./avgpar.sh ${workDir}" >> ${doneFile}

echo "./avgsum.sh ${workDir}" >> ${doneFile}

echo "./avgdist.sh ${workDir}" >> ${doneFile}
 
echo "rm ${doneFile}" >> ${doneFile}