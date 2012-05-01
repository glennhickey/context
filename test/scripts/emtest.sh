#!/bin/bash

outNameBase=${1}
numRuns=${2}
numSeq=${3}
minLen=${4}
maxLen=${5}
numEMIts=${6}
fprob=${7}
win=${8}
ecr=3
ecrth=1e-5
optits=500

noRate=0.0
lowRate=0.01
medRate=0.1
highRate=0.5

pClose=0.5
pCont=0.75
pEnd=`echo "${minLen} ${maxLen}" | awk '{print 2 / ($1 + $2)}'`

# toggle for cluster / desktop
subCmd=""
#subCmd="qsub -p -100"

exeFile="../../src/traintest"
outDir="../results/em/${outNameBase}"
mkdir "${outDir}" 2> /dev/null

doneFile="${outNameBase}_em_finish.sh"
echo '#!/bin/bash' > ${doneFile}
chmod u+x ${doneFile}

currently=0
for iRate in ${noRate} ${lowRate} ${medRate} ${highRate}
#for iRate in 0.0000001
do
#    for dRate in ${noRate} ${lowRate} ${medRate} ${highRate}
    for dRate in ${iRate}
    do
        for ciRate in ${noRate} ${lowRate} ${medRate} ${highRate}
        do
#            for cdRate in ${noRate} ${lowRate} ${medRate} ${highRate}
            for cdRate in ${ciRate}
            do
                for sRate in ${lowRate} ${medRate} ${highRate}
                do
                    for xRate in ${lowRate} ${medRate} ${highRate}
                    do

                        workDir="${outDir}/out_i${iRate}_d${dRate}_ci${ciRate}_cd${cdRate}_s${sRate}_x${xRate}"
                        mkdir ${workDir} 2> /dev/null
                        parFile="${workDir}/params.txt"
                        echo "<*t= 1 mu= ${sRate} ga= ${xRate} RD= ${dRate} RI= ${iRate} RMD= ${cdRate} RMI= ${ciRate} PE= ${pEnd} PCD= ${pClose} PCI= ${pClose} KD= ${pCont} KI= ${pCont} >" > ${parFile}

                        for (( i=1; i<=${numRuns}; i++ ))
                        do
                            currently=$(($currently+1))
                            rnum=$RANDOM
                            echo "${currently}) ${workDir}"  
                          
                            sfile=${outNameBase}_em_i${iRate}_d${dRate}_ci${ciRate}_cd${cdRate}_s${sRate}_x${xRate}_${i}.sh
                            echo '#!/bin/bash' > ${sfile}
                            chmod u+x ${sfile}
                            echo "pushd . > /dev/null" >> ${sfile}
                            echo "cd ~/Documents/devel/context/test/scripts" >> ${sfile}

  
                            echo ${exeFile} fpfile=${parFile} n=${numSeq} minlen=${minLen} maxlen=${maxLen} win=${win} optits=${optits} emthreshold=-100 emits=${numEMIts} ecr=${ecr} ecrthreshold=${ecrth} seed=${rnum}  fprob=${fprob} sym=true outfile="${workDir}/emtab${i}.txt" ">" "${workDir}/emout${i}.txt" >> ${sfile}

                              
                            echo ${exeFile} fpfile=${parFile} n=${numSeq} minlen=${minLen} maxlen=${maxLen} win=2 optits=${optits} emthreshold=-100 emits=${numEMIts} ecr=${ecr} ecrthreshold=${ecrth} seed=${rnum} context=false  fprob=${fprob} sym=true outfile="${workDir}/emtab${i}_nc.txt" ">" "${workDir}/emout${i}_nc.txt" >> ${sfile}

                            echo "popd > /dev/null" >> ${sfile}
                            echo "rm ${sfile}" >> ${doneFile}
                            ${subCmd} ./${sfile}

                        done
                        
                        echo "R --vanilla --slave --args ${workDir}/emtab ${numRuns} < lrtemtab.R > lrtemtab.Rout" >> ${doneFile}

                        echo "R --vanilla --slave --args ${workDir}/emtab ${numRuns} < avgemtab.R > avgemtab.Rout" >> ${doneFile}

                        echo "./difftab.sh ${workDir} ${numRuns}" >> ${doneFile}

                        echo "./avgpar.sh ${workDir}" >> ${doneFile}

                        echo "./avgsum.sh ${workDir}" >> ${doneFile}

                    done
                done
            done
        done
    done
done

echo "./emsum.sh ${outNameBase}" >> ${doneFile}

echo "rm ${doneFile}" >> ${doneFile}