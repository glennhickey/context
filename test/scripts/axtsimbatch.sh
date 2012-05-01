numRuns=${1}
numSeq=${2}
maxLenIn=${3}
minLen=${4}
maxLen=${5}
numEMIts=${6}
fprob=${7}
win=${8}

./axtsim.sh mouse_${numRuns}_${numSeq}_${maxLenIn}_${numEMIts}_${win} mouse_${numRuns}_${numSeq}_${minLen}_${maxLen}_${numEMIts}_${fprob}_${win} ${numRuns} ${numSeq} ${minLen} ${maxLen} ${numEMIts} ${fprob} ${win}

./axtsim.sh dog_${numRuns}_${numSeq}_${maxLenIn}_${numEMIts}_${win} dog_${numRuns}_${numSeq}_${minLen}_${maxLen}_${numEMIts}_${fprob}_${win} ${numRuns} ${numSeq} ${minLen} ${maxLen} ${numEMIts} ${fprob} ${win}

./axtsim.sh tarsier_${numRuns}_${numSeq}_${maxLenIn}_${numEMIts}_${win} tarsier_${numRuns}_${numSeq}_${minLen}_${maxLen}_${numEMIts}_${fprob}_${win} ${numRuns} ${numSeq} ${minLen} ${maxLen} ${numEMIts} ${fprob} ${win}

./axtsim.sh marmoset_${numRuns}_${numSeq}_${maxLenIn}_${numEMIts}_${win} marmoset_${numRuns}_${numSeq}_${minLen}_${maxLen}_${numEMIts}_${fprob}_${win} ${numRuns} ${numSeq} ${minLen} ${maxLen} ${numEMIts} ${fprob} ${win}

./axtsim.sh macaque_${numRuns}_${numSeq}_${maxLenIn}_${numEMIts}_${win} macaque_${numRuns}_${numSeq}_${minLen}_${maxLen}_${numEMIts}_${fprob}_${win} ${numRuns} ${numSeq} ${minLen} ${maxLen} ${numEMIts} ${fprob} ${win}

./axtsim.sh gorilla_${numRuns}_${numSeq}_${maxLenIn}_${numEMIts}_${win} gorilla_${numRuns}_${numSeq}_${minLen}_${maxLen}_${numEMIts}_${fprob}_${win} ${numRuns} ${numSeq} ${minLen} ${maxLen} ${numEMIts} ${fprob} ${win}

./axtsim.sh chimp_${numRuns}_${numSeq}_${maxLenIn}_${numEMIts}_${win} chimp_${numRuns}_${numSeq}_${minLen}_${maxLen}_${numEMIts}_${fprob}_${win} ${numRuns} ${numSeq} ${minLen} ${maxLen} ${numEMIts} ${fprob} ${win}




