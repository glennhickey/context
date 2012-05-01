numRuns=${1}
numSeq=${2}
maxLen=${3}
numEMIts=${4}
win=${5}

./axttest.sh ../../data/chr22.hg19.mm9.net.axt mouse_${numRuns}_${numSeq}_${maxLen}_${numEMIts}_${win} ${numRuns} ${numSeq} ${maxLen} ${numEMIts} ${win}

./axttest.sh ../../data/chr22.hg19.canFam2.net.axt dog_${numRuns}_${numSeq}_${maxLen}_${numEMIts}_${win} ${numRuns} ${numSeq} ${maxLen} ${numEMIts} ${win}

#./axttest.sh ../../data/chr22.hg19.tarSyr1.net.axt tarsier_${numRuns}_${numSeq}_${maxLen}_${numEMIts}_${win} ${numRuns} ${numSeq} ${maxLen} ${numEMIts} ${win}

#./axttest.sh ../../data/chr22.hg19.calJac3.net.axt marmoset_${numRuns}_${numSeq}_${maxLen}_${numEMIts}_${win} ${numRuns} ${numSeq} ${maxLen} ${numEMIts} ${win}

./axttest.sh ../../data/chr22.hg19.rheMac2.net.axt macaque_${numRuns}_${numSeq}_${maxLen}_${numEMIts}_${win} ${numRuns} ${numSeq} ${maxLen} ${numEMIts} ${win}

./axttest.sh ../../data/chr22.hg19.gorGor1.net.axt gorilla_${numRuns}_${numSeq}_${maxLen}_${numEMIts}_${win} ${numRuns} ${numSeq} ${maxLen} ${numEMIts} ${win}

./axttest.sh ../../data/chr22.hg19.panTro2.net.axt chimp_${numRuns}_${numSeq}_${maxLen}_${numEMIts}_${win} ${numRuns} ${numSeq} ${maxLen} ${numEMIts} ${win}




