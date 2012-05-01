#!/bin/bash

targetName=${1}

tempFile="mergevtemp.txt"

echo "Name    " `head -1 ../results/axt/${2}/${targetName}` > ${tempFile}

shift 1
for i in $*
do
    seekDir="../results/axt/${i}"
    animal=${i%%_*}
    cat ${seekDir}/${targetName} | awk -v a=${animal} 'NR > 1{print a "  " $0}' >> ${tempFile}
done 

cat $tempFile

rm -f ${tempFile}


   