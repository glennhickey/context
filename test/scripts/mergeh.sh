#!/bin/bash

targetName=${1}

tempFile=mergehtemp.txt
tempFile2=mergehtemp2.txt
tempFile3=mergehtemp3.txt
echo " " > ${tempFile}
echo " " > ${tempFile2}

shift 1
for i in $*
do
    seekDir="../results/axt/${i}"
    animal=${i%%_*}
    nf=`cat ${seekDir}/${targetName} | awk 'NR == 1 {print NF}'`

    rm -f ${tempFile2}
    for ((j = 0 ; j < $nf ; j++))
    do
        printf "${animal}  " >> ${tempFile2} 
    done
    printf "\n" >> ${tempFile2}

    cat ${seekDir}/${targetName} >> ${tempFile2}
    
    paste ${tempFile2} ${tempFile} > ${tempFile3}
    cp ${tempFile3} ${tempFile}
done 

cat $tempFile

rm -f $tempFile $tempFile2 $tempFile3


   