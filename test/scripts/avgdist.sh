#!/bin/bash

outPathBase=${1}
flist="gldist[0-9]*.txt"
#flist="gldist[0-9]*_nc.txt"

outFile=${outPathBase}/glDistTot.txt
tempFile=${outPathBase}/glDistTemp.txt

files=`ls ${outPathBase}/${flist} | grep ${flist}`

currently=0

for i in $files
do
    if [ $currently -eq 0 ]
    then
        cp $i $tempFile
    else
        head -1 $i > ${outFile}
        paste $i $tempFile | awk 'NR == 1 {next} {print $1 "  " $2 + $7 "  " $3 + $8 "  " $4 + $9 "  " $5 + $10}' >> ${outFile}
	cp ${outFile} ${tempFile}
    fi
    currently=$(($currently+1))
done

echo "len  num  tot  frloc  frglob  numNC  totNC  frlocNC  frglobNC" > ${outFile}

cat ${tempFile} | awk 'NR == 1 {next} {print $1 "  " $2 "  " $3 "  " $2 / $3 "  " $2 / ($3 + $5) "  " $4 "  " $5 "  " $4 / $5 "  "  $4 / ($3 + $5)}' >> ${outFile}

rm -f ${tempFile}

