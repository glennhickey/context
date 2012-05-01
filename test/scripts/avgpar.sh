#!/bin/bash

outPathBase=${1}
flist="emout[0-9]*.txt"
#flist="emout[0-9]*_nc.txt"

count=0
t=0
mu=0
ga=0
rd=0
ri=0
rmd=0
rmi=0
pe=0
pcd=0
pci=0
kd=0
ki=0

data=`ls ${outPathBase}/${flist} | grep ${flist}`

t=`tail -1 ${data} | awk '$1 == "<*t=" {sum+=$2; r+=1} END {print sum / r}'`

mu=`tail -1 ${data} | awk '$1 == "<*t=" {sum += $4; r+=1} END {print sum / r}'`

ga=`tail -1 ${data} | awk '$1 == "<*t=" {sum += $6; r+=1} END {print sum / r}'`

rd=`tail -1 ${data} | awk '$1 == "<*t=" {sum += $8; r+=1} END {print sum / r}'`

ri=`tail -1 ${data} | awk '$1 == "<*t=" {sum += $10; r+=1} END {print sum / r}'`

rmd=`tail -1 ${data} | awk '$1 == "<*t=" {sum += $12; r+=1} END {print sum / r}'`

rmi=`tail -1 ${data} | awk '$1 == "<*t=" {sum += $14; r+=1} END {print sum / r}'`

pe=`tail -1 ${data} | awk '$1 == "<*t=" {sum += $16; r+=1} END {print sum / r}'`

pcd=`tail -1 ${data} | awk '$1 == "<*t=" {sum += $18; r+=1} END {print sum / r}'`

pci=`tail -1 ${data} | awk '$1 == "<*t=" {sum += $20; r+=1} END {print sum / r}'`

kd=`tail -1 ${data} | awk '$1 == "<*t=" {sum += $22; r+=1} END {print sum / r}'`

ki=`tail -1 ${data} | awk '$1 == "<*t=" {sum += $24; r+=1} END {print sum / r}'`

echo "<*t= $t mu= $mu ga= $ga RD= $rd RI= $ri RMD= $rmd RMI= $rmi PE= $pe PCD= $pcd PCI= $pci KD= $kd KI= $ki >" > ${outPathBase}/paramsEst.txt


count=0
t=0
mu=0
ga=0
rd=0
ri=0
rmd=0
rmi=0
pe=0
pcd=0
pci=0
kd=0
ki=0
mse=0

data=`ls ${outPathBase}/${flist} | grep ${flist}`

t=`grep -h "Par_Delta" ${data} | sed -e "s/Par_Delta=//" | awk '$1 == "<*t=" {sum+=$2; r+=1} END {print sum / r}'`

mu=`grep -h "Par_Delta" ${data} | sed -e "s/Par_Delta=//" | awk '$1 == "<*t=" {sum += $4; r+=1} END {print sum / r}'`

ga=`grep -h "Par_Delta" ${data} | sed -e "s/Par_Delta=//" | awk '$1 == "<*t=" {sum += $6; r+=1} END {print sum / r}'`

rd=`grep -h "Par_Delta" ${data} | sed -e "s/Par_Delta=//" | awk '$1 == "<*t=" {sum += $8; r+=1} END {print sum / r}'`

ri=`grep -h "Par_Delta" ${data} | sed -e "s/Par_Delta=//" | awk '$1 == "<*t=" {sum += $10; r+=1} END {print sum / r}'`

rmd=`grep -h "Par_Delta" ${data} | sed -e "s/Par_Delta=//" | awk '$1 == "<*t=" {sum += $12; r+=1} END {print sum / r}'`

rmi=`grep -h "Par_Delta" ${data} | sed -e "s/Par_Delta=//" | awk '$1 == "<*t=" {sum += $14; r+=1} END {print sum / r}'`

pe=`grep -h "Par_Delta" ${data} | sed -e "s/Par_Delta=//" | awk '$1 == "<*t=" {sum += $16; r+=1} END {print sum / r}'`

pcd=`grep -h "Par_Delta" ${data} | sed -e "s/Par_Delta=//" | awk '$1 == "<*t=" {sum += $18; r+=1} END {print sum / r}'`

pci=`grep -h "Par_Delta" ${data} | sed -e "s/Par_Delta=//" | awk '$1 == "<*t=" {sum += $20; r+=1} END {print sum / r}'`

kd=`grep -h "Par_Delta" ${data} | sed -e "s/Par_Delta=//" | awk '$1 == "<*t=" {sum += $22; r+=1} END {print sum / r}'`

ki=`grep -h "Par_Delta" ${data} | sed -e "s/Par_Delta=//" | awk '$1 == "<*t=" {sum += $24; r+=1} END {print sum / r}'`

mse=`grep -h "Par_MSE" ${data} | sed -e "s/Par_MSE=//" | awk 'NR == 1 {sum += $1; r+=1} END {print sum / r}'`

echo "<*t= $t mu= $mu ga= $ga RD= $rd RI= $ri RMD= $rmd RMI= $rmi PE= $pe PCD= $pcd PCI= $pci KD= $kd KI= $ki >  MSE= $mse" > ${outPathBase}/paramsDelta.txt

