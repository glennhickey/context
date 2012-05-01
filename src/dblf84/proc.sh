#!/bin/bash

for i in 1 2 3 4
do
    for j in 1 2 3 4
    do
        for k in 1 2 3 4
        do
           ii=$((i - 1))
           jj=$((j - 1))
           kk=$((k - 1))
           
           fn="dblf84_${i}${j}${k}.cpp"
           
           sed -i -e "s/_mat\[${ii}\]\[${jj}\]\[${kk}\]/double val/" $fn
           sed -i -e "s/}//" $fn
           echo "_tab.setMD((DNA)${ii}, (DNA)${kk}, (DNA)${jj}, val);" >> $fn
           echo "_tab.setMI((DNA)${ii}, (DNA)${jj}, (DNA)${kk}, val);" >> $fn
           echo "}" >> $fn
        done
    done
done
