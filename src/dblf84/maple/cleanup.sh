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

            fn="dblhky_${i}${j}${k}.cpp"
            cp "map_${i}${j}${k}.txt" ${fn}
            sed -i -e "s/(double) //g" ${fn}
            sed -i -e "s/(int) //g" ${fn}
            sed -i -e "s/a\[/_a[/g" ${fn}
            sed -i -e "s/c\[/_c[/g" ${fn}
            sed -i -e "s/p\[/_p[/g" ${fn}
            sed -i -e "s/q\[/_q[/g" ${fn}
            sed -i -e "s/b/_b/g" ${fn}
            sed -i -e "s/d/_d/g" ${fn}
            sed -i -e "s/r/_r/g" ${fn}
            sed -i -e "s/mat${i}${j}${k} \= /void DblHKY::mat${i}${j}${k}() newline/" ${fn}
 
           sed -i -e 's/newline/\
{\
matrix/g' ${fn}
           
           sed -i -e "s/matrix/_mat[${ii}][${jj}][${kk}] = /g" ${fn}

           sed -i -e '1i\
#include <cmath>\
#include "dblhky.h"\
\
using namespace std;\
\

' ${fn}
            printf "\n}\n" >> ${fn}
            rm ${fn}-e
        done
    done
done
