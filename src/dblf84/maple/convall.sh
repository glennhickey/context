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
            cp "mout/map2f84_${i}${j}${k}.txt" ${fn}_pre

            sed -i -e "s/\[1\]/[ZERO]/g" ${fn}_pre
            sed -i -e "s/\[2\]/[ONE]/g" ${fn}_pre
            sed -i -e "s/\[3\]/[TWO]/g" ${fn}_pre
            sed -i -e "s/\[4\]/[THREE]/g" ${fn}_pre

            sed -i -e "s/ZERO/0/g" ${fn}_pre
            sed -i -e "s/ONE/1/g" ${fn}_pre
            sed -i -e "s/TWO/2/g" ${fn}_pre
            sed -i -e "s/THREE/3/g" ${fn}_pre

            sed -i -e "s/a\[/_a[/g" ${fn}_pre
            sed -i -e "s/c\[/_c[/g" ${fn}_pre
            sed -i -e "s/p\[/_p[/g" ${fn}_pre
            sed -i -e "s/q\[/_q[/g" ${fn}_pre
            sed -i -e "s/b/_b/g" ${fn}_pre
            sed -i -e "s/d/_d/g" ${fn}_pre
            sed -i -e "s/r/_r/g" ${fn}_pre

            ./map2c ${fn}_pre ${fn}
echo            ./map2c ${fn}_pre ${fn}

            
            sed -i -e "s/funname/mat${i}${j}${k}/" ${fn}
				sed -i -e "s/varname/_mat[${ii}][${jj}][${kk}]/" ${fn}
            
            rm -f ${fn}-e ${fn}_pre-e

   #         rm ${fn}_pre
        done
    done
done
