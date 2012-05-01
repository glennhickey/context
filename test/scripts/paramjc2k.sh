#!/bin/bash

infile=$1
outside=$2
inside=$3

epsilon=0.0001

tempfile=${1}.temp
sed -e "s/\*//g" $1 > $tempfile

ina=`cat $1 | awk '{print $6}'`
inb=`cat $1 | awk '{print $8}'`
inc=`cat $1 | awk '{print $20}'`
ind=`cat $1 | awk '{print $22}'`

inrd=`cat $1 | awk '{print $32}'`
inri=`cat $1 | awk '{print $34}'`
inrmd=`cat $1 | awk '{print $36}'`
inrmi=`cat $1 | awk '{print $38}'`
inpcd=`cat $1 | awk '{print $42}'`
inpci=`cat $1 | awk '{print $44}'`
inkd=`cat $1 | awk '{print $46}'`
inki=`cat $1 | awk '{print $48}'`

a=$epsilon
b=`echo $inb $epsilon | awk '{print $1 - $2}'`

c=$epsilon
d=`echo $ind $epsilon | awk '{print $1 - $2}'`

if [ $outside -ne "0" ]
then
    sed -i -e "s/a\= ${ina}/a= $a/" $tempfile
    sed -i -e "s/b\= ${inb}/b= $b/" $tempfile
fi

if [ $inside -ne "0" ]
then
    sed -i -e  "s/c\= ${inc}/c= $c/" $tempfile
    sed -i -e  "s/d\= ${ind}/d= $d/" $tempfile
fi

cat $tempfile

rm $tempfile
