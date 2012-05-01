#!/bin/bash

./configure CXXFLAGS="-g -Wall" LDFLAGS="-L/opt/local/lib/" CPPFLAGS="-I/opt/local/include"
mkdir ./temp232323
cp src/dblf84_*.o ./temp232323/
make clean
cp ./temp232323/*.o src/
rm -rf ./temp232323
make
