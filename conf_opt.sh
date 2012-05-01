#!/bin/bash

./configure CXXFLAGS="-O3 -DNDEBUG -Wall" LDFLAGS="-L/opt/local/lib/" CPPFLAGS="-I/opt/local/include" --prefix=/Users/hickey/Documents/devel/context
mkdir ./temp232323
cp src/dblf84_*.o ./temp232323/
make clean
cp ./temp232323/*.o src/
rm -rf ./temp232323
make
