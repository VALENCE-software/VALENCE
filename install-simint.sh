#!/bin/bash
#Bash script to install simint

git clone https://github.com/simint-chem/simint-generator simint-generator
cd simint-generator
mkdir build; cd build
cmake ../
make
cd ..
python ./create.py -g build/generator/ostei -l 3 -p 3 outdir
mkdir outdir/build; cd outdir/build
cmake -DSIMINT_C_FLAGS=-g  -DENABLE_FORTRAN=ON -DSIMINT_VECTOR=scalar -DCMAKE_INSTALL_PREFIX=../../../simint ../
make;make install
