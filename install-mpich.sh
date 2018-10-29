#!/bin/bash

wget --no-check-certificate http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
tar -xf mpich-3.2.tar.gz
rm mpich-3.2.tar.gz
cd mpich-3.2
./configure --enable-fast=all  --enable-shared --prefix=/home/travis/
make -j4
make install
cd ..
