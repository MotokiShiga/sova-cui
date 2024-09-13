#!/bin/bash

mkdir external_lib
cd external_lib

#Clone the GR repo
git clone https://github.com/sciapp/gr.git -b v0.73.7 --depth 1 

#Build the GR
cd gr
mkdir build
cd build
cmake ..
cmake --build .

#Copy a source file
cd ..
cp lib/gr3/gr3_mc.c ../../sova/computation/cavity_calculation/extension/

#Build libalgorithm.so
cd ../../sova/computation/cavity_calculation/extension/
make
