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
cd ../../../../

# histogram.so
cd ./sova/computation/histogram
python3 setup.py build_ext -i
rm -rf build

# calculate_domain_centers.so                                                                       
cd ../split_and_merge/domain_centers
python3 setup.py build_ext -i
rm -rf  build

# calculate_domain_centers.so                                                                        
cd ../util/numpy_extension
python3 setup.py build_ext -i
rm -rf  build

# cif2cell                                                                                            
cd ../../../../libs/cif2cell
python3 setup.py build_ext -i
rm -rf  build
