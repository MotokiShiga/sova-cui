#!/bin/bash

# histogram.so
cd ./sova/computation/histogram
python setup.py build_ext -i
rm -rf build

# calculate_domain_centers.so                                                                       
cd ../split_and_merge/domain_centers
python setup.py build_ext -i
rm -rf  build

# calculate_domain_centers.so                                                                        
cd ../util/numpy_extension
python setup.py build_ext -i
rm -rf  build

# cif2cell                                                                                            
cd ../../../../libs/cif2cell
python setup.py build_ext -i
rm -rf  build

cd ../../../
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

#Build libalgorithm.so using built files (macOS)
# cd ./sova/computation/cavity_calculation/extension/
# xattr -rc gr_mac/lib/*
# cp gr_mac/lib/*dylib ./
# python setup_mac.py build_ext -i
# rm -rf  build
# cd ../../../../