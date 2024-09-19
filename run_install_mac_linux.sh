#!/bin/bash

# histogram.so
cd ./sovapy/computation/histogram
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

#Build libalgorithm.so
cd ../../../sovapy/computation/cavity_calculation/extension/
make
cd ../../../../

#Build libalgorithm.so using built files (macOS)
# cd ./sovapy/computation/cavity_calculation/extension/
# xattr -rc gr_mac/lib/*
# cp gr_mac/lib/*dylib ./
# python setup_mac.py build_ext -i
# rm -rf  build
# cd ../../../../