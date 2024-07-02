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

