rem histogram.pyd
cd sova\computation\histogram
python setup.py build_ext -i
rmdir /s /q build

rem calculate_domain_centers.pyd
cd ..\split_and_merge\domain_centers
python setup.py build_ext -i
rmdir /s /q  build

rem calculate_domain_centers.pyd
cd ..\util\numpy_extension
python setup.py build_ext -i
rmdir /s /q  build

rem cif2cell
cd ..\..\..\..\libs\cif2cell
python setup.py build_ext -i
rmdir /s /q  build
