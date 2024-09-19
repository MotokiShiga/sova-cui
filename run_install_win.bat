rem Run using x64 Native Tools Command Prompt for Visual Studio 2022

rem Build histogram.so
cd ./sovapy/computation/histogram
python setup.py build_ext -i
rd /s build

rem Build calculate_domain_centers.so                                                                       
cd ../split_and_merge/domain_centers
python setup.py build_ext -i
rd /s build

rem Build calculate_domain_centers.so                                                                        
cd ../util/numpy_extension
python setup.py build_ext -i
rd /s build

rem Build cif2cell.so                                                                                      
cd ../../../../libs/cif2cell
python setup.py build_ext -i
rd /s build

cd ../../computation/cavity_calculation/extension

rem Run using  x64 Native Tools Command Prompt for Visual Studio 2022
cl /c /nologo /O2 /W3 /GL /DNDEBUG /MD -I./gr_win/include  /Tcalgorithm_win.c /Foalgorithm.obj

link /nologo /INCREMENTAL:NO /LTCG /DLL /MANIFEST:EMBED,ID=2 /MANIFESTUAC:NO /LIBPATH:./gr_win/lib libGR3static.lib algorithm.obj /OUT:libalgorithm.dll 

rem copy .\gr_win\dll\vcruntime140.dll .\
del /f *.obj *.lib *.exp

rem If you meet problems of dll dependencies, run "dumpbin /DEPENDENTS  xxx.dll".
rem Command dumpbin is useful to check dependent dll files.