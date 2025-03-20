from setuptools import setup, Extension
#from numpy.distutils.misc_util import get_numpy_include_dirs #deprecated from Python 3.12
import numpy

setup(
      name = 'find_index_of_first_element_not_equivalent', 
      version = '1.0.0',  \
      ext_modules = [
          Extension('find_index_of_first_element_not_equivalent', 
                    ['find_index_of_first_element_not_equivalent_macos.c'],
                    # include_dirs=[] + get_numpy_include_dirs(),
                    include_dirs=[] + [numpy.get_include()],
                    library_dirs=[],   
                    libraries=[],  
                    extra_compile_args=[], 
                    extra_link_args=[]
                    ) 
          ] 
      )
    
