from setuptools import setup, Extension
from numpy.distutils.misc_util import get_numpy_include_dirs

setup(
      name = 'find_index_of_first_element_not_equivalent', 
      version = '1.0.0',  \
      ext_modules = [
          Extension('find_index_of_first_element_not_equivalent', 
                    ['find_index_of_first_element_not_equivalent_macos.c'],
                    include_dirs=[] + get_numpy_include_dirs(),  
                    library_dirs=[],   
                    libraries=[],  
                    extra_compile_args=[], 
                    extra_link_args=[]
                    ) 
          ] 
      )
    
