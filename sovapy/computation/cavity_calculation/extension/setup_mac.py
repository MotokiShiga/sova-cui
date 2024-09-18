from setuptools import setup, Extension

setup(name = 'libalgorithm', version = '1.0.0',  \
      ext_modules = [
	      Extension('libalgorithm', 
						      sources=['algorithm.c'],
						      include_dirs=['./gr_mac/include'],
                              library_dirs=['./gr_mac/lib'],
                              libraries=['GR','GR3'],
						    #   runtime_library_dirs=['./gr/lib'],
						      extra_compile_args=["-O3", "-Wall", "-Wextra", "-fPIC"])
	    ])