from setuptools import setup, Extension

setup(name = 'histogram', version = '1.0.0',  \
      ext_modules = [Extension('histogram', ['histogram_py3.c'])])
    
