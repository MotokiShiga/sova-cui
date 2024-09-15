from setuptools import setup, Extension

setup(name = 'calculate_domain_centers', version = '1.0.0',  \
      ext_modules = [Extension('calculate_domain_centers', ['calculate_domain_centers_macos.c'])])
    
