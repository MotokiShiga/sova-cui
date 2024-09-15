from setuptools import setup, Extension
setup(name='calculate_atomsdata',
        version='1.0',
        ext_modules=[Extension('calculate_atomsdata', ['calculate_atomsdata.c'])]
)
