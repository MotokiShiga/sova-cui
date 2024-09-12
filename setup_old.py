# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 18:16:52 2024

@author: H. Morita and M. Shiga
"""

from setuptools import setup, find_packages
from setuptools.command.install import install
import shutil
import os
import subprocess
import platform

class sova_install(install):
  description = "install SOVA"

  def initialize_options(self): 
      #self.my_data_dir = '/opt/myapp-data'
      install.initialize_options(self)

  def _pre_install(self):
      #os.mkdir('data1')
      if platform.system() == 'Windows':
          path = os.path.join(os.path.dirname(__file__), 'make.bat')
          p = subprocess.Popen(path)
          p.wait()
      else:
          p = subprocess.Popen('bash make.sh', shell=True)
          p.wait()

  def run(self): 
      self._pre_install()
      install.run(self)

  def get_outputs(self): 
      return install.get_outputs(self) + [self.my_data_dir]
        
setup(
    name='sova',
    version='0.5.0',
    install_requires=[
        'ase',
        'numpy',
        'scipy',
        'matplotlib',
        'h5py',
        'networkx==3.1',
        'numba',
        'tqdm',
        'spglib',
        'PyCifRW',
    ],
    author='Tohoku Univ.',
    author_email='motoki.shiga.b4@tohoku.ac.jp',
    url='https://www.shiga-lab.org/sova',
    
    include_package_data=True,
    packages=find_packages(),   
    cmdclass={'install': sova_install},
)
