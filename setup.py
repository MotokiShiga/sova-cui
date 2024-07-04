# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 18:16:52 2024

@author: H. Morita
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
          subprocess.Popen(path)
      else:
          os.system("./make.sh")

  def run(self): 
      self._pre_install()
      install.run(self)

  def get_outputs(self): 
      return install.get_outputs(self) + [self.my_data_dir]
        
setup(
    name='sova',
    version='0.1.0',
    install_requires=[
        'numpy==1.23.5',
        'scipy==1.8.1',
        'ase==3.22.1',
        'matplotlib',
        'h5py',
    ],
    author='Tohoku Univ.',
    author_email='motoki.shiga.b4@tohoku.ac.jp',
    url='https://www.shiga-lab.org/sova',
    
    include_package_data=True,
    packages=find_packages(),   
    cmdclass={'install': sova_install},
)
