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

class my_install(install):
  description = "install myapp"

  # 自分の処理に必要なオプションを追加
  #user_options = install.user_options + [
  #  ('my-data-dir=', 'd', "base directory for installing my data files." ),
  #]
  def initialize_options(self): 
      self.my_data_dir = '/opt/myapp-data'
      install.initialize_options(self)

  def _pre_install(self):
      #os.mkdir('data1')
      path = os.path.join(os.path.dirname(__file__), 'make.bat')
      subprocess.Popen(path)
      #shutil.copytree('./data', self.my_data_dir)

  def run(self): 
      self._pre_install()
      install.run(self)

  def get_outputs(self): 
      # get_outputsは--recordオプション指定時に呼ばれ、install時に作成したファイルとディレクトリのリストを返す。 
      # pip uninstall をした時に削除してほしいものがあれば、ここで返すようにする。 
      return install.get_outputs(self) + [self.my_data_dir]
        
setup(
    name='sova',
    version='0.1.0',
    
    author='Tohoku Univ.',
    author_email='motoki.shiga.b4@tohoku.ac.jp',
    url='https://www.shiga-lab.org/sova',
    
    include_package_data=True,
    packages=find_packages(),   
    cmdclass={'install': my_install},
)