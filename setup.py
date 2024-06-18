# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 18:16:52 2024

@author: H. Morita
"""

from setuptools import setup, find_packages

setup(
    name='sova',
    version='0.1.0',
    
    author='Tohoku Univ.',
    author_email='motoki.shiga.b4@tohoku.ac.jp',
    url='https://www.shiga-lab.org/sova',
    
    include_package_data=True,
    packages=find_packages(),      
)