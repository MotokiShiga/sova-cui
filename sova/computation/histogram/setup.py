#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 09:31:54 2023

@author: morita
"""

from distutils.core import setup, Extension

setup(name = 'histogram', version = '1.0.0',  \
      ext_modules = [Extension('histogram', ['histogram_py3.c'])])
    
