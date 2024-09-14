#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 09:31:54 2023

@author: morita
"""

from distutils.core import setup, Extension

setup(name = 'calculate_domain_centers', version = '1.0.0',  \
      ext_modules = [Extension('calculate_domain_centers', ['calculate_domain_centers_macos.c'])])
    
