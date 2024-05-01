#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 28 12:42:00 2024

@author: morita
"""

import math
import numpy as np

try:
    from computation.histogram import histogram as hist
    has_hist = True
except ImportError:
    import computation.calc_histogram as hist
    has_hist = False
    print('Warning calculate histogram is slow')

def histogram(centres,metric,dr,symbols=None,truncated=False):
    
    print(centres)