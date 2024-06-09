# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 15:21:55 2024

@author: H. Morita
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from core.file import File

path = "../data/amorphous_rmc/sio.cfg"
f = File.open(path)
elements = ['Si','O']
atoms = f.getatoms(0, elements)

print(dir(atoms))

#atoms.bonds
#for i in range(1000,1100):
#    print(len(atoms.bonds[i]))
