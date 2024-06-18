# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

import os, sys
sys.path.append(os.path.join(os.path.dirname('__file__'), '..'))

from core.file import File
from computation.cavity import Cavity

path = "../data/amorphous_rmc/sio.cfg"
f = File.open(path)
elements = ['Si','O']
atoms = f.getatoms(0, elements)

cavity = Cavity(atoms)
cavity.calculate()

domains = cavity.domains

print("Calculate domains")
print('Found {:d} domains'.format(len(domains.volumes)))
index =  0
print('Domain volume of index {} : {}'.format(index, domains.volumes[index]))
print('Found {:d} critical domains in file {}. Domain indices: {}'.format(
    len(domains.critical_domains), os.path.basename(path), 
    domains.critical_domains))
