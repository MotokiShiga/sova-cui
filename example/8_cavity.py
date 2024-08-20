# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita and M. Shiga
"""
import os
from sova.core.file import File
from sova.computation.cavity import Cavity

# Load structural information from a xyz file
structure_file = "./data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(structure_file)

# Get atomic and cell (simulation box) data
atoms = f.getatoms()

print("Atom symbols:", atoms.symbols)

# Cavity analysis
cavity = Cavity(atoms)

# Calculate cavities
cavity.calculate()


# Print caluculation results
print("Calculate domains")
print('Found {:d} domains'.format(len(cavity.domains.volumes)))
index =  0
print('Domain volume of index {} : {}'.format(index, cavity.domains.volumes[index]))
print('Found {:d} critical domains in file {}. Domain indices: {}'.format(
    len(cavity.domains.critical_domains), os.path.basename(structure_file), 
    cavity.domains.critical_domains))
