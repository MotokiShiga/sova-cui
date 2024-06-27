# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 13:26:37 2024

@author: H. Morita
"""

from sova.core.data import ResultsFile

path = "./data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(path)
atoms = f.getatoms()

f = ResultsFile(atoms, filepath=path)
path = "./a_SiO2_speed1e11K_rand.hdf5"
f.write(filepath=path)

f.read(filepath=path)
atoms = f.atoms
