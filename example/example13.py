# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 13:26:37 2024

@author: H. Morita
"""

from sova.core.file import File
from sova.core.data import Results

path = "../data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(path)
atoms = f.getatoms()

results = Results(atoms=atoms, filepath=path)
