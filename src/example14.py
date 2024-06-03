# -*- coding: utf-8 -*-
"""
Created on Fri May  3 09:13:22 2024

@author: H. Morita
"""

import os
from core.file import File
from core import data
from core.calculation.algorithm import CavityCalculation, DomainCalculation, FakeDomainCalculation
from core.calculation.discretization import DiscretizationCache, AtomDiscretization

path = "./data/amorphous_rmc/sio.cfg"

f = File.open(path)
elements = ['Si','O']
atoms = f.getatoms(0, elements)

resolution = 64
cachedir = './cache'
cachepath = os.path.join(cachedir, 'discretization_cache.hdf5')
discretization_cache = DiscretizationCache(cachepath)
with DiscretizationCache(cachepath) as discretization_cache:
    discretization = discretization_cache.get_discretization(atoms.volume, resolution)
atom_discretization = AtomDiscretization(atoms, discretization)
domain_calculation = DomainCalculation(discretization, atom_discretization)

print("Calculating domains")
domains = data.Domains(domain_calculation)
print('Found {:d} domains'.format(len(domain_calculation.domain_volumes)))
index =  0
print('Domain volume of index {} : {}'.format(index, domain_calculation.domain_volumes[index]))
print('Found {:d} critical domains in file {}. Domain indices: {}'.format(
    len(domain_calculation.critical_domains), os.path.basename(path), 
    domain_calculation.critical_domains))

print()
print("Calculating surface-based cavities")
gyration_tensor_parameters = False
cavity_calculation = CavityCalculation(domain_calculation, use_surface_points=True,
                                       gyration_tensor_parameters=gyration_tensor_parameters)
surface_cavities = data.Cavities(cavity_calculation)
print('Found {:d} surface-based cavities'.format(len(cavity_calculation.cavity_volumes)))

print()
print("Calculating center-based cavities")
results = data.Results(path, 0, resolution, atoms,
                       domains=domains, 
                       surface_cavities=surface_cavities)

domain_calculation = FakeDomainCalculation(discretization, atom_discretization, results)
cavity_calculation = CavityCalculation(domain_calculation, use_surface_points=False,
                                       gyration_tensor_parameters=gyration_tensor_parameters)
results.center_cavities = data.Cavities(cavity_calculation)
print('Found {:d} surface-based cavities'.format(len(cavity_calculation.cavity_volumes)))