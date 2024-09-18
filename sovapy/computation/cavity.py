# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 15:08:59 2024

@author: H. Morita
"""

import os
from ..core import data
from .cavity_calculation.algorithm import CavityCalculation, DomainCalculation, FakeDomainCalculation
from .cavity_calculation.discretization import DiscretizationCache, AtomDiscretization


class Cavity(object):
    def __init__(self, atoms):
        self.atoms = atoms
        self.results = None
        
    @property
    def domains(self):
        return self.results.domains

    @property
    def surface_cavities(self):
        return self.results.surface_cavities

    @property
    def center_cavities(self):
        return self.results.center_cavities
        
    def calculate(self, resolution=64):
        self.resolution = resolution
        cachedir = './cache'
        cachepath = os.path.join(cachedir, 'discretization_cache.hdf5')
        discretization_cache = DiscretizationCache(cachepath)
        with DiscretizationCache(cachepath) as discretization_cache:
            discretization = discretization_cache.get_discretization(self.atoms.volume, resolution)
        atom_discretization = AtomDiscretization(self.atoms, discretization)
        domain_calculation = DomainCalculation(discretization, atom_discretization)

        #print("Calculating domains")
        domains = data.Domains(domain_calculation)
        #print('Found {:d} domains'.format(len(domain_calculation.domain_volumes)))
        #index =  0
        #print('Domain volume of index {} : {}'.format(index, domain_calculation.domain_volumes[index]))
        #print('Found {:d} critical domains in file {}. Domain indices: {}'.format(
        #    len(domain_calculation.critical_domains), os.path.basename(path), 
        #    domain_calculation.critical_domains))

        #print()
        #print("Calculating surface-based cavities")
        gyration_tensor_parameters = False
        cavity_calculation = CavityCalculation(domain_calculation, use_surface_points=True,
                                               gyration_tensor_parameters=gyration_tensor_parameters)
        surface_cavities = data.Cavities(cavity_calculation)
        #print('Found {:d} surface-based cavities'.format(len(cavity_calculation.cavity_volumes)))

        #print()
        #print("Calculating center-based cavities")
        path = ''
        self.results = data.Results(path, 0, self.resolution, self.atoms,
                               domains=domains, 
                               surface_cavities=surface_cavities)

        domain_calculation = FakeDomainCalculation(discretization, atom_discretization, self.results)
        cavity_calculation = CavityCalculation(domain_calculation, use_surface_points=False,
                                               gyration_tensor_parameters=gyration_tensor_parameters)
        self.results.center_cavities = data.Cavities(cavity_calculation)
        #print('Found {:d} surface-based cavities'.format(len(cavity_calculation.cavity_volumes)))