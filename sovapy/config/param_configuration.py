# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 16:57:33 2022

@author: H. Morita
"""

from .configobj import ConfigObj
from distutils.util import strtobool

class ParameterConfiguration(object):
    """
    edit configuration class
    """
    class atoms(object):
        """
            edit atoms parameter class
        """
        def __init__(self):
            self.opacity = 1.0
            self.scale_factor = 0.0
            self.elems_color = {}              
        def clear(self):
            self.elems_color.clear()
            
    class bonds(object):
        """
            edit bonds parameter class
        """
        def __init__(self):
            self.elems_param = {}
            self.opacity = 1.0
            self.scale_factor = 0.0
            self.color_mode = False
            self.color = [0,0,0]            
        
        def add_bond(self,pair,param):
            self.elems_param[pair] = param
        
        def clear(self):
            self.elems_param.clear()
            
        def get_max_length(self, symbol1, symbol2):
            """
            get max length

            Parameters
            ----------
            symbol1 : string
                Atomic name.
            symbol2 : string
                Atomic name.

            Returns
            -------
            float
                max radius.

            """
                        
            try:
                r = float(self.elems_param[symbol1 + '-' + symbol2][1]) 
            except:
                pass
            try:
                r = float(self.elems_param[symbol2 + '-' + symbol1][1]) 
            except :
                print('Warning : wrong bond pair (max) ', symbol1, symbol2)
            return r
        
        def get_min_length(self, symbol1, symbol2):
            try:
                r = float(self.elems_param[symbol1 + '-' + symbol2][0]) 
            except:
                pass
            try:
                r = float(self.elems_param[symbol2 + '-' + symbol1][0]) 
            except :
                print('Warning : wrong bond pair (min) ', symbol1, symbol2)
            return r
        
        def get_maximum_length(self):
            length = []
            for _min, _max, flag in self.elems_param:                
                length.append(_max)            
            return max(length)
        
        def get_minimum_length(self):
            length = []
            for _min, _max, flag in self.elems_param:                
                length.append(_min)
            return min(length)

    class polyhedron(object):
        """
            edit polyhedron parameter class
        """
        def __init__(self):
            self.pairs = []
            self.opacity = 1.0
            self.color = [255,255,255]
    
    
    
    def __init__(self, path):
        self.atoms = None
        self.bonds = None
        self.polyhedron = None
        self.path = path
        self.is_exist = False
        
    def init_atoms(self):
        self.atoms = ParameterConfiguration.atoms()
    
    def init_bonds(self):
        self.bonds = ParameterConfiguration.bonds()
        
    def init_polyhedron(self):
        self.polyhedron = ParameterConfiguration.polyhedron()
    
    def read(self):
        config = ConfigObj(self.path)
        self.is_exist = True
        try:
            config['atoms'] # open check
            if self.atoms is None:
                self.atoms = ParameterConfiguration.atoms()
            for pair, color in config['atoms']['color'].items():
                self.atoms.elems_color[pair] = [int(c) for c in color]
            self.atoms.opacity = float(config['atoms']['opacity'])
            self.atoms.scale_factor = float(config['atoms']['factor'])            
        except Exception as e:
            print('param ini. atoms : ', e)        
        try:
            config['bonds']
            if self.bonds is None:
                self.bonds = ParameterConfiguration.bonds()
            for pair, params in config['bonds']['elems'].items():
                check = True if strtobool(params[2]) > 0 else False
                self.bonds.elems_param[pair] = [float(params[0]),float(params[1]),check]                
            self.bonds.opacity = float(config['bonds']['opacity'])
            self.bonds.scale_factor = float(config['bonds']['factor'])
            self.bonds.color_mode = True if strtobool(config['bonds']['mode']) == 1 else False
            self.bonds.color = [int(c) for c in config['bonds']['color']]
        except Exception as e:
            print('param ini. bonds : ', e)        
        
        try:
            config['poly']
            if self.polyhedron is None:
                self.polyhedron = ParameterConfiguration.polyhedron()
            self.polyhedron.pairs = []
            for pair in config['poly']['pairs']:
                self.polyhedron.pairs.append(pair)            
            self.polyhedron.opacity = float(config['poly']['opacity'])
            self.polyhedron.color = [int(c) for c in config['poly']['color']]
        except Exception as e:
            print('param ini. poly : ', e)        
            
    
    def write(self):        
        config = ConfigObj(self.path)
        config['atoms'] = {}
        config['atoms']['color'] = {}
        for pair, color in self.atoms.elems_color.items():
            config['atoms']['color'][pair] = color
        config['atoms']['opacity'] = self.atoms.opacity
        config['atoms']['factor'] = self.atoms.scale_factor
        
        config['bonds'] = {}
        config['bonds']['elems'] = {}
        for pair, params in self.bonds.elems_param.items():
            config['bonds']['elems'][pair] = params
        config['bonds']['opacity'] = self.bonds.opacity
        config['bonds']['factor'] = self.bonds.scale_factor
        config['bonds']['mode'] = self.bonds.color_mode
        config['bonds']['color'] = self.bonds.color
        
        config['poly'] = {}
        config['poly']['pairs'] = self.polyhedron.pairs
        config['poly']['opacity'] = self.polyhedron.opacity
        config['poly']['color'] = self.polyhedron.color
        config.write()