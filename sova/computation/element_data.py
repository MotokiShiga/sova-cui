# -*- coding: utf-8 -*-

import os
import numpy as np
import math

try:    
    from _path import basedir
except ImportError:    
    basedir = os.path.dirname(__file__)    

class ElementData(object):
    """
    Element Data class
    """
    def __init__(self, number, symbol, name, mass=0, spin='-'):
        self.number = number
        self.symbol = symbol
        self.name = name
        self.mass = mass
        self.spin = spin

class ElementsDict(object):
    """
    Ordered dict of Elements with lookup by number, symbol, and name.
    """
    def __init__(self, *elements):
        self._dict = {}
        self._list = []
        for element in elements:
            if element.number > len(self._list) + 1:
                raise ValueError("Elements must be added in order")
            if element.number <= len(self._list):
                self._list[element.number - 1] = element
            else:
                self._list.append(element)
            self._dict[element.number] = element
            self._dict[element.symbol] = element
            self._dict[element.symbol.upper()] = element
            self._dict[element.name] = element
            
    def __getitem__(self, key):
        try:
            return self._dict[key]
        except KeyError:
            try:
                start, stop, step = key.indices(len(self._list))
                return self._list[slice(start - 1, stop - 1, step)]
            except:
                raise KeyError


class XrayData(ElementData):
    """
    Xray Data class
    """
    def __init__(self, number, symbol, name):
        super(XrayData, self).__init__(number, symbol, name)
        self.a = np.zeros(5)
        self.b = np.zeros(5)
        self.c = 0.0


class XrayDict(ElementsDict):
    """
    Xray Dictionary with Xray Data
    """
    def __init__(self, *elements):
        super(XrayDict, self).__init__(*elements)

    def read(self, path):
        self.path = path
        try:
            f = open(self.path, "r")
            lines = f.readlines()
            for line in lines:
                if line[0] == '#': continue
                v = line.split()
                atom = v[0].upper()
                xray = self._dict[atom]
                for i in range(5):
                    xray.a[i] = float(v[2*i+1])
                    xray.b[i] = float(v[2*(i+1)])
                xray.c = float(v[-1])
            f.close()
        except IOError:
            raise

Xray = XrayDict(
    XrayData(1, 'H', 'Hydrogen'),
    XrayData(2, 'He', 'Helium'),
    XrayData(3, 'Li', 'Lithium'),
    XrayData(4, 'Be', 'Beryllium'),
    XrayData(5, 'B', 'Boron'),
    XrayData(6, 'C', 'Carbon'),
    XrayData(7, 'N', 'Nitrogen'),
    XrayData(8, 'O', 'Oxygen'),
    XrayData(9, 'F', 'Fluorine'),
    XrayData(10, 'Ne', 'Neon'),
    XrayData(11, 'Na', 'Sodium'),
    XrayData(12, 'Mg', 'Magnesium'),
    XrayData(13, 'Al', 'Aluminium'),
    XrayData(14, 'Si', 'Silicon'),
    XrayData(15, 'P', 'Phosphorus'),
    XrayData(16, 'S', 'Sulfur'),
    XrayData(17, 'Cl', 'Chlorine'),
    XrayData(18, 'Ar', 'Argon'),
    XrayData(19, 'K', 'Potassium'),
    XrayData(20, 'Ca', 'Calcium'),
    XrayData(21, 'Sc', 'Scandium'),
    XrayData(22, 'Ti', 'Titanium'),
    XrayData(23, 'V', 'Vanadium'),
    XrayData(24, 'Cr', 'Chromium'),
    XrayData(25, 'Mn', 'Manganese'),
    XrayData(26, 'Fe', 'Iron'),
    XrayData(27, 'Co', 'Cobalt'),
    XrayData(28, 'Ni', 'Nickel'),
    XrayData(29, 'Cu', 'Copper'),
    XrayData(30, 'Zn', 'Zinc'),
    XrayData(31, 'Ga', 'Gallium'),
    XrayData(32, 'Ge', 'Germanium'),
    XrayData(33, 'As', 'Arsenic'),
    XrayData(34, 'Se', 'Selenium'),
    XrayData(35, 'Br', 'Bromine'),
    XrayData(36, 'Kr', 'Krypton'),
    XrayData(37, 'Rb', 'Rubidium'),
    XrayData(38, 'Sr', 'Strontium'),
    XrayData(39, 'Y', 'Yttrium'),
    XrayData(40, 'Zr', 'Zirconium'),
    XrayData(41, 'Nb', 'Niobium'),
    XrayData(42, 'Mo', 'Molybdenum'),
    XrayData(43, 'Tc', 'Technetium'),
    XrayData(44, 'Ru', 'Ruthenium'),
    XrayData(45, 'Rh', 'Rhodium'),
    XrayData(46, 'Pd', 'Palladium'),
    XrayData(47, 'Ag', 'Silver'),
    XrayData(48, 'Cd', 'Cadmium'),
    XrayData(49, 'In', 'Indium'),
    XrayData(50, 'Sn', 'Tin'),
    XrayData(51, 'Sb', 'Antimony'),
    XrayData(52, 'Te', 'Tellurium'),
    XrayData(53, 'I', 'Iodine'),
    XrayData(54, 'Xe', 'Xenon'),
    XrayData(55, 'Cs', 'Caesium'),
    XrayData(56, 'Ba', 'Barium'),
    XrayData(57, 'La', 'Lanthanum'),
    XrayData(58, 'Ce', 'Cerium'),
    XrayData(59, 'Pr', 'Praseodymium'),
    XrayData(60, 'Nd', 'Neodymium'),
    XrayData(61, 'Pm', 'Promethium'),
    XrayData(62, 'Sm', 'Samarium'),
    XrayData(63, 'Eu', 'Europium'),
    XrayData(64, 'Gd', 'Gadolinium'),
    XrayData(65, 'Tb', 'Terbium'),
    XrayData(66, 'Dy', 'Dysprosium'),
    XrayData(67, 'Ho', 'Holmium'),
    XrayData(68, 'Er', 'Erbium'),
    XrayData(69, 'Tm', 'Thulium'),
    XrayData(70, 'Yb', 'Ytterbium'),
    XrayData(71, 'Lu', 'Lutetium'),
    XrayData(72, 'Hf', 'Hafnium'),
    XrayData(73, 'Ta', 'Tantalum'),
    XrayData(74, 'W', 'Tungsten'),
    XrayData(75, 'Re', 'Rhenium'),
    XrayData(76, 'Os', 'Osmium'),
    XrayData(77, 'Ir', 'Iridium'),
    XrayData(78, 'Pt', 'Platinum'),
    XrayData(79, 'Au', 'Gold'),
    XrayData(80, 'Hg', 'Mercury'),
    XrayData(81, 'Tl', 'Thallium'),
    XrayData(82, 'Pb', 'Lead'),
    XrayData(83, 'Bi', 'Bismuth'),
    XrayData(84, 'Po', 'Polonium'),
    XrayData(85, 'At', 'Astatine'),
    XrayData(86, 'Rn', 'Radon'),
    XrayData(87, 'Fr', 'Francium'),
    XrayData(88, 'Ra', 'Radium'),
    XrayData(89, 'Ac', 'Actinium'),
    XrayData(90, 'Th', 'Thorium'),
    XrayData(91, 'Pa', 'Protactinium'),
    XrayData(92, 'U', 'Uranium'),
    XrayData(93, 'Np', 'Neptunium'),
    XrayData(94, 'Pu', 'Plutonium'),
    XrayData(95, 'Am', 'Americium'),
    XrayData(96, 'Cm', 'Curium'),
    XrayData(97, 'Bk', 'Berkelium'),
    XrayData(98, 'Cf', 'Californium'),
    XrayData(99, 'Es', 'Einsteinium'),
    XrayData(100, 'Fm', 'Fermium'),
    XrayData(101, 'Md', 'Mendelevium'),
    XrayData(102, 'No', 'Nobelium'),
    XrayData(103, 'Lr', 'Lawrencium'),
    XrayData(104, 'Rf', 'Rutherfordium'),
    XrayData(105, 'Db', 'Dubnium'),
    XrayData(106, 'Sg', 'Seaborgium'),
    XrayData(107, 'Bh', 'Bohrium'),
    XrayData(108, 'Hs', 'Hassium'),
    XrayData(109, 'Mt', 'Meitnerium')
)

Xray.read(os.path.join(basedir, "xparm_ff.txt"))

class NeutronData(ElementData):
    """
    Neutron Data class
    """
    def __init__(self, number, symbol, name, mass=0, spin='-', complex='Re'):
        super(NeutronData, self).__init__(number, symbol, name, mass, spin)
        self.complex = complex
        self.bc = 0.0
        self.S_t = 0.0
        self.S_a = 0.0
        self.c = 0.0
        self.d = 0.0
        self.e = 0.0

class NeutronDict(ElementsDict):
    """
    Neutron Dictionary with Neutron Data
    """
    def __init__(self, *elements):
        #super(NeutronDict, self).__init__(*elements)
        self._dict = {}
        self._list = []
        for element in elements:
            self._list.append(element)
            key = element.symbol.upper() + '/' + \
                  str(abs(element.mass)) + '/' + \
                  element.complex
            self._dict[key] = element

    def append(self, *elements):
        for element in elements:
            self._list.append(element)
            key = element.symbol.upper() + '/' + \
                  str(abs(element.mass)) + '/' + \
                  element.complex
            self._dict[key] = element
    
    def read(self, path):
        self.path = path
        try:
            f = open(self.path, "r")
            lines = f.readlines()
            
            for line in lines:
                v = line.split()
                complex = 'Re'
                if float(v[-1]) == 0.0 and float(v[-2]) == 0.0 and \
                   float(v[-3]) == 0.0 and float(v[-4]) == 0.0:
                    complex = 'Im'

                mass = int(v[1])
                key = v[0].upper() + '/' + str(abs(mass)) + '/' + complex
                neutron = self._dict[key]
                
                neutron.bc = float(v[2])
                neutron.S_t = float(v[-2])
                neutron.S_a = float(v[-1])

            f.close()

        except IOError:
            raise

    def first(self, symbol):
        for neutron in self._list:
            if neutron.symbol.upper() == symbol.upper():
                return neutron

    def __getitem__(self, key):
        try:
            return self._dict[key]
        except KeyError:
            try:
                start, stop, step = key.indices(len(self._list))
                return self._list[slice(start - 1, stop - 1, step)]
            except:
                raise KeyError

Neutron = NeutronDict(
    NeutronData(1, 'H', 'Hydrogen'),
    NeutronData(1, 'H', 'Hydrogen', 1, '1/2(+)'),
    NeutronData(1, 'H', 'Hydrogen', 2, '1(+)'),
    NeutronData(1, 'H', 'Hydrogen', 3, '1/2(+)'),
    NeutronData(2, 'He', 'Helium'),
    NeutronData(2, 'He', 'Helium', 3, '1/2(+)'),
    NeutronData(2, 'He', 'Helium', 3, '1/2(+)', 'Im'),
    NeutronData(2, 'He', 'Helium', 4, '0(+)'),
    NeutronData(3, 'Li', 'Lithium'),
    NeutronData(3, 'Li', 'Lithium', 6, '1(+)'),
    NeutronData(3, 'Li', 'Lithium', 6, '1(+)', 'Im'),
    NeutronData(3, 'Li', 'Lithium', 7, '3/2(-)'),
    NeutronData(4, 'Be', 'Beryllium', 9, '3/2(-)'),
    NeutronData(5, 'B', 'Boron'),
    NeutronData(5, 'B', 'Boron', complex='Im'),
    NeutronData(5, 'B', 'Boron', 10, '3(+)'),
    NeutronData(5, 'B', 'Boron', 10, '3(+)', 'Im'),
    NeutronData(5, 'B', 'Boron', 11, '3/2(-)'),    
    NeutronData(6, 'C', 'Carbon', 0, '-'),
    NeutronData(6, 'C', 'Carbon', 12, '0(+)'),
    NeutronData(6, 'C', 'Carbon', 13, '1/2(-)'),
    NeutronData(7, 'N', 'Nitrogen'),
    NeutronData(7, 'N', 'Nitrogen', 14, '1(+)'),
    NeutronData(7, 'N', 'Nitrogen', 15, '1/2(-)'),
    NeutronData(8, 'O', 'Oxygen', 0, '-'),
    NeutronData(8, 'O', 'Oxygen', 16, '0(+)'),
    NeutronData(8, 'O', 'Oxygen', 17, '5/2(+)'),
    NeutronData(8, 'O', 'Oxygen', 18, '0(+)'),
    NeutronData(9, 'F', 'Fluorine', 19, '1/2(+)'),
    NeutronData(10, 'Ne', 'Neon', 0, '-'),
    NeutronData(10, 'Ne', 'Neon', 20, '0(+)'),
    NeutronData(10, 'Ne', 'Neon', 21, '3/2(+)'),
    NeutronData(10, 'Ne', 'Neon', 22, '0(+)'),
    NeutronData(11, 'Na', 'Sodium', 23, '3/2(+)'),
    NeutronData(12, 'Mg', 'Magnesium', 0),
    NeutronData(12, 'Mg', 'Magnesium', 24, '0(+)'),
    NeutronData(12, 'Mg', 'Magnesium', 25, '5/2(+)'),
    NeutronData(12, 'Mg', 'Magnesium', 26, '0(+)'),
    NeutronData(13, 'Al', 'Aluminium', 27, '5/2(+)'),
    NeutronData(14, 'Si', 'Silicon', 0),
    NeutronData(14, 'Si', 'Silicon', 28, '0(+)'),
    NeutronData(14, 'Si', 'Silicon', 29, '1/2(+)'),
    NeutronData(14, 'Si', 'Silicon', 30, '0(+)'),
    NeutronData(15, 'P', 'Phosphorus', 31, '1/2(+)'),
    NeutronData(16, 'S', 'Sulfur', 0),
    NeutronData(16, 'S', 'Sulfur', 32, '0(+)'),
    NeutronData(16, 'S', 'Sulfur', 33, '3/2(+)'),
    NeutronData(16, 'S', 'Sulfur', 34, '0(+)'),
    NeutronData(16, 'S', 'Sulfur', 36, '0(+)'),
    NeutronData(17, 'Cl', 'Chlorine', 0),
    NeutronData(17, 'Cl', 'Chlorine', 35, '3/2(+)'),
    NeutronData(17, 'Cl', 'Chlorine', 37, '3/2(+)'),
    NeutronData(18, 'Ar', 'Argon', 0),
    NeutronData(18, 'Ar', 'Argon', 36, '0(+)'),
    NeutronData(18, 'Ar', 'Argon', 38, '0(+)'),
    NeutronData(18, 'Ar', 'Argon', 40, '0(+)'),
    NeutronData(19, 'K', 'Potassium'),
    NeutronData(19, 'K', 'Potassium', 39, '3/2(+)'),
    NeutronData(19, 'K', 'Potassium', 40, '4(-)'),
    NeutronData(19, 'K', 'Potassium', 41, '3/2(+)'),
    NeutronData(20, 'Ca', 'Calcium'),
    NeutronData(20, 'Ca', 'Calcium', 40, '0(+)'),
    NeutronData(20, 'Ca', 'Calcium', 42, '0(+)'),
    NeutronData(20, 'Ca', 'Calcium', 43, '7/2(-)'),
    NeutronData(20, 'Ca', 'Calcium', 44, '0(+)'),
    NeutronData(20, 'Ca', 'Calcium', 46, '0(+)'),
    NeutronData(20, 'Ca', 'Calcium', 48, '0(+)'),
    NeutronData(21, 'Sc', 'Scandium', 45, '7/2(-)'),
    NeutronData(22, 'Ti', 'Titanium'),
    NeutronData(22, 'Ti', 'Titanium', 46, '0(+)'),
    NeutronData(22, 'Ti', 'Titanium', 47, '5/2(-)'),
    NeutronData(22, 'Ti', 'Titanium', 48, '0(+)'),
    NeutronData(22, 'Ti', 'Titanium', 49, '7/2(-)'),
    NeutronData(22, 'Ti', 'Titanium', 50, '0(+)'),
    NeutronData(23, 'V', 'Vanadium', ),
    NeutronData(23, 'V', 'Vanadium', 50, '6(+)'),
    NeutronData(23, 'V', 'Vanadium', 51, '7/2(-)'),
    NeutronData(24, 'Cr', 'Chromium'),
    NeutronData(24, 'Cr', 'Chromium', 50, '0(+)'),
    NeutronData(24, 'Cr', 'Chromium', 52, '0(+)'),
    NeutronData(24, 'Cr', 'Chromium', 53, '3/2(-)'),
    NeutronData(24, 'Cr', 'Chromium', 54, '0(+)'),
    NeutronData(25, 'Mn', 'Manganese', 55, '5/2(-)'),
    NeutronData(26, 'Fe', 'Iron'),
    NeutronData(26, 'Fe', 'Iron', 54, '0(+)'),
    NeutronData(26, 'Fe', 'Iron', 56, '0(+)'),
    NeutronData(26, 'Fe', 'Iron', 57, '1/2(-)'),
    NeutronData(26, 'Fe', 'Iron', 58, '0(+)'),
    NeutronData(27, 'Co', 'Cobalt', 59, '7/2(-)'),
    NeutronData(28, 'Ni', 'Nickel'),
    NeutronData(28, 'Ni', 'Nickel', 58, '0(+)'),
    NeutronData(28, 'Ni', 'Nickel', 60, '0(+)'),
    NeutronData(28, 'Ni', 'Nickel', 61, '3/2(-)'),
    NeutronData(28, 'Ni', 'Nickel', 62, '0(+)'),
    NeutronData(28, 'Ni', 'Nickel', 64, '0(+)'),
    NeutronData(29, 'Cu', 'Copper'),
    NeutronData(29, 'Cu', 'Copper', 63, '3/2(-)'),
    NeutronData(29, 'Cu', 'Copper', 65, '3/2(-)'),
    NeutronData(30, 'Zn', 'Zinc'),
    NeutronData(30, 'Zn', 'Zinc', 64, '0(+)'),
    NeutronData(30, 'Zn', 'Zinc', 66, '0(+)'),
    NeutronData(30, 'Zn', 'Zinc', 67, '5/2(-)'),
    NeutronData(30, 'Zn', 'Zinc', 68, '0(+)'),
    NeutronData(30, 'Zn', 'Zinc', 70, '0(+)'),
    NeutronData(31, 'Ga', 'Gallium'),
    NeutronData(31, 'Ga', 'Gallium', 69, '3/2(-)'),
    NeutronData(31, 'Ga', 'Gallium', 71, '3/2(-)'),
    NeutronData(32, 'Ge', 'Germanium', ),
    NeutronData(32, 'Ge', 'Germanium', 70, '0(+)'),
    NeutronData(32, 'Ge', 'Germanium', 72, '0(+)'),
    NeutronData(32, 'Ge', 'Germanium', 73, '9/2(+)'),
    NeutronData(32, 'Ge', 'Germanium', 74, '0(+)'),
    NeutronData(32, 'Ge', 'Germanium', 76, '0(+)'),
    NeutronData(33, 'As', 'Arsenic', 75, '3/2(-)'),
    NeutronData(34, 'Se', 'Selenium'),
    NeutronData(34, 'Se', 'Selenium', 74, '0(+)'),
    NeutronData(34, 'Se', 'Selenium', 76, '0(+)'),
    NeutronData(34, 'Se', 'Selenium', 77, '1/2(-)'),
    NeutronData(34, 'Se', 'Selenium', 78, '0(+)'),
    NeutronData(34, 'Se', 'Selenium', 80, '0(+)'),
    NeutronData(34, 'Se', 'Selenium', 82, '0(+)'),
    NeutronData(35, 'Br', 'Bromine'),
    NeutronData(35, 'Br', 'Bromine', 79, '3/2(-)'),
    NeutronData(35, 'Br', 'Bromine', 81, '3/2(-)'),
    NeutronData(36, 'Kr', 'Krypton'),
    NeutronData(36, 'Kr', 'Krypton', 78, '0(+)'),
    NeutronData(36, 'Kr', 'Krypton', 80, '0(+)'),
    NeutronData(36, 'Kr', 'Krypton', 82, '0(+)'),
    NeutronData(36, 'Kr', 'Krypton', 83, '9/2(+)'),
    NeutronData(36, 'Kr', 'Krypton', 84, '0(+)'),
    NeutronData(36, 'Kr', 'Krypton', 86, '0(+)'),
    NeutronData(37, 'Rb', 'Rubidium'),
    NeutronData(37, 'Rb', 'Rubidium', 85, '5/2(-)'),
    NeutronData(37, 'Rb', 'Rubidium', 87, '3/2(-)'),
    NeutronData(38, 'Sr', 'Strontium'),
    NeutronData(38, 'Sr', 'Strontium', 84, '0(+)'),
    NeutronData(38, 'Sr', 'Strontium', 86, '0(+)'),
    NeutronData(38, 'Sr', 'Strontium', 87, '9/2(+)'),
    NeutronData(38, 'Sr', 'Strontium', 88, '0(+)'),
    NeutronData(39, 'Y', 'Yttrium', 89, '1/2(-)'),
    NeutronData(40, 'Zr', 'Zirconium'),
    NeutronData(40, 'Zr', 'Zirconium', 90, '0(+)'),
    NeutronData(40, 'Zr', 'Zirconium', 91, '5/2(+)'),
    NeutronData(40, 'Zr', 'Zirconium', 92, '0(+)'),
    NeutronData(40, 'Zr', 'Zirconium', 94, '0(+)'),
    NeutronData(40, 'Zr', 'Zirconium', 96, '0(+)'),
    NeutronData(41, 'Nb', 'Niobium', 93, '9/2(+)'),
    NeutronData(42, 'Mo', 'Molybdenum'),
    NeutronData(42, 'Mo', 'Molybdenum', 92, '0(+)'),
    NeutronData(42, 'Mo', 'Molybdenum', 94, '0(+)'),
    NeutronData(42, 'Mo', 'Molybdenum', 95, '5/2(+)'),
    NeutronData(42, 'Mo', 'Molybdenum', 96, '0(+)'),
    NeutronData(42, 'Mo', 'Molybdenum', 97, '5/2(+)'),
    NeutronData(42, 'Mo', 'Molybdenum', 98, '0(+)'),
    NeutronData(42, 'Mo', 'Molybdenum', 100, '0(+)'),
    NeutronData(43, 'Tc', 'Technetium', 99, '9/2(+)'),
    NeutronData(44, 'Ru', 'Ruthenium'),
    NeutronData(44, 'Ru', 'Ruthenium', 96, '0(+)'),
    NeutronData(44, 'Ru', 'Ruthenium', 98, '0(+)'),
    NeutronData(44, 'Ru', 'Ruthenium', 99, '5/2(+)'),
    NeutronData(44, 'Ru', 'Ruthenium', 100, '0(+)'),
    NeutronData(44, 'Ru', 'Ruthenium', 101, '5/2(+)'),
    NeutronData(44, 'Ru', 'Ruthenium', 102, '0(+)'),
    NeutronData(44, 'Ru', 'Ruthenium', 104, '0(+)'),
    NeutronData(45, 'Rh', 'Rhodium', 103, '1/2(-)'),
    NeutronData(46, 'Pd', 'Palladium'),
    NeutronData(46, 'Pd', 'Palladium', 102, '0(+)'),
    NeutronData(46, 'Pd', 'Palladium', 104, '0(+)'),
    NeutronData(46, 'Pd', 'Palladium', 105, '5/2(+)'),
    NeutronData(46, 'Pd', 'Palladium', 106, '0(+)'),
    NeutronData(46, 'Pd', 'Palladium', 108, '0(+)'),
    NeutronData(46, 'Pd', 'Palladium', 110, '0(+)'),
    NeutronData(47, 'Ag', 'Silver'),
    NeutronData(47, 'Ag', 'Silver', 107, '1/2(-)'),
    NeutronData(47, 'Ag', 'Silver', 109, '1/2(-)'),
    NeutronData(48, 'Cd', 'Cadmium'),
    NeutronData(48, 'Cd', 'Cadmium', complex='Im'),
    NeutronData(48, 'Cd', 'Cadmium', 106, '0(+)'),
    NeutronData(48, 'Cd', 'Cadmium', 108, '0(+)'),
    NeutronData(48, 'Cd', 'Cadmium', 110, '0(+)'),
    NeutronData(48, 'Cd', 'Cadmium', 111, '1/2(+)'),
    NeutronData(48, 'Cd', 'Cadmium', 112, '0(+)'),
    NeutronData(48, 'Cd', 'Cadmium', 113, '1/2(+)'),
    NeutronData(48, 'Cd', 'Cadmium', 113, '1/2(+)', 'Im'),
    NeutronData(48, 'Cd', 'Cadmium', 114, '0(+)'),
    NeutronData(48, 'Cd', 'Cadmium', 116, '0(+)'),
    NeutronData(49, 'In', 'Indium'),
    NeutronData(49, 'In', 'Indium', complex='Im'),
    NeutronData(49, 'In', 'Indium', 113, '9/2(+)'),
    NeutronData(49, 'In', 'Indium', 115, '9/2(+)'),
    NeutronData(49, 'In', 'Indium', 115, '9/2(+)', 'Im'),
    NeutronData(50, 'Sn', 'Tin'),
    NeutronData(50, 'Sn', 'Tin', 112, '0(+)'),
    NeutronData(50, 'Sn', 'Tin', 114, '0(+)'),
    NeutronData(50, 'Sn', 'Tin', 115, '1/2(+)'),
    NeutronData(50, 'Sn', 'Tin', 116, '0(+)'),
    NeutronData(50, 'Sn', 'Tin', 117, '1/2(+)'),
    NeutronData(50, 'Sn', 'Tin', 118, '0(+)'),
    NeutronData(50, 'Sn', 'Tin', 119, '1/2(+)'),
    NeutronData(50, 'Sn', 'Tin', 120, '0(+)'),
    NeutronData(50, 'Sn', 'Tin', 122, '0(+)'),
    NeutronData(50, 'Sn', 'Tin', 124, '0(+)')
)

Neutron.append(
    NeutronData(51, 'Sb', 'Antimony'),
    NeutronData(51, 'Sb', 'Antimony', 121, '7/2(+)'),
    NeutronData(51, 'Sb', 'Antimony', 123, '5/2(+)'),
    NeutronData(52, 'Te', 'Tellurium'),
    NeutronData(52, 'Te', 'Tellurium', 120, '0(+)'),
    NeutronData(52, 'Te', 'Tellurium', 122, '0(+)'),
    NeutronData(52, 'Te', 'Tellurium', 123, '1/2(+)'),
    NeutronData(52, 'Te', 'Tellurium', 123, '1/2(+)', 'Im'),
    NeutronData(52, 'Te', 'Tellurium', 124, '0(+)'),
    NeutronData(52, 'Te', 'Tellurium', 125, '1/2(+)'),
    NeutronData(52, 'Te', 'Tellurium', 126, '0(+)'),
    NeutronData(52, 'Te', 'Tellurium', 128, '0(+)'),
    NeutronData(52, 'Te', 'Tellurium', 130, '0(+)'),
    NeutronData(53, 'I', 'Iodine', 127, '5/2(+)'),
    NeutronData(54, 'Xe', 'Xenon'),
    NeutronData(54, 'Xe', 'Xenon', 124, '0(+)'),
    NeutronData(54, 'Xe', 'Xenon', 126, '0(+)'),
    NeutronData(54, 'Xe', 'Xenon', 128, '0(+)'),
    NeutronData(54, 'Xe', 'Xenon', 129, '0(+)'),
    NeutronData(54, 'Xe', 'Xenon', 130, '0(+)'),
    NeutronData(54, 'Xe', 'Xenon', 131, '3/2(+)'),
    NeutronData(54, 'Xe', 'Xenon', 132, '0(+)'),
    NeutronData(54, 'Xe', 'Xenon', 134, '0(+)'),
    NeutronData(54, 'Xe', 'Xenon', 136, '0(+)'),
    NeutronData(55, 'Cs', 'Caesium', 133, '7/2(+)'),
    NeutronData(56, 'Ba', 'Barium'),
    NeutronData(56, 'Ba', 'Barium', 130, '0(+)'),
    NeutronData(56, 'Ba', 'Barium', 132, '0(+)'),
    NeutronData(56, 'Ba', 'Barium', 134, '0(+)'),
    NeutronData(56, 'Ba', 'Barium', 135, '3/2(+)'),
    NeutronData(56, 'Ba', 'Barium', 136, '0(+)'),
    NeutronData(56, 'Ba', 'Barium', 137, '3/2(+)'),
    NeutronData(56, 'Ba', 'Barium', 138, '0(+)'),
    NeutronData(57, 'La', 'Lanthanum'),
    NeutronData(57, 'La', 'Lanthanum', 138, '5(+)'),
    NeutronData(57, 'La', 'Lanthanum', 139, '7/2(+)'),
    NeutronData(58, 'Ce', 'Cerium'),
    NeutronData(58, 'Ce', 'Cerium', 136, '0(+)'),
    NeutronData(58, 'Ce', 'Cerium', 138, '0(+)'),
    NeutronData(58, 'Ce', 'Cerium', 140, '0(+)'),
    NeutronData(58, 'Ce', 'Cerium', 142, '0(+)'),
    NeutronData(59, 'Pr', 'Praseodymium', 141, '5/2(+)'),
    NeutronData(60, 'Nd', 'Neodymium'),
    NeutronData(60, 'Nd', 'Neodymium', 142, '0(+)'),
    NeutronData(60, 'Nd', 'Neodymium', 143, '7/2(-)'),
    NeutronData(60, 'Nd', 'Neodymium', 144, '0(+)'),
    NeutronData(60, 'Nd', 'Neodymium', 145, '7/2(-)'),
    NeutronData(60, 'Nd', 'Neodymium', 146, '0(+)'),
    NeutronData(60, 'Nd', 'Neodymium', 148, '0(+)'),
    NeutronData(60, 'Nd', 'Neodymium', 150, '0(+)'),
    NeutronData(61, 'Pm', 'Promethium', 147, '7/2(+)'),
    NeutronData(62, 'Sm', 'Samarium'),
    NeutronData(62, 'Sm', 'Samarium', complex='Im'),
    NeutronData(62, 'Sm', 'Samarium', 144, '0(+)'),
    NeutronData(62, 'Sm', 'Samarium', 147, '7/2(-)'),
    NeutronData(62, 'Sm', 'Samarium', 148, '0(+)'),
    NeutronData(62, 'Sm', 'Samarium', 149, '7/2(-)'),
    NeutronData(62, 'Sm', 'Samarium', 149, '7/2(-)', 'Im'),
    NeutronData(62, 'Sm', 'Samarium', 150, '0(+)'),
    NeutronData(62, 'Sm', 'Samarium', 152, '0(+)'),
    NeutronData(62, 'Sm', 'Samarium', 154, '0(+)'),
    NeutronData(63, 'Eu', 'Europium'),
    NeutronData(63, 'Eu', 'Europium', complex='Im'),
    NeutronData(63, 'Eu', 'Europium', 151, '5/2(+)'),
    NeutronData(63, 'Eu', 'Europium', 151, '5/2(+)', 'Im'),
    NeutronData(63, 'Eu', 'Europium', 153, '5/2(+)'),
    NeutronData(64, 'Gd', 'Gadolinium'),
    NeutronData(64, 'Gd', 'Gadolinium', complex='Im'),
    NeutronData(64, 'Gd', 'Gadolinium', 152, '0(+)'),
    NeutronData(64, 'Gd', 'Gadolinium', 154, '0(+)'),
    NeutronData(64, 'Gd', 'Gadolinium', 155, '3/2(-)'),
    NeutronData(64, 'Gd', 'Gadolinium', 155, '3/2(-)', 'Im'),
    NeutronData(64, 'Gd', 'Gadolinium', 156, '0(+)'),
    NeutronData(64, 'Gd', 'Gadolinium', 157, '3/2(-)'),
    NeutronData(64, 'Gd', 'Gadolinium', 157, '3/2(-)', 'Im'),
    NeutronData(64, 'Gd', 'Gadolinium', 158, '0(+)'),
    NeutronData(64, 'Gd', 'Gadolinium', 160, '0(+)'),
    NeutronData(65, 'Tb', 'Terbium', 159, '3/2(+)'),
    NeutronData(66, 'Dy', 'Dysprosium'),
    NeutronData(66, 'Dy', 'Dysprosium', complex='Im'),
    NeutronData(66, 'Dy', 'Dysprosium', 156, '0(+)'),
    NeutronData(66, 'Dy', 'Dysprosium', 158, '0(+)'),
    NeutronData(66, 'Dy', 'Dysprosium', 160, '0(+)'),
    NeutronData(66, 'Dy', 'Dysprosium', 161, '5/2(+)'),
    NeutronData(66, 'Dy', 'Dysprosium', 162, '0(+)'),
    NeutronData(66, 'Dy', 'Dysprosium', 163, '5/2(+)'),
    NeutronData(66, 'Dy', 'Dysprosium', 164, '0(+)'),
    NeutronData(66, 'Dy', 'Dysprosium', 164, '0(+)', 'Im'),
    NeutronData(67, 'Ho', 'Holmium', 165, '7/2(-)'),
    NeutronData(68, 'Er', 'Erbium'),
    NeutronData(68, 'Er', 'Erbium', 162, '0(+)'),
    NeutronData(68, 'Er', 'Erbium', 164, '0(+)'),
    NeutronData(68, 'Er', 'Erbium', 166, '0(+)'),
    NeutronData(68, 'Er', 'Erbium', 167, '7/2(+)'),
    NeutronData(68, 'Er', 'Erbium', 168, '0(+)'),
    NeutronData(68, 'Er', 'Erbium', 170, '0(+)'),
    NeutronData(69, 'Tm', 'Thulium', 169, '1/2(+)'),
    NeutronData(70, 'Yb', 'Ytterbium'),
    NeutronData(70, 'Yb', 'Ytterbium', 168, '0(+)'),
    NeutronData(70, 'Yb', 'Ytterbium', 168, '0(+)', 'Im'),
    NeutronData(70, 'Yb', 'Ytterbium', 170, '0(+)'),
    NeutronData(70, 'Yb', 'Ytterbium', 171, '1/2(-)'),
    NeutronData(70, 'Yb', 'Ytterbium', 172, '0(+)'),
    NeutronData(70, 'Yb', 'Ytterbium', 173, '5/2(-)'),
    NeutronData(70, 'Yb', 'Ytterbium', 174, '0(+)'),
    NeutronData(70, 'Yb', 'Ytterbium', 176, '0(+)'),
    NeutronData(71, 'Lu', 'Lutetium'),
    NeutronData(71, 'Lu', 'Lutetium', 175, '7/2(+)'),
    NeutronData(71, 'Lu', 'Lutetium', 176, '7(-)'),
    NeutronData(71, 'Lu', 'Lutetium', 176, '7(-)', 'Im'),
    NeutronData(72, 'Hf', 'Hafnium'),
    NeutronData(72, 'Hf', 'Hafnium', 174, '0(+)'),
    NeutronData(72, 'Hf', 'Hafnium', 176, '0(+)'),
    NeutronData(72, 'Hf', 'Hafnium', 177, '7/2(-)'),
    NeutronData(72, 'Hf', 'Hafnium', 178, '0(+)'),
    NeutronData(72, 'Hf', 'Hafnium', 179, '9/2(+)'),
    NeutronData(72, 'Hf', 'Hafnium', 180, '0(+)'),
    NeutronData(73, 'Ta', 'Tantalum'),
    NeutronData(73, 'Ta', 'Tantalum', 180, '9(-)'),
    NeutronData(73, 'Ta', 'Tantalum', 181, '7/2(+)'),
    NeutronData(74, 'W', 'Tungsten'),
    NeutronData(74, 'W', 'Tungsten', 180, '0(+)'),
    NeutronData(74, 'W', 'Tungsten', 182, '0(+)'),
    NeutronData(74, 'W', 'Tungsten', 183, '1/2(-)'),
    NeutronData(74, 'W', 'Tungsten', 184, '0(+)'),
    NeutronData(74, 'W', 'Tungsten', 186, '0(+)'),
    NeutronData(75, 'Re', 'Rhenium'),
    NeutronData(75, 'Re', 'Rhenium', 185, '5/2(+)'),
    NeutronData(75, 'Re', 'Rhenium', 187, '5/2(+)'),
    NeutronData(76, 'Os', 'Osmium'),
    NeutronData(76, 'Os', 'Osmium', 184, '0(+)'),
    NeutronData(76, 'Os', 'Osmium', 186, '0(+)'),
    NeutronData(76, 'Os', 'Osmium', 187, '1/2(-)'),
    NeutronData(76, 'Os', 'Osmium', 188, '0(+)'),
    NeutronData(76, 'Os', 'Osmium', 189, '3/2(-)'),
    NeutronData(76, 'Os', 'Osmium', 190, '0(+)'),
    NeutronData(76, 'Os', 'Osmium', 192, '0(+)'),
    NeutronData(77, 'Ir', 'Iridium'),
    NeutronData(77, 'Ir', 'Iridium', 191, '3/2(+)'),
    NeutronData(77, 'Ir', 'Iridium', 193, '3/2(+)'),
    NeutronData(78, 'Pt', 'Platinum'),
    NeutronData(78, 'Pt', 'Platinum', 190, '0(+)'),
    NeutronData(78, 'Pt', 'Platinum', 192, '0(+)'),
    NeutronData(78, 'Pt', 'Platinum', 194, '0(+)'),
    NeutronData(78, 'Pt', 'Platinum', 195, '1/2(-)'),
    NeutronData(78, 'Pt', 'Platinum', 196, '0(+)'),
    NeutronData(78, 'Pt', 'Platinum', 198, '0(+)'),
    NeutronData(79, 'Au', 'Gold', 197, '3/2(+)'),
    NeutronData(80, 'Hg', 'Mercury'),
    NeutronData(80, 'Hg', 'Mercury', 196, '0(+)'),
    NeutronData(80, 'Hg', 'Mercury', 198, '0(+)'),
    NeutronData(80, 'Hg', 'Mercury', 199, '1/2(-)'),
    NeutronData(80, 'Hg', 'Mercury', 200, '0(+)'),
    NeutronData(80, 'Hg', 'Mercury', 201, '3/2(-)'),
    NeutronData(80, 'Hg', 'Mercury', 202, '0(+)'),
    NeutronData(80, 'Hg', 'Mercury', 204, '0(+)'),
    NeutronData(81, 'Tl', 'Thallium'),
    NeutronData(81, 'Tl', 'Thallium', 203, '1/2(+)'),
    NeutronData(81, 'Tl', 'Thallium', 205, '1/2(+)'),
    NeutronData(82, 'Pb', 'Lead'),
    NeutronData(82, 'Pb', 'Lead', 204, '0(+)'),
    NeutronData(82, 'Pb', 'Lead', 206, '0(+)'),
    NeutronData(82, 'Pb', 'Lead', 207, '1/2(-)'),
    NeutronData(82, 'Pb', 'Lead', 208, '0(+)'),
    NeutronData(83, 'Bi', 'Bismuth', 209, '9/2(-)'),
    NeutronData(88, 'Ra', 'Radium', 226, '0(+)'),
    NeutronData(90, 'Th', 'Thorium', 232, '0(+)'),
    NeutronData(91, 'Pa', 'Protactinium', 231, '3/2(-)'),
    NeutronData(92, 'U', 'Uranium'),
    NeutronData(92, 'U', 'Uranium', 233, '5/2(+)'),
    NeutronData(92, 'U', 'Uranium', 234, '0(+)'),
    NeutronData(92, 'U', 'Uranium', 235, '7/2(-)'),
    NeutronData(92, 'U', 'Uranium', 238, '0(+)'),
    NeutronData(93, 'Np', 'Neptunium', 237, '5/2(+)'),
    NeutronData(94, 'Pu', 'Plutonium', 238, '0(+)'),
    NeutronData(94, 'Pu', 'Plutonium', 239, '1/2(+)'),
    NeutronData(94, 'Pu', 'Plutonium', 240, '0(+)'),
    NeutronData(94, 'Pu', 'Plutonium', 242, '0(+)'),
    NeutronData(95, 'Am', 'Americium', 243, '5/2(-)'),
    NeutronData(96, 'Cm', 'Curium', 244, '0(+)'),
    NeutronData(96, 'Cm', 'Curium', 246, '0(+)'),
    NeutronData(96, 'Cm', 'Curium', 248, '0(+)')
)

Neutron.read(os.path.join(basedir, "nparma_bb.txt"))


if __name__ == '__main__':
    #print Xray['Li'].c
    #Xray.read("./xparm_ff.txt")
    #print Xray['Ne'].a[0]
    
    #print Neutron['H'].number
    pass
