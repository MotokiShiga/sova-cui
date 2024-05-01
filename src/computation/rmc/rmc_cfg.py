# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 15:12:43 2023

@author: H. Morita
"""

import numpy as np
import collections

def write_rmc_config(path,title,ngen,ntried,nacc,nsaved,
                     truncated,vectors,ni,nsites,sites,
                     neuler,centres):

    # Writes a configuration data file
    nmol = np.sum(ni)
    nmoltypes = len(ni)
    nsites_max = 1    
    for i in range(nmoltypes):
       nsites_max = max(nsites_max, nsites[i])
    if type(centres) == list:
        nvar = len(centres)
    elif type(centres) == np.ndarray:
        nvar = centres.shape[0]
    
    with open(path, 'w') as f:
        f.write("(Version 3 format configuration file)\n")
        f.write(title + '\n')
        f.write('\n\n')            
        f.write(" ")
        f.write("{:>10d}".format(ngen))
        f.write("{:>10d}".format(ntried))
        f.write("{:>10d}".format(nacc))
        f.write(" moves generated, tried, accepted\n")
        f.write(" ")
        f.write("{:>10d}".format(nsaved))
        f.write("                     configurations saved\n\n")
        f.write(" ")
        f.write("{:>10d}".format(nmol))
        f.write(" molecules (of all types)\n")
        f.write(" ")
        f.write("{:>10d}".format(nmoltypes))
        f.write(" types of molecule\n")
        f.write(" ")
        f.write("{:>10d}".format(nsites_max))
        f.write(" is the largest number of atoms in a molecule\n")
        f.write(" ")
        f.write("{:>10d}".format(neuler))
        f.write(" Euler angles are provided\n\n")
        f.write("          ")
        if truncated == True:
            f.write('T')
            f.write(" (Box is truncated octahedral)\n")
        else:
            f.write('F')
            f.write(" (Box is not truncated octahedral)\n")
        f.write("           ")
        f.write(" Defining vectors are:\n")            
        for i in range(3):
            f.write("           ")
            f.write(' {:10.06f}'.format(vectors[i][0]))
            f.write(' {:10.06f}'.format(vectors[i][1]))
            f.write(' {:10.06f}\n'.format(vectors[i][2]))
        f.write('\n')

        for i in range(nmoltypes):
            f.write(" {:>10d}".format(ni[i]))
            f.write(" molecules of type ")
            f.write("{:>2d}\n".format(i+1))
            f.write(" {:>10d}".format(nsites[i]))                
            f.write(" atomic sites\n")
            f.write("           ")
            f.write(' {:10.06f}'.format(sites[i][0]))
            f.write(' {:10.06f}'.format(sites[i][1]))
            f.write(' {:10.06f}'.format(sites[i][2]))
            f.write('\n\n')

        for i in range(nmol):
            f.write(' {:10.07f}'.format(centres[0][i]))
            f.write(' {:10.07f}'.format(centres[1][i]))
            f.write(' {:10.07f}'.format(centres[2][i]))
            for j in range(nvar-3):
                f.write(' {:10.08f}'.format(centres[3+j][i]))                
            f.write('\n')
    
def write(path,title,vectors,elements,positions):    
    ngen = 0
    ntried = 0
    nacc = 0
    nsaved = 0
    ec = [k for k, v in collections.Counter(elements).items() if v > 0]    
    c = collections.Counter(elements)
    n = len(elements)
    nmol_types = len(ec)
    nsites_max = 1
    neuler = 0
    truncated = False
    ni = []
    for e in ec:
        ni.append(c[e])
    nsites = np.ones(nmol_types, dtype=int)
    sites = np.zeros((nmol_types,3))
    matrix = vectors/2
    minv = np.linalg.inv(matrix)    
    
    with open(path, mode='w') as f:
        f.write("(Version 3 format configuration file)\n")
        f.write(title + '\n')
        f.write('\n\n')            
        f.write(" ")
        f.write("{:>10d}".format(ngen))
        f.write("{:>10d}".format(ntried))
        f.write("{:>10d}".format(nacc))
        f.write(" moves generated, tried, accepted\n")
        f.write(" ")
        f.write("{:>10d}".format(nsaved))
        f.write("                     configurations saved\n\n")
        f.write(" ")
        f.write("{:>10d}".format(n))
        f.write(" molecules (of all types)\n")
        f.write(" ")
        f.write("{:>10d}".format(nmol_types))
        f.write(" types of molecule\n")
        f.write(" ")
        f.write("{:>10d}".format(nsites_max))
        f.write(" is the largest number of atoms in a molecule\n")
        f.write(" ")
        f.write("{:>10d}".format(neuler))
        f.write(" Euler angles are provided\n\n")
        f.write("          ")
        if truncated == True:
            f.write('T')
        else:
            f.write('F')
        f.write(" (Box is not truncated octahedral)\n")
        f.write("           ")
        f.write(" Defining vectors are:\n")            
        for i in range(3):
            f.write("           ")
            f.write('{:11.06f}'.format(matrix[i][0]))
            f.write('{:11.06f}'.format(matrix[i][1]))
            f.write('{:11.06f}\n'.format(matrix[i][2]))
        f.write('\n')
        
        for i in range(nmol_types):
            f.write("{:>11d}".format(ni[i]))
            f.write(" molecules of type")
            f.write("{:>3d}\n".format(i+1))
            f.write("{:>11d}".format(nsites[i]))                
            f.write(" atomic sites\n")
            f.write("           ")
            f.write('{:11.06f}'.format(sites[i][0]))
            f.write('{:11.06f}'.format(sites[i][1]))
            f.write('{:11.06f}'.format(sites[i][2]))
            f.write('\n\n')
        
        for i in range(n):
            position = minv.dot(np.array(positions[i]))
            f.write('{:11.07f}'.format(position[0]))
            f.write('{:11.07f}'.format(position[1]))
            f.write('{:11.07f}\n'.format(position[2]))