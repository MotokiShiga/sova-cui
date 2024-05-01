import os
import sys
import math

from .rmc_configuration import RmcConfiguration

import numpy as np
import matplotlib.pyplot as plt

class Grid(object):
    """
    :class RMC Grid class
    """
    def __init__(self, cfg):
        self.cfg = cfg
        self.grid_size = 30
        self.grid = [[[[] for i in range(self.grid_size+1)] \
                          for j in range(self.grid_size+1)] \
                          for k in range(self.grid_size+1)]
        
        self.grid_co = np.zeros((self.cfg.nmol, 3), dtype=int)

    def make_grid(self):
 
        # atoms positions
        centres = self.cfg.atoms.positions

        for i in range(self.cfg.nmol):
            ix = int((centres[i][0]+1.0)*self.grid_size/2.0)+1
            iy = int((centres[i][1]+1.0)*self.grid_size/2.0)+1
            iz = int((centres[i][2]+1.0)*self.grid_size/2.0)+1
            ix = min(ix, self.grid_size)
            iy = min(iy, self.grid_size)
            iz = min(iz, self.grid_size)

            self.grid[ix][iy][iz].append(i)
            self.grid_co[i][0] = ix
            self.grid_co[i][1] = iy
            self.grid_co[i][2] = iz

        self.cell_width = 2.0*self.cfg.d/self.grid_size

    def neighbours(self, ic, rmax, n1, n2):
        """
        rmax : search in max radius
        """
        centres = self.cfg.atoms.positions
        metric = self.cfg.metric
        self.inei = []
        self.coords = [[] for i in range(3)]
        self.d = []

        neigh = 0
        ng = int(rmax/self.cell_width)+1
        gridx = self.grid_co[ic][0]
        gridy = self.grid_co[ic][1]
        gridz = self.grid_co[ic][2]
        for ix in range(gridx-ng, gridx+ng+1):
            iix = ix
            if (iix <= 0): iix = iix+self.grid_size
            if (iix > self.grid_size): iix = iix-self.grid_size

            for iy in range(gridy-ng, gridy+ng+1):
                iiy = iy
                if (iiy <= 0): iiy = iiy+self.grid_size
                if (iiy > self.grid_size): iiy = iiy-self.grid_size

                for iz in range(gridz-ng, gridz+ng+1):
                    iiz = iz
                    if (iiz <= 0): iiz = iiz+self.grid_size
                    if (iiz > self.grid_size): iiz = iiz-self.grid_size

                    for ino in range(len(self.grid[iix][iiy][iiz])):
                        ig = self.grid[iix][iiy][iiz][ino]
                        if(ig >= n1 and ig <= n2 and ig != ic):
                            x = centres[ig][0]-centres[ic][0]+3.
                            y = centres[ig][1]-centres[ic][1]+3.
                            z = centres[ig][2]-centres[ic][2]+3.
                            x = 2.*(x/2.-int(x/2.))-1.
                            y = 2.*(y/2.-int(y/2.))-1.
                            z = 2.*(z/2.-int(z/2.))-1.
                    
                            if (self.cfg.truncated == True and \
                                    math.fabs(x)+math.fabs(y)+math.fabs(z) > 1.5):
                                
                                x = x-sign(1.,x)
                                y = y-sign(1.,y)
                                z = z-sign(1.,z)
                            
                            dd = metric[0][0]*x*x+metric[1][1]*y*y+metric[2][2]*z*z \
                                + 2.0*(metric[0][1]*x*y+metric[0][2]*x*z+metric[1][2]*y*z)

                            dd = math.sqrt(dd)
                                
                            if(dd < rmax):
                                neigh = neigh+1
                                self.inei.append(ig)
                                self.coords[0].append(x)
                                self.coords[1].append(y)
                                self.coords[2].append(z)
                                self.d.append(dd)
        
        # Now sort into order
        for i in range(neigh):
            imin = i
            for j in range(i+1, neigh):
                if (self.d[j] < self.d[imin]): 
                    imin = j
                
            ini = self.inei[imin]
            xd = self.coords[0][imin]
            yd = self.coords[1][imin]
            zd = self.coords[2][imin]
            
            dd = self.d[imin]
            self.inei[imin] = self.inei[i]
            self.coords[0][imin] = self.coords[0][i]
            self.coords[1][imin] = self.coords[1][i]
            self.coords[2][imin] = self.coords[2][i]

            self.d[imin] = self.d[i]
            self.inei[i] = ini
            self.coords[0][i] = xd
            self.coords[1][i] = yd
            self.coords[2][i] = zd
            self.d[i] = dd

                            
def sign(a, b):
    if(b >= 0.0):
        return math.fabs(a)
    else:
        return -math.fabs(a)

if __name__ == '__main__':
    cfg = RmcConfiguration("../sio.cfg")
    cfg.read()

    grid = Grid(cfg)
    grid.make_grid()

    rmax = 5.0
    ic = 0
    n1 = 1000
    n2 = cfg.nmol
    grid.neighbours(ic, rmax, n1, n2)

    #x = np.random.normal(50, 10, 1000)
    plt.hist(grid.d)
    plt.show()
    #plt.hist(grid.d, bins=16)

    #print grid.inei
    #print grid.d    
    #print grid.coords[0]
