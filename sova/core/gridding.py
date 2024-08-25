
import math
import numpy as np
import sys

def metric(vectors):
    _metric = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                _metric[i][j] = _metric[i][j]+vectors[k][i]*vectors[k][j]
    return _metric

def volume(vectors, truncated=False):
    triprod = vectors[0][0]*vectors[1][1]*vectors[2][2] \
            + vectors[1][0]*vectors[2][1]*vectors[0][2] \
            + vectors[2][0]*vectors[0][1]*vectors[1][2] \
            - vectors[2][0]*vectors[1][1]*vectors[0][2] \
            - vectors[1][0]*vectors[0][1]*vectors[2][2] \
            - vectors[0][0]*vectors[2][1]*vectors[1][2]
    _volume = 8.0*abs(triprod)
    if truncated:
        _volume /= 2.0

    return _volume

def d(vectors, truncated=False):
    triprod = vectors[0][0]*vectors[1][1]*vectors[2][2] \
            + vectors[1][0]*vectors[2][1]*vectors[0][2] \
            + vectors[2][0]*vectors[0][1]*vectors[1][2] \
            - vectors[2][0]*vectors[1][1]*vectors[0][2] \
            - vectors[1][0]*vectors[0][1]*vectors[2][2] \
            - vectors[0][0]*vectors[2][1]*vectors[1][2]

    axb1=vectors[1][0]*vectors[2][1]-vectors[2][0]*vectors[1][1]
    axb2=vectors[2][0]*vectors[0][1]-vectors[0][0]*vectors[2][1]
    axb3=vectors[0][0]*vectors[1][1]-vectors[1][0]*vectors[0][1]
    bxc1=vectors[1][1]*vectors[2][2]-vectors[2][1]*vectors[1][2]
    bxc2=vectors[2][1]*vectors[0][2]-vectors[0][1]*vectors[2][2]
    bxc3=vectors[0][1]*vectors[1][2]-vectors[1][1]*vectors[0][2]
    cxa1=vectors[1][2]*vectors[2][0]-vectors[2][2]*vectors[1][0]
    cxa2=vectors[2][2]*vectors[0][0]-vectors[0][2]*vectors[2][0]
    cxa3=vectors[0][2]*vectors[1][0]-vectors[1][2]*vectors[0][0]
    d1 = triprod/math.sqrt(axb1**2+axb2**2+axb3**2)
    d2 = triprod/math.sqrt(bxc1**2+bxc2**2+bxc3**2)
    d3 = triprod/math.sqrt(cxa1**2+cxa2**2+cxa3**2)
    _d = min(d1,d2,d3)
    if (truncated):        
        d1 = 1.5*triprod/math.sqrt( \
            (axb1+bxc1+cxa1)**2+(axb2+bxc2+cxa2)**2+(axb3+bxc3+cxa3)**2)
        d2 = 1.5*triprod/math.sqrt( \
            (axb1-bxc1+cxa1)**2+(axb2-bxc2+cxa2)**2+(axb3-bxc3+cxa3)**2)
        d3 = 1.5*triprod/math.sqrt( \
            (axb1+bxc1-cxa1)**2+(axb2+bxc2-cxa2)**2+(axb3+bxc3-cxa3)**2)
        d4 = 1.5*triprod/math.sqrt( \
            (axb1-bxc1-cxa1)**2+(axb2-bxc2-cxa2)**2+(axb3-bxc3-cxa3)**2)
        _d = min(_d,d1,d2,d3,d4)
    
    return _d
    
class Grid(object):
    """
    :class RMC Grid class
    """
    def __init__(self,centres,vectors):
        # atoms positions
        self.centres = centres
        self.vectors = vectors
        self.metric  = metric(self.vectors)
        self.truncated = False
        self.nmol = self.centres.shape[0]
        self.grid_size = 30
        self.grid = [[[[] for i in range(self.grid_size+1)] \
                          for j in range(self.grid_size+1)] \
                          for k in range(self.grid_size+1)]
        self.grid_co = np.zeros((self.nmol, 3), dtype=int)
        self.make_grid()
        
    def make_grid(self):
        for i in range(self.nmol):
            ix = int((self.centres[i][0]+1.0)*self.grid_size/2.0)+1
            iy = int((self.centres[i][1]+1.0)*self.grid_size/2.0)+1
            iz = int((self.centres[i][2]+1.0)*self.grid_size/2.0)+1
            ix = min(ix, self.grid_size)
            iy = min(iy, self.grid_size)
            iz = min(iz, self.grid_size)            

            self.grid[ix][iy][iz].append(i)
            self.grid_co[i][0] = ix
            self.grid_co[i][1] = iy
            self.grid_co[i][2] = iz

        dis = d(self.vectors)
        self.cell_width = 2.0*dis/self.grid_size
        self.rmax = dis
        
    def update_grid(self, imove, xold, yold, zold, xnew, ynew, znew):
        ix = int((xold+1.0)*self.grid_size/2.0)+1
        iy = int((yold+1.0)*self.grid_size/2.0)+1
        iz = int((zold+1.0)*self.grid_size/2.0)+1
        ix = min(ix, self.grid_size)
        iy = min(iy, self.grid_size)
        iz = min(iz, self.grid_size)
        
        if imove in self.grid[ix][iy][iz]:
            self.grid[ix][iy][iz].remove(imove)
        else:
            print('Not found {} in grid : Grid.update_grid'.format(imove))
            sys.exit()

        ix = int((xnew+1.0)*self.grid_size/2.0)+1
        iy = int((ynew+1.0)*self.grid_size/2.0)+1
        iz = int((znew+1.0)*self.grid_size/2.0)+1
        ix = min(ix, self.grid_size)
        iy = min(iy, self.grid_size)
        iz = min(iz, self.grid_size)
        self.grid[ix][iy][iz].append(imove)
        self.grid_co[imove,0] = ix
        self.grid_co[imove,1] = iy
        self.grid_co[imove,2] = iz
        
    def neighbours(self, ic, rmax, n1, n2):
        """
        rmax : search in max radius
        """
        #metric = utils.metric(self.vectors)
        self.inei = []
        self.coords = [[] for i in range(3)]
        self.d = []
        self.shifts = [[] for i in range(3)]
                
        #rmax = min(self.rmax, rmax)
        
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
                            x = self.centres[ig][0]-self.centres[ic][0]+3.
                            y = self.centres[ig][1]-self.centres[ic][1]+3.
                            z = self.centres[ig][2]-self.centres[ic][2]+3.
                            shiftx = int(x/2.)-1
                            shifty = int(y/2.)-1
                            shiftz = int(z/2.)-1
                            x = 2.*(x/2.-int(x/2.))-1.
                            y = 2.*(y/2.-int(y/2.))-1.
                            z = 2.*(z/2.-int(z/2.))-1.
                    
                            if (self.truncated == True and \
                                    math.fabs(x)+math.fabs(y)+math.fabs(z) > 1.5):
                                
                                x = x-self.sign(1.,x)
                                y = y-self.sign(1.,y)
                                z = z-self.sign(1.,z)
                            
                            d2 = self.metric[0][0]*x*x+self.metric[1][1]*y*y+self.metric[2][2]*z*z \
                                + 2.0*(self.metric[0][1]*x*y+self.metric[0][2]*x*z+self.metric[1][2]*y*z)

                            dis = math.sqrt(d2)
                            if(dis <= rmax):
                                if ig in self.inei: continue
                                neigh = neigh+1
                                self.inei.append(ig)
                                self.coords[0].append(x)
                                self.coords[1].append(y)
                                self.coords[2].append(z)
                                self.d.append(dis)
                                self.shifts[0].append(shiftx)
                                self.shifts[1].append(shifty)
                                self.shifts[2].append(shiftz)
                
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
            sx = self.shifts[0][imin]
            sy = self.shifts[1][imin]
            sz = self.shifts[2][imin]
            
            dd = self.d[imin]
            self.inei[imin] = self.inei[i]
            self.coords[0][imin] = self.coords[0][i]
            self.coords[1][imin] = self.coords[1][i]
            self.coords[2][imin] = self.coords[2][i]
            self.shifts[0][imin] = self.shifts[0][i]
            self.shifts[1][imin] = self.shifts[1][i]
            self.shifts[2][imin] = self.shifts[2][i]

            self.d[imin] = self.d[i]
            self.inei[i] = ini
            self.coords[0][i] = xd
            self.coords[1][i] = yd
            self.coords[2][i] = zd
            self.d[i] = dd
            self.shifts[0][i] = sx
            self.shifts[1][i] = sy
            self.shifts[2][i] = sz

        # convert numpy array
        #self.shifts = np.array(self.shifts).T
        #self.coords = np.array(self.coords).T        

    @staticmethod
    def sign(a, b):
        if(b >= 0.0):
            return math.fabs(a)
        else:
            return -math.fabs(a)
