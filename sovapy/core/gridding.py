
import math
import numpy as np
import sys

def metric(vectors):
    """Metric in inner-product space

    Metric in the inner-product space
    computed by lattice vectors 

    Parameters
    ----------
    vectors : numpy.ndarray (3,3)
        Lattice vectors (in each column)

    Returns
    -------
    _metric : numpy.ndarray (3,3)
        Metric (inner-products of basis vectors)
    """
    _metric = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                _metric[i][j] = _metric[i][j]+vectors[k][i]*vectors[k][j]
    return _metric

# def volume(vectors, truncated=False):
#     triprod = vectors[0][0]*vectors[1][1]*vectors[2][2] \
#             + vectors[1][0]*vectors[2][1]*vectors[0][2] \
#             + vectors[2][0]*vectors[0][1]*vectors[1][2] \
#             - vectors[2][0]*vectors[1][1]*vectors[0][2] \
#             - vectors[1][0]*vectors[0][1]*vectors[2][2] \
#             - vectors[0][0]*vectors[2][1]*vectors[1][2]
#     _volume = 8.0*abs(triprod)
#     if truncated:
#         _volume /= 2.0

#     return _volume

# def d(vectors, truncated=False):
def shortest_cell_length(vectors):
    """Shortest length of the cell

    Compute shortest length of the orhogonal vectors
    to a base of the parallelepiped (cell).

    Parameters
    ----------
    vectors : numpy.ndarray (3,3)
        Three lattice vectors (Lattice constants)
    truncated : bool, optional
        truncate?, by default False

    Returns
    -------
    h_min : float
        The shortest length of the cell 
    """
    # Compute the determinant of lattice matrix (1/8 of the cell volume) 
    triprod = vectors[0][0]*vectors[1][1]*vectors[2][2] \
            + vectors[1][0]*vectors[2][1]*vectors[0][2] \
            + vectors[2][0]*vectors[0][1]*vectors[1][2] \
            - vectors[2][0]*vectors[1][1]*vectors[0][2] \
            - vectors[1][0]*vectors[0][1]*vectors[2][2] \
            - vectors[0][0]*vectors[2][1]*vectors[1][2]

    # Outer product of vector a and b
    axb1=vectors[1][0]*vectors[2][1]-vectors[2][0]*vectors[1][1]
    axb2=vectors[2][0]*vectors[0][1]-vectors[0][0]*vectors[2][1]
    axb3=vectors[0][0]*vectors[1][1]-vectors[1][0]*vectors[0][1]
    # Outer product of vector b and c
    bxc1=vectors[1][1]*vectors[2][2]-vectors[2][1]*vectors[1][2]
    bxc2=vectors[2][1]*vectors[0][2]-vectors[0][1]*vectors[2][2]
    bxc3=vectors[0][1]*vectors[1][2]-vectors[1][1]*vectors[0][2]
    # Outer product of vector c and a
    cxa1=vectors[1][2]*vectors[2][0]-vectors[2][2]*vectors[1][0]
    cxa2=vectors[2][2]*vectors[0][0]-vectors[0][2]*vectors[2][0]
    cxa3=vectors[0][2]*vectors[1][0]-vectors[1][2]*vectors[0][0]

    # Orthogonal vectors to a base of the parallelepiped
    d1 = triprod/math.sqrt(axb1**2+axb2**2+axb3**2) # volume / (area of a x b)
    d2 = triprod/math.sqrt(bxc1**2+bxc2**2+bxc3**2) # volume / (area of b x c)
    d3 = triprod/math.sqrt(cxa1**2+cxa2**2+cxa3**2) # volume / (area of c x a)
    # The shortest length of orthogonal vector orthogonal to a base of the parallelepiped
    h_min = min(d1,d2,d3)

    # if truncated:        
    #     d1t = 1.5*triprod/math.sqrt( \
    #         (axb1+bxc1+cxa1)**2+(axb2+bxc2+cxa2)**2+(axb3+bxc3+cxa3)**2)
    #     d2t = 1.5*triprod/math.sqrt( \
    #         (axb1-bxc1+cxa1)**2+(axb2-bxc2+cxa2)**2+(axb3-bxc3+cxa3)**2)
    #     d3t = 1.5*triprod/math.sqrt( \
    #         (axb1+bxc1-cxa1)**2+(axb2+bxc2-cxa2)**2+(axb3+bxc3-cxa3)**2)
    #     d4t = 1.5*triprod/math.sqrt( \
    #         (axb1-bxc1-cxa1)**2+(axb2-bxc2-cxa2)**2+(axb3-bxc3-cxa3)**2)
    #     h_min = min(h_min,d1t,d2t,d3t,d4t)
        
    return h_min
    
class Grid(object):
    """
    :class RMC Grid class
    """
    def __init__(self, centres, vectors, grid_size=30):
        """Initialization

        Parameters
        ----------
        centres : numpy.ndarray (#atoms, 3)
            Normalized atomic positions whose range is [-1.0, +1.0]
        vectors : numpy.ndarray (3,3)
            Lattice (or unit cell) vectors
        grid_size : int, optional
            The number of grid points along each axis, by default 30
        """
        # Atom positions
        self.centres = centres
        # Lattice (or unit cell) vectors
        self.vectors = vectors
        # Metric in the inner-product of the lattice vectors
        self.metric  = metric(self.vectors)
        # self.truncated = False
        # The number of atoms
        # self.nmol = self.centres.shape[0]
        self.num_atoms = self.centres.shape[0]
        # The number of grid points along each axis
        self.grid_size = int(grid_size)
        # The shortest length of the cell (lattice)
        # dis = d(self.vectors)
        dis = shortest_cell_length(self.vectors)
        # The interval of the grid points
        self.cell_width = 2.0*dis/self.grid_size
        # Radius (Maximum distance)
        self.rmax = dis

        # Initialize grid data
        self.init_grid_data()
        
    def init_grid_data(self):
        """Initialize grid data
            Generate grid data using normalized atomic positions ([-1,+1])
            and lattice vectors
        """

        # Grid data taht contains atom ids in each grid cell
        self.grid = [[[[] for i in range(self.grid_size+1)] \
                          for j in range(self.grid_size+1)] \
                          for k in range(self.grid_size+1)]
        # Grid indices of atoms to refere from atom id
        self.grid_co = np.zeros((self.num_atoms, 3), dtype=int)

        for i in range(self.num_atoms):
            # Calculate the grid index of each atom 
            ix = int((self.centres[i][0]+1.0)*self.grid_size/2.0)+1
            iy = int((self.centres[i][1]+1.0)*self.grid_size/2.0)+1
            iz = int((self.centres[i][2]+1.0)*self.grid_size/2.0)+1
            ix = min(ix, self.grid_size)
            iy = min(iy, self.grid_size)
            iz = min(iz, self.grid_size)
            
            # Store the index in the grid cell
            self.grid[ix][iy][iz].append(i)
            # Store grid index to refere from the atom id
            self.grid_co[i][0] = ix
            self.grid_co[i][1] = iy
            self.grid_co[i][2] = iz
        
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
        """Search neighbors of the atom (index: ic)

        Search neighbors of the atom. 
        The neighbors, inter-atomic differences, distances and shifted grids are
        stored in self.inei, self.coords, self.d and self.shifts.

        Parameters
        ----------
        ic : int
            The index of the atom to seach neighbor atoms
        rmax : float
            The max radius to search
        n1 : int
            The first intex to be searched as neighbors
        n2 : int
            The last intex to be searched as neighbors
        """
        #metric = utils.metric(self.vectors)
        # The list of neighbors
        self.inei = []
        # The list of The inter-atomic difference in the cell
        self.coords = [[] for i in range(3)]
        # The list of the atomic distance in the cell
        self.d = []
        # The list of the difference of shifted grids
        self.shifts = [[] for i in range(3)]
                
        #rmax = min(self.rmax, rmax)
        
        neigh = 0        
        # The half number of grid points 
        #  to search radius rmax from center atom ic
        ng = int(rmax/self.cell_width)+1


        grid_x = self.grid_co[ic][0]
        grid_y = self.grid_co[ic][1]
        grid_z = self.grid_co[ic][2]
        
        for ix in range(grid_x-ng, grid_x+ng+1):
            # To refrect periodic boundary condition (PBC)
            iix = ix            
            if (iix <= 0): iix = iix+self.grid_size
            if (iix > self.grid_size): iix = iix-self.grid_size            
            
            for iy in range(grid_y-ng, grid_y+ng+1):
                # To refrect periodic boundary condition (PBC)
                iiy = iy
                if (iiy <= 0): iiy = iiy+self.grid_size
                if (iiy > self.grid_size): iiy = iiy-self.grid_size                

                for iz in range(grid_z-ng, grid_z+ng+1):
                    # To refrect periodic boundary condition (PBC)
                    iiz = iz
                    if (iiz <= 0): iiz = iiz+self.grid_size
                    if (iiz > self.grid_size): iiz = iiz-self.grid_size                    
                    
                    # Main loop for atomic positions in the grid cell
                    for ino in range(len(self.grid[iix][iiy][iiz])):
                        # Atom index in the cell
                        ig = self.grid[iix][iiy][iiz][ino]
                        if(ig >= n1 and ig <= n2 and ig != ic):
                            # (Difference between atomic positions) plus 3.0                        
                            x = self.centres[ig][0]-self.centres[ic][0]+3.
                            y = self.centres[ig][1]-self.centres[ic][1]+3.
                            z = self.centres[ig][2]-self.centres[ic][2]+3.
                            # The number of shifted grids  
                            shift_x = int(x/2.)-1
                            shift_y = int(y/2.)-1
                            shift_z = int(z/2.)-1
                            # The inter-atomic difference in the cell
                            x = 2.*(x/2.-int(x/2.))-1.
                            y = 2.*(y/2.-int(y/2.))-1.
                            z = 2.*(z/2.-int(z/2.))-1.
                    
                            # if (self.truncated == True and \
                            #         math.fabs(x)+math.fabs(y)+math.fabs(z) > 1.5):
                                
                            #     x = x-self.sign(1.,x)
                            #     y = y-self.sign(1.,y)
                            #     z = z-self.sign(1.,z)
                            
                            # Compute the square of the distance
                            d2 = self.metric[0][0]*x*x+self.metric[1][1]*y*y+self.metric[2][2]*z*z \
                                + 2.0*(self.metric[0][1]*x*y+self.metric[0][2]*x*z+self.metric[1][2]*y*z)
                            dis = math.sqrt(d2)

                            # Add the atom as a neighbor if the distance is smaller than rmax 
                            if(dis <= rmax):
                                # Skip the following procedure if the atom is in the current list
                                if ig in self.inei: continue

                                # Add data if the atom is not in the current list
                                neigh = neigh+1 # Count the number of neighbors
                                self.inei.append(ig)     #Add the atom index
                                self.coords[0].append(x) #Add the inter-atomic difference in the cell
                                self.coords[1].append(y)
                                self.coords[2].append(z)
                                self.d.append(dis)             #Add the atomic distance
                                self.shifts[0].append(shift_x) #Add the difference of shifted grids
                                self.shifts[1].append(shift_y)
                                self.shifts[2].append(shift_z)
                
        # Sort indexes according atomic distance (by the bubble sort algorithm)
        for i in range(neigh):
            # Find the index whose distance is minimum
            imin = i
            for j in range(i+1, neigh):
                if (self.d[j] < self.d[imin]): 
                    imin = j
                
            # Copy the data of the index whose distance is minimum
            ini = self.inei[imin]
            xd = self.coords[0][imin]
            yd = self.coords[1][imin]
            zd = self.coords[2][imin]
            sx = self.shifts[0][imin]
            sy = self.shifts[1][imin]
            sz = self.shifts[2][imin]
            
            # Save the data of atom i
            dd = self.d[imin]
            self.inei[imin] = self.inei[i]
            self.coords[0][imin] = self.coords[0][i]
            self.coords[1][imin] = self.coords[1][i]
            self.coords[2][imin] = self.coords[2][i]
            self.shifts[0][imin] = self.shifts[0][i]
            self.shifts[1][imin] = self.shifts[1][i]
            self.shifts[2][imin] = self.shifts[2][i]

            # Save the data of the index whose distance is minimum
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

    # @staticmethod
    # def sign(a, b):
    #     if(b >= 0.0):
    #         return math.fabs(a)
    #     else:
    #         return -math.fabs(a)
