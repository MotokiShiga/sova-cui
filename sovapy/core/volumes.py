from math import ceil, sin, cos, pi, sqrt, acos
cot = lambda alpha: cos(alpha) / sin(alpha)
import itertools
import numpy as np
import numpy.linalg as la

try:
    import numexpr as ne
    NUMEXPR = True
except ImportError:
    NUMEXPR = False

VOLUME_TYPES = ["HEX","CUB","ORT","MON","TET","TRI","RHO"]

class CellInfo(object):
    """Abstract class for cell information

    Cell information that contain crystal system, space group, etc.

    Attributes
    ----------
    crystal_system : string
        The name of the crystal system
        (cubic, orthorhombic, monoclinic, tetragonal, triclinic, rhombohedral)
    space_group_number : int
        Space group number
    Hall_symbol : string
        Hall notation of the space group
    Hermann_Mauguin_symbol : string
        Hermann_Mauguin of the space group
    origin : numpy.ndarray
        The origin of the coordinate system
    periodic : bool
        If the system is periodic?
    periodic_boundary : list of float
        The minimum and maximum of the coordinates
    """
    def __init__(self):
        self.crystal_system = None
        self.space_group_number = None
        self.Hall_symbol = None
        self.Hermann_Mauguin_symbol = None
        self.origin = np.zeros(3)
        self.periodic = True
        self.periodic_boundary = [-0.5,0.5]

class TriclinicVolume(CellInfo):
    """ Class for triclinic crystal system

    Triclinic crystal system (sub-class of CellInfo)
    The system can be described lattice constants (six variables)

    Attributes
    ----------
    crystal_system : string
        The name of the crystal system
        (cubic, orthorhombic, monoclinic, tetragonal, triclinic, rhombohedral)
    space_group_number : int
        Space group number
    Hall_symbol : string
        Hall notation of the space group
    Hermann_Mauguin_symbol : string
        Hermann_Mauguin of the space group
    origin : numpy.ndarray
        The origin of the coordinate system
    periodic : bool
        If the system is periodic?
    periodic_boundary : list of float
        The minimum and maximum of the coordinates
    a : float
        The length of the 1st lattice vector
    b : float
        The length of the 2nd lattice vector
    c : float
        The length of the 3rd lattice vector
    alpha : float
        The angle between the 1st and 2nd lattice vector, given in radian
    beta : float
        The angle between the 1st and 3rd lattice vector, given in radian
    gamma : float
        The angle between the 2nd and 3rd lattice vector, given in radian
    volume : float
        The volume of the cell (simulation box)
    V : float
        The value of 'volume' devided by (a b c)
    M : numpy.ndarray
        The cartesian-to-fractional transformation matrix
    Minv : numpy.ndarray
        The fractional-to-cartesian transformation matrix
    side_lengths : list of float
        The side lengths of an axis-aligned bounding box
    translation_vectors : type
        The lattice system translation vectors (`right`, `up`, `forward`)
    edges : list of float
        List of edges of the cell (bounding box)
    _vectors : numpy.ndarray
        (Vectors computed by Minv?)
    truncated : bool
        (Related with the coordinate definition?)
    """
    def __init__(self, *args):
        super().__init__()
        
        if len(args) == 6:
            # Lattice constants
            # Lattice vectors
            a, b, c = args[0], args[1], args[2]
            # Lattice angles
            alpha, beta, gamma = args[3], args[4], args[5]
            # Save lattice constants
            self.a, self.b, self.c = a, b, c
            self.alpha, self.beta, self.gamma = alpha, beta, gamma            
            # The volume of the cell divided by a b c
            self.V = (1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+2*cos(alpha)*cos(beta)*cos(gamma))**0.5
            # The cell volume
            self.volume = self.V*self.a*self.b*self.c
            # The cartesian-to-fractional transformation matrix
            self.M = np.matrix([[1/a, 1/a * (-cos(gamma)/sin(gamma)), 1/a * (cos(alpha)*cos(gamma)-cos(beta))/(self.V*sin(gamma))],
                                [  0, 1/b *             1/sin(gamma), 1/b * (cos(beta)*cos(gamma)-cos(alpha))/(self.V*sin(gamma))],
                                [  0,                              0, 1/c * sin(gamma)/self.V]])
            eps = 1.e-6
            self.M[np.abs(self.M) < eps] = 0.0
            # The fracional-to-cartesian transformation matrix
            self.Minv = la.inv(self.M)            
        else:
            # Lattice constants
            # Lattice vectors
            v1 = np.array([float(f) for f in args[0:3]])
            v2 = np.array([float(f) for f in args[3:6]])
            v3 = np.array([float(f) for f in args[6:]])
            self.a = la.norm(v1)
            self.b = la.norm(v2)
            self.c = la.norm(v3)
            # Lattice angles
            alpha = acos((v1.dot(v3))/(self.a * self.c))
            beta = acos((v2.dot(v3))/(self.b * self.c))
            gamma = acos((v1.dot(v2))/(self.a * self.b))
            self.alpha = alpha
            self.beta = beta
            self.gamma = gamma

            # The volume of the cell divided by a b c
            self.V = (1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+2*cos(alpha)*cos(beta)*cos(gamma))**0.5
            # The cell volume
            self.volume = self.V*self.a*self.b*self.c
            # The fracional-to-cartesian transformation matrix
            self.Minv = np.matrix(np.array([v1, v2, v3]).T)
            # The cartesian-to-fractional transformation matrix
            self.M = np.matrix(la.inv(self.Minv))

        # Optional variables
        self._vectors = None
        self.truncated = False
        
        # The lattice system translation vectors (`right`, `up`, `forward`)
        self.translation_vectors = [tuple(self.Minv.T[i].tolist()[0]) for i in range(3)]
        min_point = [float('inf')]*3
        max_point = [float('-inf')]*3
        for i, j, k in itertools.product((-0.5, 0, 0.5), repeat=3):
            point = self.Minv*np.matrix((i, j, k)).T
            point = point.T.tolist()[0]
            for l in range(3):
                min_point[l] = min(min_point[l], point[l])
                max_point[l] = max(max_point[l], point[l])

        # The side lengths of an axis-aligned bounding box
        self.side_lengths = [d-c for c, d in zip(min_point, max_point)]

        # Compute list of edges of the cell
        edges = [((-0.5, -0.5, -0.5), (-0.5, -0.5, 0.5)), \
                 ((-0.5, -0.5, -0.5), (-0.5, 0.5, -0.5)), \
                 ((-0.5, -0.5, -0.5), (0.5, -0.5, -0.5)), \
                 ((-0.5, -0.5, 0.5), (-0.5, 0.5, 0.5)), \
                 ((-0.5, -0.5, 0.5), (0.5, -0.5, 0.5)), \
                 ((-0.5, 0.5, -0.5), (-0.5, 0.5, 0.5)), \
                 ((-0.5, 0.5, -0.5), (0.5, 0.5, -0.5)), \
                 ((-0.5, 0.5, 0.5), (0.5, 0.5, 0.5)), \
                 ((0.5, -0.5, -0.5), (0.5, -0.5, 0.5)), \
                 ((0.5, -0.5, -0.5), (0.5, 0.5, -0.5)), \
                 ((0.5, -0.5, 0.5), (0.5, 0.5, 0.5)), \
                 ((0.5, 0.5, -0.5), (0.5, 0.5, 0.5))]
        new_edges = []
        for edge in edges:
            point1, point2 = edge
            point1 = self.Minv*np.matrix(point1).T
            point1 = point1.T.tolist()[0]
            point2 = self.Minv*np.matrix(point2).T
            point2 = point2.T.tolist()[0]
            new_edges.append((point1, point2))
        #: A list of edges as point pairs
        self.edges = new_edges

    def set_edges_from_periodic(self, periodic_boundary=[-0.5, 0.5]):
        """Compute and set edges of the cell

        Compute and set edges of the cell (boundaray box)
        using the minimum and maximum of the coordinates.

        Parameters
        ----------
        periodic_boundary : list, optional
            The minimum and maximum of the coordinates, by default [-0.5, 0.5]
        """
        if periodic_boundary == [-0.5, 0.5]:
            edges = [((-0.5, -0.5, -0.5), (-0.5, -0.5, 0.5)), \
                     ((-0.5, -0.5, -0.5), (-0.5, 0.5, -0.5)), \
                     ((-0.5, -0.5, -0.5), (0.5, -0.5, -0.5)), \
                     ((-0.5, -0.5, 0.5), (-0.5, 0.5, 0.5)), \
                     ((-0.5, -0.5, 0.5), (0.5, -0.5, 0.5)), \
                     ((-0.5, 0.5, -0.5), (-0.5, 0.5, 0.5)), \
                     ((-0.5, 0.5, -0.5), (0.5, 0.5, -0.5)), \
                     ((-0.5, 0.5, 0.5), (0.5, 0.5, 0.5)), \
                     ((0.5, -0.5, -0.5), (0.5, -0.5, 0.5)), \
                     ((0.5, -0.5, -0.5), (0.5, 0.5, -0.5)), \
                     ((0.5, -0.5, 0.5), (0.5, 0.5, 0.5)), \
                     ((0.5, 0.5, -0.5), (0.5, 0.5, 0.5))]
        else:
            edges = [((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), \
                     ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)), \
                     ((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)), \
                     ((0.0, 0.0, 1.0), (0.0, 1.0, 1.0)), \
                     ((0.0, 0.0, 1.0), (1.0, 0.0, 1.0)), \
                     ((0.0, 1.0, 0.0), (0.0, 1.0, 1.0)), \
                     ((0.0, 1.0, 0.0), (1.0, 1.0, 0.0)), \
                     ((0.0, 1.0, 1.0), (1.0, 1.0, 1.0)), \
                     ((1.0, 0.0, 0.0), (1.0, 0.0, 1.0)), \
                     ((1.0, 0.0, 0.0), (1.0, 1.0, 0.0)), \
                     ((1.0, 0.0, 1.0), (1.0, 1.0, 1.0)), \
                     ((1.0, 1.0, 0.0), (1.0, 1.0, 1.0))]            
        new_edges = []
        for edge in edges:
            point1, point2 = edge
            point1 = self.Minv*np.matrix(point1).T
            point1 = point1.T.tolist()[0]
            point2 = self.Minv*np.matrix(point2).T
            point2 = point2.T.tolist()[0]
            new_edges.append((point1, point2))
        # List of edges as point pairs
        self.edges = new_edges
    
    @property
    def volume_from_vectors(self):
        """Compute the cell volume from lattice vectors

        Returns
        -------
        volume : float
            The cell volume
        """
        triprod = self.vectors[0][0]*self.vectors[1][1]*self.vectors[2][2] \
                + self.vectors[1][0]*self.vectors[2][1]*self.vectors[0][2] \
                + self.vectors[2][0]*self.vectors[0][1]*self.vectors[1][2] \
                - self.vectors[2][0]*self.vectors[1][1]*self.vectors[0][2] \
                - self.vectors[1][0]*self.vectors[0][1]*self.vectors[2][2] \
                - self.vectors[0][0]*self.vectors[2][1]*self.vectors[1][2]
        _volume = 8.0*abs(triprod)
        if self.truncated == True:
            _volume = _volume/2.0

        return _volume

    @property
    def vectors(self):
        """Getter of the lattice vectors

        Returns
        -------
        _vectors : numpy.array
            The lattice vectors
        """
        if self._vectors is None:
            self._vectors = 0.5*np.array(self.Minv)
        return self._vectors

    def is_inside(self, point, periodic_boundary=[-0.5,0.5]):
        """Evaluate if the point is inside of the cell

        Evaluate if point is inside of the cell.

        Parameters
        ----------
        point : numpy.ndarray
            3D coordinate (cartesian)
        periodic_boundary : list of float, optional
            The minimum and maximum of the coordinates, by default [-0.5,0.5]

        Returns
        -------
        flag : bool
            if point is inside of the cell
        """
        if isinstance(point, np.ndarray) and len(point.shape) > 1:
            fp = np.asarray((self.M * point.T).T)
            #return np.all(np.abs(fp) < 0.5, axis=1)
            return np.all(np.abs(fp) < periodic_boundary[1], axis=1)
        else:
            fractional_point = self.M*np.matrix(point).T
            return all((periodic_boundary[0] < float(c) < periodic_boundary[1] for c in fractional_point))

    def get_equivalent_point(self, point, periodic_boundary=0.5):
        """Enumerate equivalent point

        Identify an equivalent point of a given point in inside the volume.

        Parameters
        ----------
        point : numpy.ndarray
            3D Coordinate (cartesian)
        periodic_boundary : float, optional
            The maximum coordinate, by default 0.5

        Returns
        -------
        new_point : tuple of float
            Equivalent points in the cell
        """
        eps = 1.e-8
        #fractional_point = self.M*np.matrix(point).T
        fractional_point = np.array(self.M).dot(point)
        fractional_point[np.abs(fractional_point) < eps] = 0.0
        
        #fractional_point = fractional_point.T.tolist()[0]
        for i in range(3): 
            if periodic_boundary == 1.0:
                if fractional_point[i] < 0.0:
                    fractional_point[i] -= ceil(fractional_point[i]-1.0)
                #elif fractional_point[i] >= 1.0:
                #    fractional_point[i] -= floor(fractional_point[i])
            elif periodic_boundary == 0.5:            
                #fractional_point[i] -= ceil(fractional_point[i]-0.5)
                fractional_point[i] -= ceil(fractional_point[i]-periodic_boundary)
        new_point = self.Minv*np.matrix(fractional_point).T
        new_point = tuple(new_point.T.tolist()[0])
        return new_point

    def get_distance(self, p1, p2):
        """Compute the vector between two points

        Compute the shortest distance vector between two points.

        Parameters
        ----------
        p1 : numpy.ndarray
            3D coordinate (cartesian)
        p2 : numpy.ndarray
            3D coordinate (cartesian)

        Returns
        -------
        numpy.ndarray
            3D coordinate (cartesian)
        """
        tp1 = (self.M * np.matrix(p1, copy=False).T).T
        tp2 = (self.M * np.matrix(p2, copy=False).T).T
        td = tp2 - tp1
        td -= np.ceil(td - 0.5)
        return np.array((self.Minv * td.T).T, copy=False)

    def __repr__(self):
        return "TRICLINIC a=%f b=%f c=%f alpha=%f beta=%f gamma=%f" % (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)

    def __str__(self):
        return "TRI " + " ".join([str(x) for v in self.translation_vectors for x in v])

class MonoclinicVolume(TriclinicVolume):
    """ Class for moclinic crystal system

    Monoclinic crystal system (sub-class of TriclinicVolume),
    which a special case of a triclinic volume with ``alpha=gamma=pi/2``.
    This class is initialized using four lattice constants (a, b, c, and beta)

    Attributes
    ----------
    crystal_system : string
        The name of the crystal system
        (cubic, orthorhombic, monoclinic, tetragonal, triclinic, rhombohedral)
    space_group_number : int
        Space group number
    Hall_symbol : string
        Hall notation of the space group
    Hermann_Mauguin_symbol : string
        Hermann_Mauguin of the space group
    origin : numpy.ndarray
        The origin of the coordinate system
    periodic : bool
        If the system is periodic?
    periodic_boundary : list of float
        The minimum and maximum of the coordinates
    a : float
        The length of the 1st lattice vector
    b : float
        The length of the 2nd lattice vector
    c : float
        The length of the 3rd lattice vector
    alpha : float
        The angle between the 1st and 2nd lattice vector, given in radian
    beta : float
        The angle between the 1st and 3rd lattice vector, given in radian
    gamma : float
        The angle between the 2nd and 3rd lattice vector, given in radian
    volume : float
        The volume of the cell (simulation box)
    V : float
        The value of 'volume' devided by (a b c)
    M : numpy.ndarray
        The cartesian-to-fractional transformation matrix
    Minv : numpy.ndarray
        The fractional-to-cartesian transformation matrix
    side_lengths : list of float
        The side lengths of an axis-aligned bounding box
    translation_vectors : type
        The lattice system translation vectors (`right`, `up`, `forward`)
    edges : list of float
        List of edges of the cell (bounding box)
    _vectors : numpy.ndarray
        (Vectors computed by Minv?)
    truncated : bool
        (Related with the coordinate definition?)
    """
    def __init__(self, *args):
        if len(args) == 4:
            a = args[0]
            b = args[1]
            c = args[2]
            beta = args[3]
            alpha = pi/2
            gamma = pi/2
            TriclinicVolume.__init__(self, a, b, c, alpha, beta, gamma)
        else:
            v1 = np.array([float(f) for f in args[0:3]])
            v2 = np.array([float(f) for f in args[3:6]])
            v3 = np.array([float(f) for f in args[6:]])
            TriclinicVolume.__init__(self, v1, v2, v3)

    def __repr__(self):
        return "MONOCLINIC a=%f b=%f c=%f beta=%f" % (self.a, self.b, self.c, self.beta)

    def __str__(self):
        return "MON %f %f %f %f" % (self.a, self.b, self.c, self.beta)

class OrthorhombicVolume(TriclinicVolume):
    """ Class for orthorhombic crystal system

    Orthorhombic crystal system (sub-class of TriclinicVolume),
    which a special case of a triclinic volume with ``alpha=beta=gamma=pi/2``.
    This class is initialized using three lattice constants (a, b, and c)

    Attributes
    ----------
    crystal_system : string
        The name of the crystal system
        (cubic, orthorhombic, monoclinic, tetragonal, triclinic, rhombohedral)
    space_group_number : int
        Space group number
    Hall_symbol : string
        Hall notation of the space group
    Hermann_Mauguin_symbol : string
        Hermann_Mauguin of the space group
    origin : numpy.ndarray
        The origin of the coordinate system
    periodic : bool
        If the system is periodic?
    periodic_boundary : list of float
        The minimum and maximum of the coordinates
    a : float
        The length of the 1st lattice vector
    b : float
        The length of the 2nd lattice vector
    c : float
        The length of the 3rd lattice vector
    alpha : float
        The angle between the 1st and 2nd lattice vector, given in radian
    beta : float
        The angle between the 1st and 3rd lattice vector, given in radian
    gamma : float
        The angle between the 2nd and 3rd lattice vector, given in radian
    volume : float
        The volume of the cell (simulation box)
    V : float
        The value of 'volume' devided by (a b c)
    M : numpy.ndarray
        The cartesian-to-fractional transformation matrix
    Minv : numpy.ndarray
        The fractional-to-cartesian transformation matrix
    side_lengths : list of float
        The side lengths of an axis-aligned bounding box
    translation_vectors : type
        The lattice system translation vectors (`right`, `up`, `forward`)
    edges : list of float
        List of edges of the cell (bounding box)
    _vectors : numpy.ndarray
        (Vectors computed by Minv?)
    truncated : bool
        (Related with the coordinate definition?)
    """
    def __init__(self, *args):
        if len(args) == 3:
            a = args[0]
            b = args[1]
            c = args[2]
            alpha = pi/2
            beta = pi/2
            gamma = pi/2
            TriclinicVolume.__init__(self, a, b, c, alpha, beta, gamma)
        else:
            # cell vectors
            v1 = np.array([float(f) for f in args[0:3]])
            v2 = np.array([float(f) for f in args[3:6]])
            v3 = np.array([float(f) for f in args[6:]])
            TriclinicVolume.__init__(self, v1, v2, v3)

    def __repr__(self):
        return "ORTHORHOMBIC a=%f b=%f c=%f" % (self.a, self.b, self.c)

    def __str__(self):
        return "ORT %f %f %f" % (self.a, self.b, self.c)

class TetragonalVolume(TriclinicVolume):
    """ Class for tetragonal crystal system

    Tetragonal crystal system (sub-class of TriclinicVolume),
    which a special case of a triclinic volume with ``a=b`` and ``alpha=beta=gamma=pi/2``.
    This class is initialized using two lattice constants (a and c)

    Attributes
    ----------
    crystal_system : string
        The name of the crystal system
        (cubic, orthorhombic, monoclinic, tetragonal, triclinic, rhombohedral)
    space_group_number : int
        Space group number
    Hall_symbol : string
        Hall notation of the space group
    Hermann_Mauguin_symbol : string
        Hermann_Mauguin of the space group
    origin : numpy.ndarray
        The origin of the coordinate system
    periodic : bool
        If the system is periodic?
    periodic_boundary : list of float
        The minimum and maximum of the coordinates
    a : float
        The length of the 1st lattice vector
    b : float
        The length of the 2nd lattice vector
    c : float
        The length of the 3rd lattice vector
    alpha : float
        The angle between the 1st and 2nd lattice vector, given in radian
    beta : float
        The angle between the 1st and 3rd lattice vector, given in radian
    gamma : float
        The angle between the 2nd and 3rd lattice vector, given in radian
    volume : float
        The volume of the cell (simulation box)
    V : float
        The value of 'volume' devided by (a b c)
    M : numpy.ndarray
        The cartesian-to-fractional transformation matrix
    Minv : numpy.ndarray
        The fractional-to-cartesian transformation matrix
    side_lengths : list of float
        The side lengths of an axis-aligned bounding box
    translation_vectors : type
        The lattice system translation vectors (`right`, `up`, `forward`)
    edges : list of float
        List of edges of the cell (bounding box)
    _vectors : numpy.ndarray
        (Vectors computed by Minv?)
    truncated : bool
        (Related with the coordinate definition?)
    """
    def __init__(self, *args):
        if len(args) == 2:
            a = args[0]
            c = args[1]
            b = a
            alpha = pi/2
            beta = pi/2
            gamma = pi/2
            TriclinicVolume.__init__(self, a, b, c, alpha, beta, gamma)
        else:
            # cell vectors
            v1 = np.array([float(f) for f in args[0:3]])
            v2 = np.array([float(f) for f in args[3:6]])
            v3 = np.array([float(f) for f in args[6:]])
            TriclinicVolume.__init__(self, v1, v2, v3)

    def __repr__(self):
        return "TETRAGONAL a=%f c=%f" % (self.a, self.c)

    def __str__(self):
        return "TET %f %f" % (self.a, self.c)

class RhombohedralVolume(TriclinicVolume):
    """ Class for rhombohedral crystal system

    Rhombohedral crystal system (sub-class of TriclinicVolume),
    which a special case of a triclinic volume with ``a=b=c`` and ``alpha=beta=gamma``.
    This class is initialized using two lattice constants (a and alpha)

    Attributes
    ----------
    crystal_system : string
        The name of the crystal system
        (cubic, orthorhombic, monoclinic, tetragonal, triclinic, rhombohedral)
    space_group_number : int
        Space group number
    Hall_symbol : string
        Hall notation of the space group
    Hermann_Mauguin_symbol : string
        Hermann_Mauguin of the space group
    origin : numpy.ndarray
        The origin of the coordinate system
    periodic : bool
        If the system is periodic?
    periodic_boundary : list of float
        The minimum and maximum of the coordinates
    a : float
        The length of the 1st lattice vector
    b : float
        The length of the 2nd lattice vector
    c : float
        The length of the 3rd lattice vector
    alpha : float
        The angle between the 1st and 2nd lattice vector, given in radian
    beta : float
        The angle between the 1st and 3rd lattice vector, given in radian
    gamma : float
        The angle between the 2nd and 3rd lattice vector, given in radian
    volume : float
        The volume of the cell (simulation box)
    V : float
        The value of 'volume' devided by (a b c)
    M : numpy.ndarray
        The cartesian-to-fractional transformation matrix
    Minv : numpy.ndarray
        The fractional-to-cartesian transformation matrix
    side_lengths : list of float
        The side lengths of an axis-aligned bounding box
    translation_vectors : type
        The lattice system translation vectors (`right`, `up`, `forward`)
    edges : list of float
        List of edges of the cell (bounding box)
    _vectors : numpy.ndarray
        (Vectors computed by Minv?)
    truncated : bool
        (Related with the coordinate definition?)
    """
    def __init__(self, *args):
        if len(args) == 2:
            a = args[0]
            alpha = args[1]
            b = a
            c = a
            beta = alpha
            gamma = alpha
            TriclinicVolume.__init__(self, a, b, c, alpha, beta, gamma)
        else:
            # cell vectors
            v1 = np.array([float(f) for f in args[0:3]])
            v2 = np.array([float(f) for f in args[3:6]])
            v3 = np.array([float(f) for f in args[6:]])
            TriclinicVolume.__init__(self, v1, v2, v3)

    def __repr__(self):
        return "RHOMBOHEDRAL a=%f alpha=%f" % (self.a, self.alpha)

    def __str__(self):
        return "RHO %f %f" % (self.a, self.alpha)

class CubicVolume(TriclinicVolume):
    """ Class for cubic crystal system

    Cubic crystal system (sub-class of TriclinicVolume),
    which a special case of a triclinic volume with ``a=b=c`` and ``alpha=beta=gamma=pi/2``.
    This class is initialized using one lattice constants (a)

    Attributes
    ----------
    crystal_system : string
        The name of the crystal system
        (cubic, orthorhombic, monoclinic, tetragonal, triclinic, rhombohedral)
    space_group_number : int
        Space group number
    Hall_symbol : string
        Hall notation of the space group
    Hermann_Mauguin_symbol : string
        Hermann_Mauguin of the space group
    origin : numpy.ndarray
        The origin of the coordinate system
    periodic : bool
        If the system is periodic?
    periodic_boundary : list of float
        The minimum and maximum of the coordinates
    a : float
        The length of the 1st lattice vector
    b : float
        The length of the 2nd lattice vector
    c : float
        The length of the 3rd lattice vector
    alpha : float
        The angle between the 1st and 2nd lattice vector, given in radian
    beta : float
        The angle between the 1st and 3rd lattice vector, given in radian
    gamma : float
        The angle between the 2nd and 3rd lattice vector, given in radian
    volume : float
        The volume of the cell (simulation box)
    V : float
        The value of 'volume' devided by (a b c)
    M : numpy.ndarray
        The cartesian-to-fractional transformation matrix
    Minv : numpy.ndarray
        The fractional-to-cartesian transformation matrix
    side_lengths : list of float
        The side lengths of an axis-aligned bounding box
    translation_vectors : type
        The lattice system translation vectors (`right`, `up`, `forward`)
    edges : list of float
        List of edges of the cell (bounding box)
    _vectors : numpy.ndarray
        (Vectors computed by Minv?)
    truncated : bool
        (Related with the coordinate definition?)
    """
    def __init__(self, *args):        
        if len(args) == 1:
            a = args[0]
            b = a
            c = a
            alpha = pi/2
            beta = pi/2
            gamma = pi/2            
            TriclinicVolume.__init__(self, a, b, c, alpha, beta, gamma)
        else:            
            # cell vectors
            v1 = np.array([float(f) for f in args[0:3]])
            v2 = np.array([float(f) for f in args[3:6]])
            v3 = np.array([float(f) for f in args[6:]])
            TriclinicVolume.__init__(self, v1, v2, v3)

    def __repr__(self):
        return "CUBIC a=%f" % self.a

    def __str__(self):
        return "CUB %f" % self.a

class HexagonalVolume(CellInfo):
    '''
    A hexagonal volume centered in the origin with a side length of ``a`` (for
    the 6 individual outer sides) and a height of ``c``.
    '''
    """ Class for hexagonal crystal system

    Hexagonal crystal system (sub-class of CellInfo),
    which is centered in the origin with a side length of ``a`` (for
    the 6 individual outer sides) and a height of ``c``.

    Attributes
    ----------
    crystal_system : string
        The name of the crystal system
        (cubic, orthorhombic, monoclinic, tetragonal, triclinic, rhombohedral)
    space_group_number : int
        Space group number
    Hall_symbol : string
        Hall notation of the space group
    Hermann_Mauguin_symbol : string
        Hermann_Mauguin of the space group
    origin : numpy.ndarray
        The origin of the coordinate system
    periodic : bool
        If the system is periodic?
    periodic_boundary : list of float
        The minimum and maximum of the coordinates
    a : float
        The length of the 1st lattice vector
    c : float
        The length of the 3rd lattice vector
    gamma : float
        The angle between the 2nd and 3rd lattice vector, given in radian
    volume : float
        The volume of the cell (simulation box)
    side_lengths : list of float
        The side lengths of an axis-aligned bounding box
    translation_vectors : type
        The lattice system translation vectors (`right`, `up`, `forward`)
    edges : list of float
        List of edges of the cell (bounding box)
    """
    def __init__(self, a, c):
        super().__init__()
        self.a = float(a)
        self.c = float(c)
        self.gamma = pi/3
        f = 2*self.a*sin(pi/3)

        # The lattice system translation vectors (`right`, `right-up`, `forward`)
        self.translation_vectors = [(cos(pi*i/3)*f, sin(pi*i/3)*f, 0) for i in range(2)] + [(0, 0, self.c)]

        # The side lengths of an axis-aligned bounding box
        self.side_lengths = [2*self.a*sin(pi/3), 2*self.a, self.c]

        # The cell volume, calculated as 6 times the area of an equilateral
        # triangle with side length a (3**0.5/4*a*a) times the height c.
        self.volume = 6*(3**0.5/4*self.a*self.a)*self.c

        # List of edges as point pairs
        self.edges = []
        for i in range(6):
            # edge on the lower hexagon
            p1 = (sin(pi/3*i)*self.a, cos(pi/3*i)*self.a, -0.5*self.c)
            p2 = (sin(pi/3*(i+1))*self.a, cos(pi/3*(i+1))*self.a, -0.5*self.c)
            self.edges.append((p1, p2))
            # edge on the upper hexagon
            p3 = (sin(pi/3*i)*self.a, cos(pi/3*i)*self.a, 0.5*self.c)
            p4 = (sin(pi/3*(i+1))*self.a, cos(pi/3*(i+1))*self.a, 0.5*self.c)
            self.edges.append((p3, p4))
            # edge connecting the two hexagons
            p5 = (sin(pi/3*i)*self.a, cos(pi/3*i)*self.a, -0.5*self.c)
            p6 = (sin(pi/3*i)*self.a, cos(pi/3*i)*self.a, 0.5*self.c)
            self.edges.append((p5, p6))

    def is_inside(self, point):
        """Evaluate if the point is inside of the cell

        Evaluate if point is inside of the cell.

        Parameters
        ----------
        point : numpy.ndarray
            3D coordinate (cartesian)

        Returns
        -------
        flag : bool
            if point is inside of the cell
        """
        if isinstance(point, np.ndarray) and len(point.shape) > 1:
            if NUMEXPR:
                a = self.a
                c = self.c                
                sinpi3 = sin(pi / 3)
                cotpi3 = cot(pi / 3)
                x = point[:, 0]
                y = point[:, 1]
                z = point[:, 2]
                return ne.evaluate("(abs(z) <= c / 2) & \
                                    (abs(x) <= sinpi3 * a) & \
                                    (abs(y) <= a) & ((abs(y) <= a / 2) | \
                                    (abs(y) <= a - abs(x) * cotpi3))")
            else:
                ap = np.abs(point)
                result = ap[:, 2] <= self.c / 2
                result = np.logical_and(result, ap[:, 0] <= sin(pi / 3) * self.a)
                result = np.logical_and(result, ap[:, 1] <= self.a)
                tmp = np.logical_or(ap[:, 1] <= self.a / 2,
                                    ap[:, 1] <= self.a - ap[:, 0] * cot(pi / 3))
                result = np.logical_and(result, tmp)
                return result
        else:
            x, y, z = point
            ax = abs(x)
            ay = abs(y)
            az = abs(z)
            if az > self.c/2:
                return False
            if ax > sin(pi/3)*self.a:
                return False
            if ay > self.a:
                return False
            if ay > self.a/2 and ay > self.a - ax*cot(pi/3):
                return False
            return True

    def get_equivalent_point(self, point, periodic_boundary=0.5):
        """Enumerate equivalent point

        Identify an equivalent point of a given point in inside the volume.

        Parameters
        ----------
        point : numpy.ndarray
            3D Coordinate (cartesian)
        periodic_boundary : float, optional
            The maximum coordinate, by default 0.5

        Returns
        -------
        tuple of float
            Equivalent points in the cell
        """
        equivalent_point = np.array(point)
        translation_vectors = np.array(self.translation_vectors)
        translation_vectors = np.append(translation_vectors, [translation_vectors[1]-translation_vectors[0]])
        translation_vectors.shape = (4, 3)
        translation_vector_lengths = np.array([la.norm(v) for v in translation_vectors])
        for i, translation_vector_length in enumerate(translation_vector_lengths):
            translation_vectors[i] /= translation_vector_length
        projection_matrix = np.matrix(translation_vectors)

        may_be_outside = True
        while may_be_outside:
            projected_point = (projection_matrix * np.matrix(equivalent_point).T)
            scaled_projected_point = np.array(projected_point.flat) / translation_vector_lengths
            max_index = np.argmax(np.abs(scaled_projected_point))
            if abs(scaled_projected_point[max_index]) > 0.5:
            #if abs(scaled_projected_point[max_index]) > periodic_boundary:
                may_be_outside = True
                equivalent_point -= ceil(scaled_projected_point[max_index]-0.5)*translation_vector_lengths[max_index]*translation_vectors[max_index]
                #equivalent_point -= ceil(scaled_projected_point[max_index]-periodic_boundary)*translation_vector_lengths[max_index]*translation_vectors[max_index]
            else:
                may_be_outside = False
        return tuple(equivalent_point)

    def _wrap(self, p):
        f = sqrt(3) * self.a
        M60 = np.matrix([
            [           0.5, 0.5 * sqrt(3), 0.0],
            [-0.5 * sqrt(3),           0.5, 0.0],
            [           0.0,           0.0, 1.0]])
        p[:, 2] -= np.ceil(p[:, 2] / self.c - 0.5) * self.c
        p[:, 0] -= np.ceil(p[:, 0] / f - 0.5) * f
        p = (M60 * p.T).T
        p[:, 0] -= np.ceil(p[:, 0] / f - 0.5) * f
        p = (M60 * p.T).T
        p[:, 0] -= np.ceil(p[:, 0] / f - 0.5) * f
        p = ((M60 * M60).T * p.T).T
        return p

    def get_distance(self, p1, p2):
        """Compute the vector between two points

        Compute the shortest distance vector between two points.

        Parameters
        ----------
        p1 : numpy.ndarray
            3D coordinate (cartesian)
        p2 : numpy.ndarray
            3D coordinate (cartesian)

        Returns
        -------
        numpy.ndarray
            3D coordinate (cartesian)
        """
        p1 = self._wrap(np.matrix(p1, copy=True))
        p2 = self._wrap(np.matrix(p2, copy=True))
        d = self._wrap(p2 - p1)
        return np.array(d, copy=False)

    def __repr__(self):
        return "HEXAGONAL a=%f c=%f" % (self.a, self.c)

    def __str__(self):
        return "HEX %f %f" % (self.a, self.c)

class NonVolume(CellInfo):
    """Class for non-periodic structure

    Attributes
    ----------
    periodic : bool
        If the system is periodic?
    Minv : numpy.ndarray or None
        The fractional-to-cartesian transformation matrix
    """
    def __init__(self, *args):
        super().__init__()
        # The structure is NOT periodic.
        self.periodic = False
        # The transformation matrix cannot be defined.
        self.Minv = None
        
    def __repr__(self):
        return "NON"

    def __str__(self):
        return "NON"

class Volume(object):
    """ Class to generate a class of a crystal system 
    
    This class is for generating a sub-class of a crystal system 
    from a string using class method 'fromstring'.
    """

    volumes = {
        'HEX': (HexagonalVolume, 'ff'),
        'MON': (MonoclinicVolume, 'ffff'),
        'TRI': (TriclinicVolume, 'ffffff'),
        'ORT': (OrthorhombicVolume, 'fff'),
        'TET': (TetragonalVolume, 'ff'),
        'RHO': (RhombohedralVolume, 'ff'),
        'CUB': (CubicVolume, 'f'),
        'NON': (NonVolume, '')
    }

    convert_functions = {
        'f': float,
        'i': int,
        's': str
    }

    @classmethod
    def fromstring(cls, s):
        """Generate a sub-class of a crystal system from a string
        
        Parameters
        ----------
        s : string
            String cantaining a crystal system name and lattice constants
        
        Returns
        -------
        volume : CellInfo
            Instance of a crystal system class

        Raises
        ------
        Exception
            Error if the given string is invalid 
        """
        try:
            s = s.split()
            if len(s)>0:
                t = s[0].upper()        # volume type
                if not t in VOLUME_TYPES:
                    s = []
                    t = 'NON'
                    print('Warning: NON was set as volume type beaccuse of the wrong description in the 2nd line of xyz file.')
                    raise Exception
            else:
                t = 'NON'
                print('Warning: NON was set as volume type beaccuse cell info was not found in the 2nd line of xyz file.')
                raise Exception
            cl = cls.volumes[t][0]  # volume class
            if len(s) == 10:        # cell vectors given
                param = [float(f) for f in s[1:]]
            else:
                param_list = s[1:]
                # parsing parameter
                param = [cls.convert_functions[p](param_list[i])
                         for i, p in enumerate(cls.volumes[t][1])]
            return cl(*param)
        except Exception:
            print('Not found volume type : volumes.fromstring')
            return None
