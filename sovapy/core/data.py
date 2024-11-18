import os, sys, re
import collections
import numpy as np
import h5py
import dateutil.parser
from datetime import datetime
from . import volumes
from . import elements
from . import bonds
from .volumes import HexagonalVolume, NonVolume
from ..config.configuration import config
from . import gridding
from ..computation.rings import Ring
from ..computation.cavity import Cavity

try:
    from util.logger import Logger
    USE_LOGGER = True
except ImportError:
    USE_LOGGER = False

if USE_LOGGER: 
    logger = Logger("core.data") 
    logger.setstream("default", sys.stdout, Logger.WARNING)

__all__ = ["Atoms",
           "FileInfo",
           "Domains",
           "Cavities",
           "ResultInfo",
           "CalculatedFrames",
           "Results"]

# Get the current version number:
cur_dirname =  os.path.dirname(__file__)
file_version = os.path.join(cur_dirname, '..', '__init__.py')
with open(file_version) as fv:
    VERSION = re.search("__version__ = '(.*)'", fv.read()).group(1)

def writedataset(h5group, name, data, overwrite=True):
    """
    Write a dataset to a hdf5 file.
    The `overwrite` parameter controls the behaviour when the dataset
    already exists    

    Parameters
    ----------
    h5group : HDF5 object
        DESCRIPTION.
    name : string
        DESCRIPTION.
    data : TYPE
        DESCRIPTION.
    overwrite : TYPE, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    bool
        DESCRIPTION.

    """
    if name in h5group:
        if overwrite:
            del h5group[name]
        else:
            return False
    h5group[name] = data
    return True

class TimestampList(object):
    """
    A `list`-like structure with a fixed length to store :class:`datetime`
    objects.  For each frame it contains a calculation date or `None`.

    Attributes
    ----------
    timestamps : list
        time stamp object list.
    """
    def __init__(self, *args):
        """
        Creates an empty :class:`TimestampList` with a given length
        or copies values from a `list`.
        
        The constructor can be called in two ways:

        - ``TimestampList(num_frames)`` :
            create an empty :class:`TimestampList` with the length `num_frames`

        - ``TimestampList(list)`` :
            copy the values from the given list of strings

        Parameters
        ----------
        *args : list or int
            object list or number of frames.

        Returns
        -------
        None.

        """
        if isinstance(args[0], list):
            arr = args[0]
            self.timestamps = [None] * len(arr)
            for i, s in enumerate(arr):
                if len(s) > 0:
                    self.timestamps[i] = dateutil.parser.parse(s)
        else:
            num_frames = args[0]
            self.timestamps = [None] * num_frames
    
    @property
    def num_frames(self):
        return len(self.timestamps)

    def __getitem__(self, index):
        if 0 <= index < len(self.timestamps):
            return self.timestamps[index]
        else:
            return None

    def __setitem__(self, index, value):
        if not isinstance(value, datetime):
            raise ValueError("datetime required")
        self.timestamps[index] = value

    def __len__(self):
        return len(self.timestamps)

    def __iter__(self):
        return iter(self.timestamps)

    def hasdata(self):
        """
        Test if the list contains any data.

        Returns
        -------
        boolean
            If any item is not `None`.

        """
        return any(map(lambda x: x is not None, self.timestamps))

    def tostrlist(self):
        """
        Converts the :class:`datetime` objects to strings.
        `None` values are converted to empty strings.

        Returns
        -------
        list
            A list of strings.

        """
        return ["" if x is None else str(x) for x in self.timestamps]

    def prettystrings(self):
        """
        Converts the :class:`datetime` objects to human-readable strings.
        `None` values are converted to "X".

        Returns
        -------
        list
            A list of strings.

        """
        def fmt(t):
            if t is None:
                return "X"
            else:
                return t.strftime("%d.%m.%Y %H:%M:%S")
        return map(fmt, self.timestamps)

class CalculatedFrames(object):
    """
    Contains information what results are calculated for what frame.
    The information about the three kinds of results are stored in
    these attributes of the type :class:`TimestampList`:

        `domains`

        `surface_cavities`

        `center_cavities`
    """

    def __init__(self, *args):
        """
        Creates an empty object or reads it from a hdf5 file.

        The constructor can be called in two ways:

        - ``CalculatedFrames(num_frames)`` :
            create an empty `CalculatedFrames` which can store information
            about `num_frames` frames

        - ``CalculatedFrames(hdf5group)`` :
            read the data from this hdf5 group
        """
        if isinstance(args[0], h5py.Group):
            h5group = args[0]
            dom_ts = list(h5group["domains"])
            sur_ts = list(h5group["surface_cavities"])
            cen_ts = list(h5group["center_cavities"])
            num_frames = len(dom_ts)
            self.domains = TimestampList(dom_ts)
            self.surface_cavities = TimestampList(sur_ts)
            self.center_cavities = TimestampList(cen_ts)
        else:
            num_frames = args[0]
            self.domains = TimestampList(num_frames)
            self.surface_cavities = TimestampList(num_frames)
            self.center_cavities = TimestampList(num_frames)

    @property
    def num_frames(self):
        return self.domains.num_frames

    def hasdata(self):
        """
        Test if member contains any data.

        **Returns:**

            If the `hasdata` method of any of the :class:`TimestampList`
            attributes returns `True`
        """
        return (self.domains.hasdata() or
                self.surface_cavities.hasdata() or
                self.center_cavities.hasdata())

    def tohdf(self, h5group, overwrite=True):
        """
        Write the data to a hdf5 Group.

        **Parameters:**
            `h5group` :
                the hdf5 group in which the data will be written

            `overwrite` :

                specifies if existing data should be overwritten
        """
        if self.hasdata():
            writedataset(h5group, "domains",
                         self.domains.tostrlist(), overwrite)
            writedataset(h5group, "surface_cavities",
                         self.surface_cavities.tostrlist(), overwrite)
            writedataset(h5group, "center_cavities",
                         self.center_cavities.tostrlist(), overwrite)

class FileInfo(object):
    """
    Contains information about input files.

    Attributes
    ----------
    id : int
        id.
    """
    def __init__(self, num_frames=None, volumestr=None, volume_guessed=False):
        """
        Parameters
        ----------
        num_frames : int
            number of frames which are stored in the file
        volumestr : string
            representative string of the volume the atoms are in
        """
        self.num_frames = num_frames
        self.volumestr = volumestr
        self._volume = None
        self.volume_guessed = volume_guessed
        
    @property
    def volume(self):
        if self._volume is None and self.volumestr is not None:
            self._volume = volumes.Volume.fromstring(self.volumestr)
        return self._volume

class ResultInfo(FileInfo):
    """
    Information about a :class:`core.file.ResultFile`:
        
    Attributes
    ----------
    num_frames : int
        number of frames which are stored in the file.
    volumestr : string
        representative string of the volume the atoms are in.
    volume : object
        the corresponding `Volume` object itself.
    sourcefilepath : string
        path of the file which contains the original input data.
    calculatedframes : object
        dictionary which associates a resolution with a
        `CalculatedFrames` object
    """
    def __init__(self, *args):
        """
        The constructor can be called in two ways:

        - ``ResultInfo()`` :
            create an empty :class:`ResultInfo` object

        - ``ResultInfo(hdf5group)`` :
            read the data from this hdf5 group
        """
        super().__init__()
        self.sourcefilepath = None
        self.calculatedframes = dict()
        if len(args) > 0 and isinstance(args[0], h5py.Group):
            h5group = args[0]
            self.num_frames = int(h5group.attrs["num_frames"])
            self.volumestr = str(h5group.attrs["volume"])
            if "sourcefile" in h5group.attrs:
                self.sourcefilepath = str(h5group.attrs["sourcefile"])
            else:
                self.sourcefilepath = None
            for name, subgroup in h5group.iteritems():
                if not name.startswith("resolution"):
                    continue
                resolution = int(name[10:])
                self.calculatedframes[resolution] = CalculatedFrames(subgroup)

    def __getitem__(self, resolution):
        """
        Get a :class:`CalculatedFrames` object for the specified resulotion.

        **Parameters:**
            `resolution`:
                the resolution for which the information will be queried

        **Returns:**
             A :class:`CalculatedFrames` object with information about
             existing calculations.
        """
        if resolution not in self.calculatedframes:
            self.calculatedframes[resolution] = CalculatedFrames(self.num_frames)
        return self.calculatedframes[resolution]

    def __contains__(self, resolution):
        """
        Check if any information about the given resolution is available.

        **Parameters:**
            `resolution`:
                the resolution for which the information will be queried

        **Returns:**
             If a :class:`CalculatedFrames` object for this resolution exists
        """
        return resolution in self.calculatedframes

    def resolutions(self):
        return self.calculatedframes.keys()

    def tohdf(self, h5group, overwrite=True):
        """
        Write the data to a hdf5 Group.

        **Parameters:**
            `h5group` :
                the hdf5 group in which the data will be written

            `overwrite` :
                specifies if existing data should be overwritten
        """
        h5group.attrs["num_frames"] = self.num_frames
        h5group.attrs["volume"] = self.volumestr
        if self.sourcefilepath is not None:
            h5group.attrs["sourcefile"] = self.sourcefilepath
        elif "sourcefile" in h5group.attrs:
            del h5group.attrs["sourcefile"]
        for resolution, info in self.calculatedframes.iteritems():
            if info.hasdata():
                subgroup = h5group.require_group("resolution{}".format(resolution))
                info.tohdf(subgroup, overwrite)

# This class preserves text data of structurral information loaded to Atoms objects
class OriginalStructureData(object):

    def __init__(self, filepath, text_data=None):
        self.file_name = None       # Name of original data
        self.structure_text = None  # Text of original data            

        # Extract the file name from a full path
        self.file_name = os.path.basename(filepath)

        if text_data is None:
            # Extract text data of structure information
            if os.path.exists(filepath):
                with open(filepath,'r') as f:
                    text = f.readlines()
                    self.structure_text = ''.join(text)
            else:
                print('{:} does NOT exist!'.format(filepath))
        else:
            self.structure_text = text_data

    # Export original text data to a text file.
    def export_txt(self, output_path=None):
        
        if output_path is None:
            output_path = self.file_name

        # Check the extension. Change it if necessary.
        ori_ext = os.path.splitext(self.file_name)[-1]
        prefix_name, out_ext = os.path.splitext(output_path)
        if ori_ext != out_ext:
            print('Extention is changed to "{:}" acoording to the original file.'.format(ori_ext))
            output_path = prefix_name + ori_ext

        # Save text data of structure information
        if not os.path.exists(output_path):
            with open(output_path,'w') as f:
                f.write(self.structure_text)
            print('Structure data has exported to {:}.'.format(output_path))
        else:
            print('{:} already exists! Choose another file name.'.format(output_path))
        

class Atoms(object):
    """
    This class represents a list of atoms and their properties.

    Attributes
    ----------
    number : int
        number of atoms.
    positions : list
        list of atoms positions
    volume : float
        the volume that contains the atoms
    radii : float
        cavity cutoff radius
    sorted_positions : list
        positions` sorted from largest to smallest radius
    sorted_radii` : list
        unique `radii` sorted from largest to smallest
    radii_as_indices : list
        indices to associate an atom in `sorted_positions` with an element
        of `sorted_radii
    """
    def __init__(self, *args):
        """
        The constructor can be called in three ways

        Parameters
        ----------
        *args : list
            various list patterns.

        - ``Atoms(positions, radii, elements, volume)`` :
            create the object using the given data

        - ``Atoms(molecule, volume)`` :
            read the data from this :class:`pybel.Molecule` object

        - ``Atoms(hdf5group)`` :
            read the data from this hdf5 group

        Returns
        -------
        None.

        """
        
        # OriginalStructureData object 
        # to preserve text data of structurral information
        self.original_file_data = None

        if isinstance(args[0], h5py.Group):
            h5group = args[0]
            positions = h5group["positions"]
            radii = h5group["radii"]
            if "elements" in h5group:
                elements = h5group["elements"]
            else:
                if USE_LOGGER: 
                    logger.warn("Dataset 'elements' not found. Using 'atom' as default value")
                elements = np.empty(len(radii), dtype="|U4")
                elements[:] = "atom"
            if "volume" in h5group.attrs:
                volume = h5group.attrs["volume"]
            else:
                volume = h5group.parent.attrs["volume"]
            volume = volumes.Volume.fromstring(str(volume))
            if "bond_lengths" in h5group.attrs:
                #TODO
                bls = h5group["bond_lengths"]
                bond_lengths = dict()
                for s in bls:
                    bond_lengths[(s[0].decode(),s[1].decode())] = float(s[2])
            else:
                bond_lengths = None

        elif isinstance(args[0], Atoms):
            atoms = args[0]
            volume = atoms.volume
                        
            positions = atoms.positions
            elements = atoms.elements   
            bond_lengths = atoms.bond_lengths
            if bond_lengths is not None:
                self.set_bond_lengths(bond_lengths)
            else:
                bond_lengths = None
        else:
            # In this case, atom positions may be outside of the volume
                        
            positions = args[0]
            radii = args[1]
            elements = args[2]
            volume = args[3]
            
            is_norm = False
            if len(args) > 4:
                is_norm = args[4]
                                
            if isinstance(volume, str):
                volume = volumes.Volume.fromstring(volume)                

            if volume is not None and volume.periodic == True:
                boundary_range = volume.periodic_boundary # boundary[min,max]
                positions = list(map(lambda pos: volume.get_equivalent_point(pos,boundary_range[1]),
                                positions))
            bond_lengths = None
        self.volume = volume
        self.positions = np.array(positions, dtype=np.float64)
        self.number = self.positions.shape[0]
        self.elements = np.array(elements, dtype="|U4")
        self.numbers = dict(collections.Counter(self.elements))
        self.elements_kind = sorted([k for k, v in self.numbers.items() if v > 0])
        self.elements_count = [self.numbers[e] for e in self.elements_kind]        
        self.indices = None
        self.radii = radii
        self._covalence_radii = None
        self._covalence_radii_by_element = None        
        self._bonds = None        
        self.change_bond_lengths = False
        self.shifts = None
        self._colors = None
        self.is_norm = is_norm
        self._norm_positions = None
        self.sorted_elements = None
        self.sorted_indexes = None
        self._elementsorted_positions = None
        self._grid = None
        self.symbols = None # equal self.elements_kind
        self.ni = None
        self.frac = None
        self.pairs = None
        self.trios = None
        # initialize
        self.symbol_order()
        # self.set_bond_lengths()
        self.bond_lengths = bond_lengths
        self.set_bond_lengths(bond_lengths) 
    
    def symbol_order(self,symbols=None):
        if symbols is None:
            atomic_numbers = elements.numbers
        else:
            atomic_numbers = {str.upper(symbols[i]):(i+1) for i in range(len(symbols))}
        self.symbols, self.ni = zip(*sorted(self.numbers.items(), key=lambda x: atomic_numbers.get(str.upper(x[0]))))
        self.symbols = list(self.symbols)
        self.ni = list(self.ni)
        self.frac = self.ni/np.sum(self.ni)
        self.pairs = []
        for i in range(len(self.symbols)):
            for j in range(i,len(self.symbols)):
                self.pairs.append(self.symbols[i] + '-' + self.symbols[j])
        self.trios = []
        for i1, symbol1 in enumerate(self.symbols):
            for i2, symbol2 in enumerate(self.symbols):
                for i3 in range(i1+1):
                    trio = '{0}-{1}-{2}'.format(symbol1, symbol2, self.symbols[i3])
                    self.trios.append(trio)
    
    @property
    def covalence_radii(self):
        if self._covalence_radii is None:
            covalence_radii = np.zeros(self.number, np.float64)
            for i, element in enumerate(self.elements):
                element_number = elements.numbers[element.upper()]
                covalence_radii[i] = elements.radii[element_number]
            self._covalence_radii = covalence_radii
        return self._covalence_radii

    @property
    def covalence_radii_by_element(self):
        covalence_radii = self.covalence_radii
        if self._covalence_radii_by_element is None:
            self._covalence_radii_by_element = dict((element, covalence_radius)
                                                    for element, covalence_radius
                                                    in zip(self.elements, covalence_radii))
        return self._covalence_radii_by_element

    @property
    def norm_positions(self):
        if self.volume.Minv is None:
            return None
        self.is_norm = True        
        if self.is_norm == False:
            return self.positions
        else:
            shift = np.array(self.volume.Minv).dot(np.array([0.5,0.5,0.5]))
            if self._norm_positions is None:
                inv = np.linalg.inv(self.volume.vectors)           
                self._norm_positions = np.array(self.positions, dtype=np.float64)                
                for i in range(self.number):
                    pos = self.positions[i] - self.volume.origin - shift
                    self._norm_positions[i] = np.dot(inv, pos)
            return self._norm_positions        
    
    @property
    def maximum(self):
        return np.max(self.positions, axis=0)

    @property
    def minimum(self):
        return np.min(self.positions, axis=0)
    
    @property
    def elementsorted_positions(self):
        if self._elementsorted_positions is None:
            elems_list = self.elements.tolist()
            indexes = list(range(len(elems_list)))
            positions_list = self.norm_positions.tolist()                
            _elems, _indexes, _positions = zip(*sorted(zip(elems_list,indexes,positions_list)))
            self.sorted_elements = list(_elems)
            self.sorted_indexes = list(_indexes)
            self._elementsorted_positions = np.array(_positions)
        return self._elementsorted_positions 
    
    @property
    def bonds(self):
        if self._bonds is None or self.change_bond_lengths == True:
            self._bonds = bonds.get_bonds_with_radii(self, 1.0)
            #self._bonds = core.bonds.get_bonds_symetric_indicies(self._bonds)                        
            self.change_bond_lengths = False
        return self._bonds
    
    def bonds_by_dist(self,min_dist,max_dist):
        self.indices = np.zeros(self.number)
        for i, elem in enumerate(self.elements_kind):
            self.indices[self.elements == elem] = i
        self._bonds = bonds.get_bonds_with_element_pair(self, min_dist, max_dist, 1.0)
        return self._bonds

    def set_bond_lengths(self,bond_lengths=None):
        def set_covalence_bond():
            self.bond_lengths = {}
            for pair in self.pairs:
                elems = pair.split('-')
                dis = 0.0
                _pair = (elems[0],elems[1])
                for element in elems:
                    element_number = elements.numbers[element.upper()]
                    covalence_radii = elements.radii[element_number]
                    dis += covalence_radii
                self.bond_lengths[_pair] = dis
        
        if self.bond_lengths is None:            
            set_covalence_bond()            
        else:
            if bond_lengths is None:                
                set_covalence_bond()                
            else:
                for elems, length in bond_lengths.items():
                    pair = (elems[0], elems[1])
                    if pair in self.bond_lengths.keys():
                        self.bond_lengths[pair] = length
                    else:
                        pair = (elems[1], elems[0])
                        if pair in self.bond_lengths.keys():
                            self.bond_lengths[pair] = length
                        else:
                            print('Not found bond pair : ', elems)                        
                self.change_bond_lengths = True
                self.bonds
    
    def bond_summary(self):
        print('bond lengths : ')
        for pair, length in self.bond_lengths.items():
            print('{:<2} - {:<2} : {:.2f}'.format(pair[0], pair[1], length))
    
    @property
    def angles(self):
        bonds = self.bonds
        _angles = []
        for i, bond in enumerate(bonds):
            pair = []
            for j, elem1 in enumerate(bond):
                for k, elem2 in enumerate(bond):
                    if j < k:
                        pair.append((elem1, elem2))
            _angles.append(pair)
        return _angles

    @property
    def colors(self):
        if self._colors is None:
            colors = np.zeros((self.number, 3), np.float64)
            for i, element in enumerate(self.elements):
                element_number = elements.numbers[element.upper()]
                colors[i] = elements.colors[element_number]
            self._colors = colors/255
        return self._colors

    @property
    def radii(self):
        return self._radii

    @radii.setter
    def radii(self, values):
        if values is None:
            self._radii = np.ones((self.number), dtype=np.float64) * config.Computation.std_cutoff_radius
        elif isinstance(values, collections.abc.Mapping):
            radii = [values[elem] for elem in self.elements]
            self._radii = np.array(radii, dtype=np.float64)
        elif isinstance(values, collections.abc.Iterable):
            self._radii = np.array(values, dtype=np.float64)
        else:
            self._radii = np.ones((self.number), dtype=np.float64) * values
        indices = np.argsort(-self._radii, kind="mergesort")
        self.sorted_positions = self.positions[indices]        
        unique_radii, indices = np.unique(-self._radii, return_inverse=True)
        self.sorted_radii = -unique_radii
        self.radii_as_indices = np.sort(indices)
    
    def transforrm(self, matrix,shift):
        pass
    
    @property
    def grid(self):
        if self._grid is None:
            if isinstance(self.volume, HexagonalVolume) or isinstance(self.volume, NonVolume):
                return None
            self._grid = gridding.Grid(self.norm_positions, self.volume.vectors)
        return self._grid
    
    @property
    def rho(self):
        return self.number/gridding.volume(self.volume.vectors)

    def tohdf(self, h5group, overwrite=True):
        """
        Write the data to a hdf5 Group.

        **Parameters:**
            `h5group` :
                the hdf5 group in which the data will be written

            `overwrite` :
                specifies if existing data should be overwritten
        """
        h5group.parent.attrs["volume"] = str(self.volume)
        writedataset(h5group, "positions", self.positions, overwrite)
        writedataset(h5group, "radii", self.radii, overwrite)
        if np.any(self.elements == "atom"):
            if USE_LOGGER:
                logger.warn("Atom.elements contains default values. Not writing dataset.")
        else:
            writedataset(h5group, "elements", self.elements, overwrite)
        #TODO
        print("self.bond_lengths")
        print(self.bond_lengths)
        if self.bond_lengths is not None:
            list_bond_lengths = list()
            for a in self.bond_lengths.items():
                list_bond_lengths.append(list(a[0])+[str(a[1])])
            writedataset(h5group, "bond_lengths", list_bond_lengths, overwrite)

    def totxt(self, fmt):
        bond_file_name = fmt.format(property="bonds")
        bond_angle_file_name = fmt.format(property="bond_angles")
        bond_chain_angle_file_name = fmt.format(property="bond_dihedral_angles")

        bond_angles, bond_chain_angles = bonds.calculate_bond_angles(self, self.bonds)

        with open(bond_file_name, 'w') as outfile:
            for source_index, target_indices in enumerate(self.bonds):
                for target_index in target_indices:
                    outfile.write("{} {}\n".format(source_index+1, target_index+1))
        with open(bond_angle_file_name, 'w') as outfile:
            for bond1, bond2 in bond_angles.keys():
                if bond1[0] > bond2[1]:
                    outfile.write("{} {} {} {}\n".format(bond1[0]+1, bond1[1]+1, bond2[1]+1, bond_angles[bond1, bond2]))
        with open(bond_chain_angle_file_name, 'w') as outfile:
            for bond_chain, angle in bond_chain_angles.items():
                outfile.write("{} {} {} {}".format(*[index+1 for index in bond_chain]))
                outfile.write(" {}\n".format(angle))

class CavitiesBase(object):
    """
    Base class to store multiple surface-based objects in the 3-dimensional
    space. The :class:`Domains` and :class:`Cavities` class inherit from it.
    """

    def __init__(self, *args):
        """
        The constructor can be called in two ways:

        - ``CavitiesBase(timestamp, volumes, surface_areas, triangles)`` :
            create the object using the given data

        - ``CavitiesBase(hdf5group)`` :
            read the data from this hdf5 group
        """
        if isinstance(args[0], h5py.Group):
            def getobj_from_h5group(attr):
                if attr in h5group:
                    return h5group[attr]
                else:
                    return None

            h5group = args[0]
            timestamp = dateutil.parser.parse(h5group.attrs["timestamp"])
            number = int(h5group.attrs["number"])
            volumes = getobj_from_h5group("volumes")
            surface_areas = getobj_from_h5group("surface_areas")
            triangles = [None] * number
            for i in range(number):
                triangles[i] = getobj_from_h5group("triangles{}".format(i))
            mass_centers = getobj_from_h5group("mass_centers")
            squared_gyration_radii = getobj_from_h5group("squared_gyration_radii")
            asphericities = getobj_from_h5group("asphericities")
            acylindricities = getobj_from_h5group("acylindricities")
            anisotropies = getobj_from_h5group("anisotropies")
            characteristic_radii = getobj_from_h5group("characteristic_radii")
            cyclic_area_indices = getobj_from_h5group("cyclic_area_indices")
        else:
            (timestamp, volumes, surface_areas, triangles,
             mass_centers, squared_gyration_radii, asphericities, acylindricities, anisotropies,
             characteristic_radii, cyclic_area_indices) = args[:11]            

        if not isinstance(timestamp, datetime):
            timestamp = dateutil.parser.parse(str(timestamp))
        self.timestamp = timestamp
        self.volumes = np.array(volumes, dtype=np.float64)
        self.number = len(volumes)
        self.surface_areas = np.array(surface_areas, dtype=np.float64)
        self.triangles = [np.array(triangle, dtype=np.float64) for triangle in triangles]
        self.mass_centers = np.array(mass_centers, dtype=np.float64)
        self.squared_gyration_radii = np.array(squared_gyration_radii, dtype=np.float64)
        self.asphericities = np.array(asphericities, dtype=np.float64)
        self.acylindricities = np.array(acylindricities, dtype=np.float64)
        self.anisotropies = np.array(anisotropies, dtype=np.float64)
        self.characteristic_radii = np.array(characteristic_radii, dtype=np.float64)
        self.cyclic_area_indices = (np.array(cyclic_area_indices, dtype=np.int32)
                                    if cyclic_area_indices is not None else np.array([]))

    def tohdf(self, h5group, overwrite=True):
        """
        Write the data to a hdf5 Group.

        **Parameters:**
            `h5group` :
                the hdf5 group in which the data will be written

            `overwrite` :
                specifies if existing data should be overwritten
        """
        h5group.attrs["timestamp"] = str(self.timestamp)
        h5group.attrs["number"] = self.number
        writedataset(h5group, "volumes", self.volumes, overwrite)
        writedataset(h5group, "surface_areas", self.surface_areas, overwrite)
        for index, triangles in enumerate(self.triangles):
            writedataset(h5group, "triangles{}".format(index), np.array(triangles), overwrite)
        writedataset(h5group, "mass_centers", self.mass_centers, overwrite)
        writedataset(h5group, "squared_gyration_radii", self.squared_gyration_radii, overwrite)
        writedataset(h5group, "asphericities", self.asphericities, overwrite)
        writedataset(h5group, "acylindricities", self.acylindricities, overwrite)
        writedataset(h5group, "anisotropies", self.anisotropies, overwrite)
        writedataset(h5group, "characteristic_radii", self.characteristic_radii, overwrite)
        writedataset(h5group, "cyclic_area_indices", self.cyclic_area_indices, overwrite)

    def _export_gyration_parameters(self, fmt):
        mass_centers_file_name = fmt.format(property="centers_of_mass")
        squared_gyration_radii_file_name = fmt.format(property="squared_gyration_radii")
        asphericities_file_name = fmt.format(property="asphericities")
        acylindricities_file_name = fmt.format(property="acylindricities")
        anisotropies_file_name = fmt.format(property="anisotropies")
        characteristic_radii_file_name = fmt.format(property="characteristic_radii")

        mass_centers = self.getattr_normalized('mass_centers')
        squared_gyration_radii = self.getattr_normalized('squared_gyration_radii')
        asphericities = self.getattr_normalized('asphericities')
        acylindricities = self.getattr_normalized('acylindricities')
        anisotropies = self.getattr_normalized('anisotropies')
        characteristic_radii = self.getattr_normalized('characteristic_radii')

        export_filenames = []

        if mass_centers is not None:
            with open(mass_centers_file_name, 'w') as outfile:
                for index, mass_center in enumerate(mass_centers, start=1):
                    outfile.write("{} {}\n".format(index, mass_center))
            export_filenames.append(mass_centers_file_name)
        if squared_gyration_radii is not None:
            with open(squared_gyration_radii_file_name, 'w') as outfile:
                for index, squared_gyration_radius in enumerate(squared_gyration_radii, start=1):
                    outfile.write("{} {}\n".format(index, squared_gyration_radius))
            export_filenames.append(squared_gyration_radii_file_name)
        if asphericities is not None:
            with open(asphericities_file_name, 'w') as outfile:
                for index, asphericity in enumerate(asphericities, start=1):
                    outfile.write("{} {}\n".format(index, asphericity))
            export_filenames.append(asphericities_file_name)
        if acylindricities is not None:
            with open(acylindricities_file_name, 'w') as outfile:
                for index, acylindricity in enumerate(acylindricities, start=1):
                    outfile.write("{} {}\n".format(index, acylindricity))
            export_filenames.append(acylindricities_file_name)
        if anisotropies is not None:
            with open(anisotropies_file_name, 'w') as outfile:
                for index, anisotropy in enumerate(anisotropies, start=1):
                    outfile.write("{} {}\n".format(index, anisotropy))
            export_filenames.append(anisotropies_file_name)
        if characteristic_radii is not None:
            with open(characteristic_radii_file_name, 'w') as outfile:
                for index, characteristic_radius in enumerate(characteristic_radii, start=1):
                    outfile.write("{} {}\n".format(index, characteristic_radius))
            export_filenames.append(characteristic_radii_file_name)

        return export_filenames

    def getattr_normalized(self, attr):
        if hasattr(self, attr):
            value = getattr(self, attr)
            is_numpy_array = isinstance(value, np.ndarray)
            if (is_numpy_array and len(value.shape) > 0) or (not is_numpy_array and len(value) > 0):
                return value
        return None

class Domains(CavitiesBase):
    """
    Stores the calculated data about the domains.
    """

    def __init__(self, *args):
        """
        The constructor can be called in three ways:

        - ``Domains(timestamp, volumes, surface_areas, triangles, centers)`` :
            create the object using the given data

        - ``Domains(domaincalculation)`` :
            copy the data from this
            :class:`core.calculation.algorithm.DomainCalculation` object

        - ``Domains(hdf5group)`` :
            read the data from this hdf5 group
        """
        # Import this here to avoid cyclic imports
        from ..computation.cavity_calculation import algorithm as algorithm

        if isinstance(args[0], h5py.Group):
            super(Domains, self).__init__(*args)
            h5group = args[0]
            centers = h5group["centers"]
        elif isinstance(args[0], algorithm.DomainCalculation):
            calculation = args[0]
            timestamp = datetime.now()
            volumes = calculation.domain_volumes
            surface_areas = calculation.domain_surface_areas
            triangles = calculation.domain_triangles
            centers = calculation.centers
            discretization = calculation.discretization
            mass_centers = calculation.mass_centers
            squared_gyration_radii = calculation.squared_gyration_radii
            asphericities = calculation.asphericities
            acylindricities = calculation.acylindricities
            anisotropies = calculation.anisotropies
            characteristic_radii = calculation.characteristic_radii
            cyclic_area_indices = calculation.cyclic_area_indices
            self.critical_domains = calculation.critical_domains
            super(Domains, self).__init__(timestamp, volumes, surface_areas, triangles,
                                          mass_centers, squared_gyration_radii, asphericities, acylindricities,
                                          anisotropies, characteristic_radii, cyclic_area_indices)
        else:
            super(Domains, self).__init__(*args)
            centers = args[4]

        self.centers = np.array(centers, dtype=np.int32)
        if 'discretization' in locals():
            # TODO: get discretization also from other constructor calls!
            self.continuous_centers = np.array([discretization.discrete_to_continuous(center) for center in centers], dtype=np.float64)

    def tohdf(self, h5group, overwrite=True):
        """
        Write the data to a hdf5 Group.

        **Parameters:**
            `h5group` :
                the hdf5 group in which the data will be written

            `overwrite` :
                specifies if existing data should be overwritten
        """
        super(Domains, self).tohdf(h5group, overwrite)
        writedataset(h5group, "centers", self.centers, overwrite)

    def totxt(self, fmt):
        domain_center_file_name = fmt.format(property="centers")
        domain_surface_file_name = fmt.format(property="surface_areas")
        domain_volume_file_name = fmt.format(property="volumes")
        domain_surface_to_volume_ratio_file_name = fmt.format(property="surface_area_to_volume_ratios")

        export_filenames = [domain_center_file_name, domain_surface_file_name, domain_volume_file_name,
                            domain_surface_to_volume_ratio_file_name]

        with open(domain_surface_file_name, 'w') as outfile:
            for index, surface_area in enumerate(self.surface_areas, start=1):
                outfile.write("{} {}\n".format(index, surface_area))
        with open(domain_volume_file_name, 'w') as outfile:
            for index, volume in enumerate(self.volumes, start=1):
                outfile.write("{} {}\n".format(index, volume))
        with open(domain_surface_to_volume_ratio_file_name, 'w') as outfile:
            for index, t in enumerate(zip(self.volumes, self.surface_areas), start=1):
                volume, surface_area = t
                outfile.write("{} {}\n".format(index, surface_area/volume))
        continuous_centers = self.getattr_normalized('continuous_centers')
        if continuous_centers is not None:
            with open(domain_center_file_name, 'w') as outfile:
                for index, continuous_center in enumerate(continuous_centers, start=1):
                    outfile.write("{} {} {} {}\n".format(index, continuous_center[0], continuous_center[1], continuous_center[2]))
        else:
            raise ValueError('No discretization present -> can only access discrete domain centers, no conversion to continuous space possible.')

        export_filenames.extend(self._export_gyration_parameters(fmt))

        return export_filenames

class Cavities(CavitiesBase):
    """
    Stores the calculated data about the cavities.
    """

    def __init__(self, *args):
        """
        The constructor can be called in three ways:

        - ``Cavities(timestamp, volumes, surface_areas, triangles, multicavities)`` :
            create the object using the given data

        - ``Cavities(cavitycalculation)`` :
            copy the data from this
            :class:`core.calculation.algorithm.CavityCalculation` object

        - ``Cavities(hdf5group)`` :
            read the data from this hdf5 group
        """
        # Import this here to avoid cyclic imports
        from ..computation.cavity_calculation import algorithm as algorithm

        if isinstance(args[0], h5py.Group):
            super(Cavities, self).__init__(*args)
            h5group = args[0]
            multicavities = [None] * self.number
            for i in range(self.number):
                multicavities[i] = h5group["multicavities{}".format(i)]
        elif isinstance(args[0], algorithm.CavityCalculation):
            calculation = args[0]
            timestamp = datetime.now()
            volumes = calculation.multicavity_volumes
            surface_areas = calculation.cavity_surface_areas
            triangles = calculation.cavity_triangles
            multicavities = calculation.multicavities
            mass_centers = calculation.mass_centers
            squared_gyration_radii = calculation.squared_gyration_radii
            asphericities = calculation.asphericities
            acylindricities = calculation.acylindricities
            anisotropies = calculation.anisotropies
            characteristic_radii = calculation.characteristic_radii
            cyclic_area_indices = calculation.cyclic_area_indices
            super(Cavities, self).__init__(timestamp, volumes, surface_areas, triangles,
                                           mass_centers, squared_gyration_radii, asphericities, acylindricities,
                                           anisotropies, characteristic_radii, cyclic_area_indices)
        else:
            super(Cavities, self).__init__(*args)
            multicavities = args[4]

        self.multicavities = [None] * self.number
        for index, cavities in enumerate(multicavities):
            # cavities might be a 0-dimensional ndarray of python objects
            if isinstance(cavities, np.ndarray):
                cavities = cavities.tolist()
            self.multicavities[index] = np.array(list(cavities), dtype=np.int32)

    def tohdf(self, h5group, overwrite=True):
        """
        Write the data to a hdf5 Group.

        **Parameters:**
            `h5group` :
                the hdf5 group in which the data will be written

            `overwrite` :
                specifies if existing data should be overwritten
        """
        super(Cavities, self).tohdf(h5group, overwrite)
        for index, cavities in enumerate(self.multicavities):
            writedataset(h5group, "multicavities{}".format(index),
                         cavities, overwrite)

    def totxt(self, fmt):
        cavity_surface_file_name = fmt.format(property="surface_areas")
        cavity_volume_file_name = fmt.format(property="volumes")
        cavity_domains_file_name = fmt.format(property="domain_indices")
        cavity_surface_to_volume_ratio_file_name = fmt.format(property="surface_area_to_volume_ratios")

        export_filenames = [cavity_surface_file_name, cavity_surface_file_name, cavity_volume_file_name,
                            cavity_surface_to_volume_ratio_file_name]

        with open(cavity_domains_file_name, 'w') as outfile:
            for index, multicavity in enumerate(self.multicavities, start=1):
                outfile.write("{}".format(index))
                for domain_index in multicavity:
                    outfile.write(" {}".format(domain_index+1))
                outfile.write("\n".format(index))
        with open(cavity_surface_file_name, 'w') as outfile:
            for index, surface_area in enumerate(self.surface_areas, start=1):
                outfile.write("{} {}\n".format(index, surface_area))
        with open(cavity_volume_file_name, 'w') as outfile:
            for index, volume in enumerate(self.volumes, start=1):
                outfile.write("{} {}\n".format(index, volume))
        with open(cavity_surface_to_volume_ratio_file_name, 'w') as outfile:
            for index, t in enumerate(zip(self.volumes, self.surface_areas), start=1):
                volume, surface_area = t
                outfile.write("{} {}\n".format(index, surface_area/volume))

        export_filenames.extend(self._export_gyration_parameters(fmt))

        return export_filenames
                
class Results(object):
    """
    Container class to store the calculated putput data together with its
    input data. This can be passed to the Visualization:
    
    Attributes
    ----------
    filepath : string
        path to the input file
    frame : int
        frame number
    resolution : int
        resolution of the discretization
    atoms : list
        input data
    domains : class
        the calculated domains or `None`
    surface_cavities : class
        the calculated surface-based cavities or `None`
    center_cavities : class
        the calculated center-based cavities or `None`
    
    """
    def __init__(self, filepath, frame, resolution, atoms, 
                 domains=None, 
                 surface_cavities=None, 
                 center_cavities=None,
                 rings=None, polyhedra=None, 
                 config=None):
        """
        constructor

        Parameters
        ----------
        filepath : string
            path to the input file.
        frame : int
            frame number.
        resolution : int
            resolution of the discretization.
        atoms : list
            input data.
        domains : class
            the calculated domains. The default is None.
        surface_cavities : class
            the calculated surface-based cavities. The default is None.
        center_cavities : class
            the calculated center-based cavities. The default is None.
        rings : class
            the calculated rings object. The default is None.

        Returns
        -------
        None.

        """
        self.filepath = filepath
        self.frame = frame
        self.resolution = resolution
        self.atoms = atoms
        self.domains = domains
        self.surface_cavities = surface_cavities
        self.center_cavities = center_cavities
        self.rings = rings
        self.polyhedra = polyhedra
        self.config = config
        
    def __str__(self):
        return self.description()

    def description(self, domain_volume=True, surface_cavity_volume=True, center_cavity_volume=True):
        s = "{}, frame {}, resolution {}".format(
                os.path.basename(self.filepath),
                self.frame + 1,
                self.resolution)
        if surface_cavity_volume and self.surface_cavities is not None and self.atoms.volume is not None:
            cavvolume = np.sum(self.surface_cavities.volumes)
            volpercent = 100 * cavvolume / self.atoms.volume.volume
            s += ", {:0.1f}% cavities (surface-based)".format(volpercent)
        if center_cavity_volume and self.center_cavities is not None and self.atoms.volume is not None:
            cavvolume = np.sum(self.center_cavities.volumes)
            volpercent = 100 * cavvolume / self.atoms.volume.volume
            s += ", {:0.1f}% cavities (center-based)".format(volpercent)
        if domain_volume and self.domains is not None and self.atoms.volume is not None:
            cavvolume = np.sum(self.domains.volumes)
            volpercent = 100 * cavvolume / self.atoms.volume.volume
            s += ", {:0.1f}% cavities (domains)".format(volpercent)
        return s

    # Properties to be compatible to the old CalculationResults
    @property
    def number_of_atoms(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        return self.atoms.number

    @property
    def atom_positions(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        return self.atoms.positions

    @property
    def atom_radii(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        return self.atoms.radii

    @property
    def number_of_domains(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.domains is not None:
            return self.domains.number
        else:
            return None

    @property
    def domain_volumes(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.domains is not None:
            return self.domains.volumes
        else:
            return None

    @property
    def domain_surface_areas(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.domains is not None:
            return self.domains.surface_areas
        else:
            return None

    @property
    def domain_centers(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.domains is not None:
            return self.domains.centers
        else:
            return None

    @property
    def domain_triangles(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.domains is not None:
            return self.domains.triangles
        else:
            return None

    @property
    def number_of_multicavities(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.surface_cavities is not None:
            return self.surface_cavities.number
        else:
            return None

    @property
    def multicavity_volumes(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.surface_cavities is not None:
            return self.surface_cavities.volumes
        else:
            return None

    @property
    def multicavity_surface_areas(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.surface_cavities is not None:
            return self.surface_cavities.surface_areas
        else:
            return None

    @property
    def multicavities(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.surface_cavities is not None:
            return self.surface_cavities.multicavities
        else:
            return None

    @property
    def multicavity_triangles(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.surface_cavities is not None:
            return self.surface_cavities.triangles
        else:
            return None

    @property
    def number_of_center_multicavities(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.center_cavities is not None:
            return self.center_cavities.number
        else:
            return None

    @property
    def center_multicavity_volumes(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.center_cavities is not None:
            return self.center_cavities.volumes
        else:
            return None

    @property
    def center_multicavity_surface_areas(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.center_cavities is not None:
            return self.center_cavities.surface_areas
        else:
            return None

    @property
    def center_multicavities(self):
        if USE_LOGGER: 
            logger.warn("use of deprecated property")
        if self.center_cavities is not None:
            return self.center_cavities.multicavities
        else:
            return None

    @property
    def center_multicavity_triangles(self):
        if USE_LOGGER:
            logger.warn("use of deprecated property")
        if self.center_cavities is not None:
            return self.center_cavities.triangles
        else:
            return None
       
class ResultsFile(Results):
    def __init__(self, file, mode, atoms=None,
                 cavity=None, rings_guttman=None, rings_king=None, rings_primitive=None):
                
        self.mode = mode
        self._cavity = cavity
        
        if cavity is None:
            domains = None
            surface_cavities = None
            center_cavities = None
        else:
            domains = cavity.domains
            surface_cavities = cavity.surface_cavities
            center_cavities = cavity.center_cavities
            
        self.filepath = file
        self.frame = 0
        self.resolution = 0
        self.atoms = atoms
        self.domains = domains
        self.surface_cavities = surface_cavities
        self.center_cavities = center_cavities
        self.rings_guttman = rings_guttman
        self.rings_king = rings_king
        self.rings_primitive = rings_primitive
        self.polyhedra = None
        self.config = None
        
    def __enter__(self):
        if (self.mode == 'r') or (self.mode == 'a'):
            self.read()        
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if (self.mode == 'w') or (self.mode == 'a'):
            self.write()
    
    @property
    def cavity(self):
        return self._cavity
    
    @cavity.setter
    def cavity(self, value):
        self._cavity = value
        self.domains = self._cavity.domains
        self.surface_cavities = self._cavity.surface_cavities
        self.center_cavities = self._cavity.center_cavities
    
    def flush(self):
        self.write()
    
    def read(self):
        with h5py.File(self.filepath, "r") as f:

            self.name    = f["name"][0].decode()
            self.version = f["version"][0].decode()

            if 'atoms' in f:
                group = f['atoms']
                positions = np.array(group['positions'])
                radii = np.array(group['radii'])
                elements = [s.decode() for s in group['elements'][:]]
                volume = group['volume'][0].decode()
                origin = np.array(group['volume_origin'])                                
                if "bond_lengths" in group.keys():
                    bls = group["bond_lengths"]
                    bond_lengths = dict()
                    for s in bls:
                        bond_lengths[(s[0].decode(),s[1].decode())] = float(s[2])
                else:
                    bond_lengths = None
                    
                self.atoms = Atoms(positions, radii, elements, volume)
                self.atoms.volume.origin = origin
                self.atoms.set_bond_lengths(bond_lengths)

                # Load original text data
                if "original_data" in f:
                    group = f.require_group("original_data")
                    ori_name = group["file_name"][0].decode()
                    ori_text = group["structure_text"][0].decode()
                    self.atoms.original_file_data = OriginalStructureData(ori_name, ori_text)
            
            if 'rings_guttman' in f:
                group = f['rings_guttman']
                irings = np.array(group['indexes']).tolist()
                self.rings_guttman = []
                for iring in irings:
                    try:
                        i = iring.index(-1)
                        indexes = iring[:i]
                    except ValueError:
                        indexes = iring               
                    ring = Ring(self.atoms, indexes)
                    self.rings_guttman.append(ring)

            if 'rings_king' in f:
                group = f['rings_king']
                irings = np.array(group['indexes']).tolist()
                self.rings_king = []
                for iring in irings:
                    try:
                        i = iring.index(-1)
                        indexes = iring[:i]
                    except ValueError:
                        indexes = iring               
                    ring = Ring(self.atoms, indexes)
                    self.rings_king.append(ring)

            if 'rings_primitive' in f:
                group = f['rings_primitive']
                irings = np.array(group['indexes']).tolist()
                self.rings_primitive = []
                for iring in irings:
                    try:
                        i = iring.index(-1)
                        indexes = iring[:i]
                    except ValueError:
                        indexes = iring               
                    ring = Ring(self.atoms, indexes)
                    self.rings_primitive.append(ring)
                    
            if 'domains' in f:
                self._cavity = Cavity(self.atoms)
                group = f['domains']
                self._cavity.results = Results('test.xyz', 0, -1, self.atoms) #TODO
                self._cavity.results.domains = Domains(group)

                if "critical_domains" in f:
                    self._cavity.results.domains.critical_domains = list(f["critical_domains"])

                group_setting = f["cavity_setting"]
                self._cavity.resolution = np.int64(group_setting["resolution"])

                cutoff_radii = np.array(group_setting["cutoff_radii"])
                if cutoff_radii.size == 1:
                    if cutoff_radii == -1:
                        self._cavity.cutoff_radii = None
                    else:
                        self._cavity.cutoff_radii = cutoff_radii
                else:
                        self._cavity.cutoff_radii = dict()
                        for s in cutoff_radii:
                            self._cavity.cutoff_radii[s[0].decode()] = float(s[1])
                
            if 'surface_cavities' in f:
                # group = f['surface_cavities']
                # self.surface_cavities = Cavities(group)
                group = f['surface_cavities']
                self._cavity.results.surface_cavities = Cavities(group)
                
            if 'center_cavities' in f:
                group = f['center_cavities']
                self._cavity.results.center_cavities = Cavities(group)
            
    def write(self, overwrite=True):
        with h5py.File(self.filepath, "w") as f:

            f.create_dataset('name', data=["sovapy"])
            f.create_dataset('version', data=[VERSION])

            group = f.create_group("atoms")
            group['positions'] = self.atoms.positions
            group['radii'] = self.atoms.radii
            group['elements'] = self.atoms.elements.tolist()
            group['volume'] = [str(self.atoms.volume)]
            group['volume_origin'] = self.atoms.volume.origin.tolist()

            #Original text data
            if (self.atoms.original_file_data is not None):
                group = f.create_group('original_data')
                dt = h5py.string_dtype()
                ds_name = group.create_dataset('file_name', (1, ), dtype=dt, compression="gzip")
                ds_name[0] = self.atoms.original_file_data.file_name
                ds_txt = group.create_dataset('structure_text', (1, ), dtype=dt, compression="gzip")
                ds_txt[0] = self.atoms.original_file_data.structure_text
                
            print('-----bond_lengths in Atoms')
            print(self.atoms.bond_lengths)
            if self.atoms.bond_lengths is not None:
                list_bond_lengths = list()
                for a in self.atoms.bond_lengths.items():
                    list_bond_lengths.append(list(a[0])+[str(a[1])])
                group["bond_lengths"] = list_bond_lengths
            
            if self.rings_guttman is not None:
                group = f.create_group("rings_guttman")
                # get max size
                max_size = 0
                for ring in self.rings_guttman:
                    max_size = max(max_size, len(ring.indexes))
                # ring = self.rings_guttman[0]
                irings = []
                for ring in self.rings_guttman:
                    irings.append(ring.indexes + [-1]*(max_size-len(ring.indexes)))
                group['indexes'] = irings
            
            if self.rings_king is not None:
                group = f.create_group("rings_king")
                # get max size
                max_size = 0
                for ring in self.rings_king:
                    max_size = max(max_size, len(ring.indexes))
                # ring = self.rings[0]
                irings = []
                for ring in self.rings_king:
                    irings.append(ring.indexes + [-1]*(max_size-len(ring.indexes)))
                group['indexes'] = irings

            if self.rings_primitive is not None:
                group = f.create_group("rings_primitive")
                # get max size
                max_size = 0
                for ring in self.rings_primitive:
                    max_size = max(max_size, len(ring.indexes))
                # ring = self.rings[0]
                irings = []
                for ring in self.rings_primitive:
                    irings.append(ring.indexes + [-1]*(max_size-len(ring.indexes)))
                group['indexes'] = irings
                
            if self.domains is not None:
                group_setting = f.create_group("cavity_setting")
                self.resolution = self.cavity.resolution
                group_setting["resolution"] = self.cavity.resolution
                if self.cavity.cutoff_radii is None:
                    group_setting["cutoff_radii"] = -1
                elif type(self.cavity.cutoff_radii) is dict:
                    list_radii = list()
                    for k, d in self.cavity.cutoff_radii.items():
                        list_radii.append([k,str(d)])
                    group_setting["cutoff_radii"] = list_radii
                else:
                    group_setting["cutoff_radii"] = self.cavity.cutoff_radii

                if "critical_domains" in f:
                    del f["critical_domains"]
                f["critical_domains"] = self.domains.critical_domains

                group = f.create_group("domains")
                self.domains.tohdf(group, overwrite)
            if self.surface_cavities is not None:
                group = f.create_group("surface_cavities")
                self.surface_cavities.tohdf(group, overwrite)
            if self.surface_cavities is not None:
                group = f.create_group("center_cavities")
                self.center_cavities.tohdf(group, overwrite)
