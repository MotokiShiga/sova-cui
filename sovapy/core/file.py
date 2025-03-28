import os, h5py
from . import data
from .io.io_utils import get_abspath, FileError
from .io.input_file import InputFile
from .io.xyz_file import XYZFile
from .io.exyz_file import EXYZFile
from .io.cfg_file import CFGFile
from .io.cif_file import CIFFile

class ResultFile(InputFile):
    """
    Abstract access to a file that contains both input data (atoms) and
    calculated results.
    The `info` attribute has the type `ResultInfo`.
    In addition the `readinfo` and `readatoms` methods, subclasses must
    implement `writeinfo`, `readresults` and `writeresults`.
    """

    def __init__(self, path, sourcefilepath=None):
        """
        **Parameters:**
            `path` :
                absolute path to the file

            `sourcefilepath` :
                path to the file where the input data originally came from
        """
        super().__init__(path)
        self._info = data.ResultInfo()
        self._info.sourcefilepath = sourcefilepath

    def getresults(self, frame, resolution):
        """
        Read results from this file.

        **Parameters:**
            `frame` :
                the frame number
            `resolution` :
                the resolution of the calculation

        **Returns:**
            A `Results` object, if the file contains results for the
            specified parameters. Otherwise it return `None`.

        **Raises:**
            - :class:`FileError`: if there are problems with the data in the file
            - :class:`IOError`: if the file cannot be read
        """
        if not self.info[resolution].domains[frame] is None:
            return self.readresults(frame, resolution)
        else:
            return None

    def addresults(self, results, overwrite=True):
        """
        Write calculated results into this file.

        **Parameters:**
            `results` :
                the results to write
            `overwrite` :
                specifies if existing results should be overwritten

        **Raises:**
            - :class:`FileError`: if there are problems with the data in the file
            - :class:`IOError`: if the file cannot be read or written
        """
        self.writeresults(results, overwrite=overwrite)
        resinfo = self.info[results.resolution]
        if results.domains:
            resinfo.domains[results.frame] \
                = results.domains.timestamp
        if results.surface_cavities:
            resinfo.surface_cavities[results.frame] \
                = results.surface_cavities.timestamp
        if results.center_cavities:
            resinfo.center_cavities[results.frame] \
                = results.center_cavities.timestamp
        self.writeinfo()

    def writeinfo(self):
        raise NotImplementedError

    def readresults(self, frame, resolution):
        raise NotImplementedError

    def writeresults(self, results, overwrite=True):
        raise NotImplementedError

class HDF5File(ResultFile):
    """
    Implementation on :class:`ResultFile` for 'hdf5' files.
    """
    def __init__(self, path, sourcefilepath=None):
        super(HDF5File, self).__init__(path, sourcefilepath)

    @classmethod
    def fromInputFile(cls, filepath, sourcefilepath):
        """
        Create a new :class:`HDF5File`. Copy atoms from an
        :class:`InputFile`. This is used to initially create a file
        to export results to.
        """
        inputfile = File.open(sourcefilepath)
        outputfile = cls(filepath, sourcefilepath)
        for frame in range(inputfile.info.num_frames):
            atoms = inputfile.getatoms(frame)
            results = data.Results(filepath, frame, 64, atoms, None, None, None)
            outputfile.writeresults(results)            
        outputfile.readinfo()
        outputfile.writeinfo()
        return outputfile

    def readatoms(self, frame):
        atoms = None
        if not os.path.isfile(self.path):
            raise IOError(2, "File not found.")
        try:
            with h5py.File(self.path) as f:
                group = "atoms/frame{}".format(frame)
                if group not in f:
                    raise IndexError("Frame {} not found".format(frame))
                atoms = data.Atoms(f[group])
        except (IOError, IndexError):
            raise
        except Exception as e:
            raise FileError("Cannot read atom data.", e)
        return atoms

    def readinfo(self):
        if not os.path.isfile(self.path):
            print(self.path)
            raise IOError(2, "File not found.")
        try:
            with h5py.File(self.path) as f:
                if "info" in f:
                    info = data.ResultInfo(f["info"])
                    if self._info.sourcefilepath is not None:
                        info.sourcefilepath = self._info.sourcefilepath
                    self._info = info
                    self.inforead = True
        except IOError:
            raise
        except Exception as e:
            raise FileError("Cannot read file info.", e)

        if not self.inforead \
                and self._info.sourcefilepath is not None \
                and os.path.isfile(self._info.sourcefilepath):
            try:
                sf = File.open(self._info.sourcefilepath)
                self._info.num_frames = sf.info.num_frames
                self._info.volumestr = sf.info.volumestr
                self.inforead = True
            except IOError:
                raise
            except Exception as e:
                raise FileError("Cannot read file info.", e)

        if not self.inforead:
            raise RuntimeError("No File Info in this file and the source file.")

    def writeinfo(self):
        try:
            with h5py.File(self.path) as f:
                h5group = f.require_group("info")
                self.info.tohdf(h5group)
        except IOError:
            raise
        except Exception as e:
            raise FileError("Cannot write file info.", e)

    def readresults(self, frame, resolution):
        if not os.path.isfile(self.path):
            raise IOError(2, "File not found.")
        try:
            results = None
            with h5py.File(self.path) as f:
                groupname = "results/frame{}/resolution{}".format(frame, resolution)
                if groupname in f:
                    group = f[groupname]
                    atoms = data.Atoms(f["atoms/frame{}".format(frame)])
                    domains = data.Domains(group["domains"])
                    if "surface_cavities" in group:
                        surface_cavities = data.Cavities(group["surface_cavities"])
                    else:
                        surface_cavities = None
                    if "center_cavities" in group:
                        center_cavities = data.Cavities(group["center_cavities"])
                    else:
                        center_cavities = None
                    if self.info.sourcefilepath is not None:
                        filepath = self.info.sourcefilepath
                    else:
                        filepath = self.path
                    results = data.Results(filepath, frame, resolution,
                                           atoms, domains, surface_cavities,
                                           center_cavities)
        except IOError:
            raise
        except Exception as e:
            raise FileError("Cannot read results.", e)
        return results

    def writeresults(self, results, overwrite=True):
        # TODO: results valid?
        return  # debug
        try:
            with h5py.File(self.path) as f:
                group = f.require_group("atoms/frame{}".format(results.frame))
                # TODO: is it OK to never overwrite atoms?
                results.atoms.tohdf(group, overwrite=False)
                if results.domains is not None or \
                        results.surface_cavities is not None or\
                        results.center_cavities is not None:
                    group = f.require_group("results/frame{}/resolution{}".format(
                                            results.frame, results.resolution))
                if results.domains is not None:
                    subgroup = group.require_group("domains")
                    results.domains.tohdf(subgroup, overwrite=overwrite)
                if results.surface_cavities is not None:
                    subgroup = group.require_group("surface_cavities")
                    results.surface_cavities.tohdf(subgroup, overwrite=overwrite)
                if results.center_cavities is not None:
                    subgroup = group.require_group("center_cavities")
                    results.center_cavities.tohdf(subgroup, overwrite=overwrite)
        except IOError:
            raise
        except Exception as e:
            raise FileError("Cannot write results.", e)

class File(object):
    """ Class to manage structure data file

    Provides static methods for easy access to files and directories.
    The class attribute `types` associates filename endings with
    classes to handle them.

    Attributes
    ----------
    types : dictionary
        Pair of file format and file class
    
    """   
    types = dict([
                     ("xyz", XYZFile),
                     ("exyz", EXYZFile),
                     ("extxyz", EXYZFile),
                     ("cfg", CFGFile),
                     ("cif", CIFFile),
                     ("hdf5", HDF5File)
                 ])

    @classmethod
    def open(cls, file_path):
        """Open a structure data file

        Get an associated sub-class of `InputFile` according to the given file.

        Parameters
        ----------
        file_path : string
            Path to the file
            
        Returns
        ----------
        FileClass : core.io.InputFile object 
            An object of a subclass of core.io.InputFile.
        
        Raises
        ----------
        ValueError : 
            Flag if the file format is unknown.
        """        
        e = file_path.split(".")[-1]
        if e not in cls.types:
            raise ValueError("Unknown file format")
        
        # Get associated class according to the file extension.
        FileClass = cls.types[e]
        return FileClass(file_path)
        
    @classmethod
    def exists(cls, file_path):
        """Check file existance

        Check if a file exists and if it can be opened.

        Parameters
        ----------
        file_path : string
            Path to the file
            
        Returns
        ----------
        flag : bool 
            `True` if the file exists and there is a subclass of :class:`InputFile`
            associated with the filename ending.
        """
        
        file_path = get_abspath(file_path)
        name = os.path.basename(file_path)
        flag = os.path.isfile(file_path) and name.split(".")[-1] in cls.types
        return flag
    
    @classmethod
    def list_dir(cls, directory):       
        """List all files

        List all files in the directory.

        Parameters
        ----------
        directory : string
            path to a directory
            
        Returns
        ----------
        List : object 
            A list of filenames which can be opened.
        """
       
        if not directory:
            directory = "."
        return [f for f in os.list_dir(directory)
                if cls.exists(f)]
    
    # Deprecated
    # @classmethod
    # def new(cls, file_path):
    #     return NewFile(file_path)
    
#Deprecated
# class NewFile(InputFile):
#     def __init__(self, path):
#         super().__init__(path)

#     def readinfo(self):        
#         try:
#             self._info.num_frames = 0            
#             self._info.num_frames += 1
#             if self._info.num_frames == 1:                
#                 volume_info = 'CUB %f' % 1000.0
#             self._info.volumestr = volume_info
#             self.inforead = True
#         except IOError:
#             raise
#         except Exception as e:
#             raise FileError("Cannot read file info.", e)
    
#     def readatoms(self, frame):        
#         if self.info.num_frames <= frame:
#             raise IndexError("Frame {} not found".format(frame))
        
#         positions = []
#         symbols = []
#         positions.append([0.,0.,0.])
#         symbols.append('H')
#         return data.Atoms(positions, None, symbols, self.info.volume)