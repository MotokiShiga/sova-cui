# -*- coding: utf-8 -*-
"""
With the classes in this module the rather complicated calculation process
can be started with a simple method call.
Additionally, results are stored in a cache and can be reused later.
"""

import os
import core.data as data
import core.file
from core.file import File
from .algorithm import CavityCalculation, DomainCalculation, FakeDomainCalculation
from .discretization import DiscretizationCache, AtomDiscretization
from util import message
from hashlib import sha256
from config.configuration import config
from util.logger import Logger
import sys
import core.bonds
# from gui.dialogs.element_settings_dialog import ElementSettingsDialog
# from PySide6 import QtWidgets

__all__ = ["Calculation",
           "CalculationCache",
           "CalculationSettings",
           #"calculated",
           #"count_frames",
           #"calculated_frames",
           #"calculate_cavities",
           #"getresults",
           #"delete_center_cavity_information",
           #"timestamp",
           #"calculate"
           ]

logger = Logger("core.calculation")
logger.setstream("default", sys.stdout, Logger.WARNING)


class CalculationSettings(object):
    """
    Structure to store the parameters for one or more calculation.

    **Attributes:**
        `datasets` :
            Dictionary, which contains filenames (Strings) as keys.
            Each value is a list of Integers, which contains the frames.
            The value ``[-1]`` means 'all frames'.
        `resolution` :
            resolution parameter for the discretization
        `domains` :
            calculate cavitiy domains
        `surface_cavities` :
            calculate surface-based cavities
        `center_cavities` :
            calculate center-based cavities
        `recalculate` :
            results will be calculated even if cached results exists
        `exporthdf` :
            ``True`` if the results should be written into a hdf5 file.
            If ``False``, they are stored in the cache.
        `exporttext` :
            ``True`` if the results should be written into text files.
        `exportdir` :
            Directory for the exports. Only used, when at least one of
            the export parameters is ``True``. If `exportdir` is ``None``,
            the directory of the input file is used.
        `bonds` :
            calculate bonds
        `dihedral_angles` :
            calculate dihedral angles
    """

    def __init__(self,
                 datasets,
                 resolution=config.Computation.std_resolution,
                 cutoff_radii=config.Computation.std_cutoff_radius,
                 domains=False,
                 surface_cavities=False,
                 center_cavities=False,
                 gyration_tensor=False,
                 recalculate=False,
                 exporthdf5=False,
                 exporttext=False,
                 exportdir=None):
        """
        """
        self.datasets = datasets
        self.resolution = resolution
        self.cutoff_radii = cutoff_radii
        self.domains = domains
        self.surface_cavities = surface_cavities
        self.center_cavities = center_cavities
        self.gyration_tensor = gyration_tensor
        self.recalculate = recalculate
        self.exporthdf5 = exporthdf5
        self.exporttext = exporttext
        self.exportdir = exportdir
        self.bonds = False
        self.dihedral_angles = False

    def copy(self):
        """
        **Returns:**
            A deep copy of this object.
        """
        datasets = dict()
        for filename, frames in self.datasets.items():
            datasets[filename] = [f for f in frames]
        dup = self.__class__(datasets, self.resolution)
        dup.cutoff_radii = self.cutoff_radii
        dup.domains = self.domains
        dup.surface_cavities = self.surface_cavities
        dup.center_cavities = self.center_cavities
        dup.gyration_tensor = self.gyration_tensor
        dup.recalculate = self.recalculate
        dup.exporthdf5 = self.exporthdf5
        dup.exporttext = self.exporttext
        dup.exportdir = self.exportdir
        dup.bonds = self.bonds
        dup.dihedral_angles = self.dihedral_angles
        return dup


class Calculation(object):
    """
    This class provides the methods to start calculations.
    Optionally, results can be read from a cache.
    """

    def __init__(self, cachedir=None):
        """
        **Parameters:**
            `cachedir` :
                directory to store the cache files in; if none is given,
                a default one is used
        """
        if cachedir is None:
            cachedir = os.path.expanduser(config.Path.cache_dir)
        self.cachedir = cachedir
        max_cachefiles = config.Computation.max_cachefiles
        self.cache = CalculationCache(cachedir, max_cachefiles)

    def calculatedframes(self, filepath, resolution, surface=False, center=False):
        """
        Query the cache if it contains results for the given parameters.

        Parameters
        ----------
        filepath : string
            absolute path of the input file.
        resolution : int
            resolution of the used discretization.
        surface : boolean, optional
            query calculation results for surface-based cavities. The default is False.
        center : boolean, optional
            query calculation results for center-based cavities. The default is False.

        Returns
        -------
        class
            A :class:`core.data.TimestampList` containing the dates
            of the calculation or `None`.

        """
        info = None
        inputfile = File.open(filepath)
        # TODO: error handling
        if isinstance(inputfile, core.file.ResultFile):
            info = inputfile.info
        else:
            if filepath in self.cache:
                cf = self.cache[filepath]
                info = cf.info
        if info is not None:
            if surface:
                return info[resolution].surface_cavities
            elif center:
                return info[resolution].center_cavities
            else:
                return info[resolution].domains
        else:
            return data.TimestampList(inputfile.info.num_frames)

    def timestamp(self, filepath, frame, resolution, surface=False, center=False):
        """
        Query the cache if it contains results for the given parameters.

        **Returns:**
            The date and time when cached results were calculated
            on `None` if they do not exist.
        """
        calc = self.calculatedframes(filepath, resolution, surface, center)
        return calc[frame]

    def is_calculated(self, filepath, frame, resolution, surface=False, center=False):
        """
        **Returns:**
            If the cache contains results for the given parameters.
        """
        return not self.timestamp(filepath, frame,
                                  resolution, surface, center) is None

    def get_new_results(self, filepath, frame, resolution=None):
        inputfile = File.new(filepath)
        resolution = 64
        #results = data.Results(filepath, frame, resolution, inputfile.getatoms(frame), None, None, None)
        results = data.Results(filepath, frame, resolution, None, None, None, None)
        return results
    
    def get_results(self, filepath, resolution=None, surface=False, center=False):
        """
        Get cached results for the given parameters.

        Parameters
        ----------
        filepath : string
            absolute path of the input file.
        frame : int
            the frame number.
        resolution : int, optional
            resolution of the used discretization. The default is None.
        surface : class, optional
            query calculation results for surface-based cavities. The default is False.
        center : class, optional
            query calculation results for center-based cavities. The default is False.

        Returns
        -------
        results : class
            core.data.Results object if cached results exist. or None

        """
        inputfile = File.open(filepath)
        
        # elements settings in RMC cfg file 
        if isinstance(inputfile, core.file.CFGFile) == True:
            elem_file = os.path.splitext(filepath)[0] + '.elm'
            if os.path.exists(elem_file) == True:
                try:
                    elements = []
                    with open(elem_file.encode("utf-8"), 'r') as f:
                        for i in range(inputfile.info.nmol_types):
                            elements.append(f.readline().strip())
                    inputfile.elements = elements
                except FileNotFoundError:
                    print('File Not Found : ', elem_file)
            else:
                nums = {}
                for i in range(inputfile.info.nmol_types):
                    nums['X' + str(i+1)] = inputfile.info.ni[i]
                
                # dialog = ElementSettingsDialog(nums)
                # if dialog.exec() == QtWidgets.QDialog.Accepted:
                #     inputfile.elements = dialog.symbols()
                #     if dialog.is_savefile() == True:
                #         with open(elem_file, 'w') as f:
                #             for elem in inputfile.elements:
                #                 f.write('{}\n'.format(elem))
                # else:
                #     pass
        
        info = None
        # TODO: error handling
        resultfile = None
        results = None
        if isinstance(inputfile, core.file.ResultFile):
            resultfile = inputfile
        elif filepath in self.cache:
            resultfile = self.cache[filepath]
        if resultfile is not None:
            results = []
            for frame in range(inputfile.info.num_frames):
                if resolution is None:
                    resolutions = sorted(resultfile.info.resolutions())[::-1]
                    resolution = 64
                    for res in resolutions:                        
                        if resultfile.info[res].domains[frame] is not None:
                            resolution = res
                            break
                result = resultfile.getresults(frame, resolution)
                results.append(result)                
        if results is None:
            if resolution is None:
                resolution = 64
            info = inputfile.info
            results = []
            for frame in range(inputfile.info.num_frames):
                result = data.Results(filepath, frame, resolution, inputfile.getatoms(frame), None, None, None)
                results.append(result)
        return results, info

    def get_result(self, filepath, frame, resolution=None, surface=False, center=False):
        """
        Get cached results for the given parameters.

        Parameters
        ----------
        filepath : string
            absolute path of the input file.
        frame : int
            the frame number.
        resolution : int, optional
            resolution of the used discretization. The default is None.
        surface : class, optional
            query calculation results for surface-based cavities. The default is False.
        center : class, optional
            query calculation results for center-based cavities. The default is False.

        Returns
        -------
        results : class
            core.data.Results object if cached results exist. or None

        """
        inputfile = File.open(filepath)
        
        # elements settings in RMC cfg file 
        if isinstance(inputfile, core.file.CFGFile) == True:
            elem_file = os.path.splitext(filepath)[0] + '.elm'
            if os.path.exists(elem_file) == True:
                try:
                    elements = []
                    with open(elem_file.encode("utf-8"), 'r') as f:
                        for i in range(inputfile.info.nmol_types):
                            elements.append(f.readline().strip())
                    inputfile.elements = elements
                except FileNotFoundError:
                    print('File Not Found : ', elem_file)
            else:
                nums = {}
                for i in range(inputfile.info.nmol_types):
                    nums['X' + str(i+1)] = inputfile.info.ni[i]
                
                dialog = ElementSettingsDialog(nums)
                if dialog.exec_() == QtWidgets.QDialog.Accepted:
                    inputfile.elements = dialog.symbols()
                    if dialog.is_savefile() == True:
                        with open(elem_file, 'w') as f:
                            for elem in inputfile.elements:
                                f.write('{}\n'.format(elem))
                else:
                    pass
        
        # TODO: error handling
        resultfile = None
        results = None
        if isinstance(inputfile, core.file.ResultFile):
            resultfile = inputfile
        elif filepath in self.cache:
            resultfile = self.cache[filepath]
        if resultfile is not None:
            if resolution is None:
                resolutions = sorted(resultfile.info.resolutions())[::-1]
                resolution = 64
                for res in resolutions:
                    if resultfile.info[res].domains[frame] is not None:
                        resolution = res
                        break
            results = resultfile.getresults(frame, resolution)
        if results is None:
            if resolution is None:
                resolution = 64
            results = data.Results(filepath, frame, resolution, inputfile.getatoms(frame), None, None, None)
        return results
    
    def calculate_result(self,results, frame, resolution, cutoff_radii=None, domains=False, surface=False, center=False,
                       atoms=None, gyration_tensor_parameters=False, recalculate=False, last_frame=True):
        
        filepath = results.filepath
        
        # always recalculate if gyration tensor parameters shall be computed for center or surface based cavities
        recalculate = recalculate or (gyration_tensor_parameters and (center or surface))
        message.progress(0)
        
        #if atoms is None:
        #    atoms = inputfile.getatoms(frame)
        atoms.radii = cutoff_radii
        volume = atoms.volume                
        
        if recalculate:
            results.domains = None
            results.surface_cavities = None
            results.center_cavities = None
        
        if not ((domains and results.domains is None)
                or (surface and results.surface_cavities is None)
                or (center and results.center_cavities is None)):
            message.print_message("Reusing results")
        else:
            cachepath = os.path.join(self.cachedir, 'discretization_cache.hdf5')
            discretization_cache = DiscretizationCache(cachepath)
            with DiscretizationCache(cachepath) as discretization_cache:
                discretization = discretization_cache.get_discretization(volume, resolution)
            atom_discretization = AtomDiscretization(atoms, discretization)
            message.progress(10)
            if (domains and results.domains is None) \
                    or (surface and results.surface_cavities is None):
                # CavityCalculation depends on DomainCalculation
                message.print_message("Calculating domains")
                domain_calculation = DomainCalculation(discretization, atom_discretization)
                if domain_calculation.critical_domains:
                    logger.warn('Found {:d} critical domains in file {}, frame {:d}. Domain indices: {}'.format(
                        len(domain_calculation.critical_domains), os.path.basename(filepath), frame,
                        domain_calculation.critical_domains
                    ))
                    message.log('Found {:d} critical domains in file {}, frame {:d}'.format(
                        len(domain_calculation.critical_domains), os.path.basename(filepath), frame + 1,
                    ))
            if results.domains is None:
                results.domains = data.Domains(domain_calculation)
            message.progress(40)

            if surface and results.surface_cavities is None:
                message.print_message("Calculating surface-based cavities")
                cavity_calculation = CavityCalculation(domain_calculation, use_surface_points=True,
                                                       gyration_tensor_parameters=gyration_tensor_parameters)
                results.surface_cavities = data.Cavities(cavity_calculation)
            message.progress(70)

            if center and results.center_cavities is None:
                message.print_message("Calculating center-based cavities")
                domain_calculation = FakeDomainCalculation(discretization, atom_discretization, results)
                cavity_calculation = CavityCalculation(domain_calculation, use_surface_points=False,
                                                       gyration_tensor_parameters=gyration_tensor_parameters)
                results.center_cavities = data.Cavities(cavity_calculation)
            ### debug
            ######resultfile.addresults(results, overwrite=recalculate)  

        message.progress(100)
        message.print_message("Calculation finished")
        if last_frame:
            message.finish()

    def calculate_frame(self, filepath, frame, resolution, cutoff_radii=None, domains=False, surface=False, center=False,
                       atoms=None, gyration_tensor_parameters=False, recalculate=False, last_frame=True):
        """
        Get results for the given parameters. They are either loaded from the
        cache or calculated.

        Parameters
        ----------
        filepath : string
            absolute path of the input file.
        frame : int
            the frame number.
        resolution : int
            resolution of the used discretization.
        cutoff_radii : list, optional
            dict that maps element symbols to cutoff radii. The default is None.
        domains : class, optional
            calculate cavitiy domains. The default is False.
        surface : class, optional
            calculate surface-based cavities. The default is False.
        center : class, optional
            calculate center-based cavities. The default is False.
        atoms : class, optional
            core.data.Atoms object. The default is None.
        gyration_tensor_parameters : class, optional
            gyration tensor parameters will be calculated for cavities (they are always calculated for
            cavity domains). The default is False.
        recalculate : bool, optional
            results will be calculated even if cached results exists. The default is False.
        last_frame : bool, optional
            last frame flag. The default is True.

        Returns
        -------
        results : class
            core.data.Results object.

        """
        # always recalculate if gyration tensor parameters shall be computed for center or surface based cavities
        recalculate = recalculate or (gyration_tensor_parameters and (center or surface))
        message.progress(0)
        inputfile = File.open(filepath)
        # TODO: error handling
        if isinstance(inputfile, core.file.ResultFile):
            resultfile = inputfile
        else:
            resultfile = self.cache[filepath]
        try:
            results = resultfile.getresults(frame, resolution)
        except Exception as e:
            logger.debug("error in resultfile.getresults: {}".format(e))
            results = None
        
        if atoms is None:
            atoms = inputfile.getatoms(frame)
        atoms.radii = cutoff_radii
        volume = atoms.volume
                
        if results is None:
            results = data.Results(filepath, frame, resolution, atoms, None, None, None)
        
        if recalculate:
            results.domains = None
            results.surface_cavities = None
            results.center_cavities = None

        if not ((domains and results.domains is None)
                or (surface and results.surface_cavities is None)
                or (center and results.center_cavities is None)):
            message.print_message("Reusing results")
        else:
            cachepath = os.path.join(self.cachedir, 'discretization_cache.hdf5')
            discretization_cache = DiscretizationCache(cachepath)
            with DiscretizationCache(cachepath) as discretization_cache:
                discretization = discretization_cache.get_discretization(volume, resolution)
            atom_discretization = AtomDiscretization(atoms, discretization)
            message.progress(10)
            if (domains and results.domains is None) \
                    or (surface and results.surface_cavities is None):
                # CavityCalculation depends on DomainCalculation
                message.print_message("Calculating domains")
                domain_calculation = DomainCalculation(discretization, atom_discretization)
                if domain_calculation.critical_domains:
                    logger.warn('Found {:d} critical domains in file {}, frame {:d}. Domain indices: {}'.format(
                        len(domain_calculation.critical_domains), os.path.basename(filepath), frame,
                        domain_calculation.critical_domains
                    ))
                    message.log('Found {:d} critical domains in file {}, frame {:d}'.format(
                        len(domain_calculation.critical_domains), os.path.basename(filepath), frame + 1,
                    ))
            if results.domains is None:
                results.domains = data.Domains(domain_calculation)
            message.progress(40)

            if surface and results.surface_cavities is None:
                message.print_message("Calculating surface-based cavities")
                cavity_calculation = CavityCalculation(domain_calculation, use_surface_points=True,
                                                       gyration_tensor_parameters=gyration_tensor_parameters)
                results.surface_cavities = data.Cavities(cavity_calculation)
            message.progress(70)

            if center and results.center_cavities is None:
                message.print_message("Calculating center-based cavities")
                domain_calculation = FakeDomainCalculation(discretization, atom_discretization, results)
                cavity_calculation = CavityCalculation(domain_calculation, use_surface_points=False,
                                                       gyration_tensor_parameters=gyration_tensor_parameters)
                results.center_cavities = data.Cavities(cavity_calculation)
            ### debug
            ######resultfile.addresults(results, overwrite=recalculate)  

        message.progress(100)
        message.print_message("Calculation finished")
        if last_frame:
            message.finish()
        return results

    def calculate_results(self, results, frame, calcsettings):        
        self.calculate_result(
            results,
            frame,
            calcsettings.resolution,
            calcsettings.cutoff_radii,
            domains=calcsettings.domains,
            surface=calcsettings.surface_cavities,
            center=calcsettings.center_cavities,
            atoms=results.atoms,
            gyration_tensor_parameters=calcsettings.gyration_tensor,
            recalculate=calcsettings.recalculate,
            last_frame=True)
        
    def calculate(self, calcsettings):
        """
        Calculate (or load from the cache) all results for the given settings.

        Parameters
        ----------
        calcsettings : class
            `CalculationSettings` object.

        Returns
        -------
        allresults : list
            A list of list of :class:`core.data.Results` objects. The outer list contains
            an entry for each entry in `calcsettings.filenames`; the inner
            list has a `Results` entry for each frame specified in
            `calcsettings.frames`.

        """
        allresults = []
        for filename, frames in calcsettings.datasets.items():
            filepath = core.file.get_abspath(filename)
            fileprefix = os.path.basename(filename).rsplit(".", 1)[0]
            
            if calcsettings.exportdir is not None:
                exportdir = core.file.get_abspath(calcsettings.exportdir)
                # replace asterisks with directories
                dirlist = os.path.dirname(filepath).split("/")
                while "*" in exportdir:
                    i = exportdir.rindex("*")
                    exportdir = os.path.join(exportdir[:i], dirlist.pop() + exportdir[i+1:])
                if (calcsettings.exporthdf5 or calcsettings.exporttext) \
                        and not os.path.exists(exportdir):
                    os.makedirs(exportdir)
            else:
                exportdir = os.path.dirname(filepath)
                
            if calcsettings.exporthdf5:
                efpath = os.path.join(exportdir, fileprefix + ".hdf5")
                efpath = core.file.get_abspath(efpath)
                # copy atoms into HDF5 file
                exportfile = core.file.HDF5File.fromInputFile(efpath, filepath)
                # use HDF5 file as input
                filepath = efpath

            fileresults = []
            if frames[0] == -1:
                inputfile = File.open(filepath)
                frames = range(inputfile.info.num_frames)
                
            last_frame = False
            for frame in frames:            
                # calculate single frame
                if frame is frames[-1]:
                    last_frame = True
                frameresult = self.calculate_frame(
                    filepath,
                    frame,
                    calcsettings.resolution,
                    calcsettings.cutoff_radii,
                    domains=calcsettings.domains,
                    surface=calcsettings.surface_cavities,
                    center=calcsettings.center_cavities,
                    gyration_tensor_parameters=calcsettings.gyration_tensor,
                    recalculate=calcsettings.recalculate,
                    last_frame=last_frame)
                # export to text file
                
                if calcsettings.exporttext:
                    fmt = os.path.join(exportdir, fileprefix) + "-{property}-{frame:06d}.txt"
                    if frameresult.atoms is not None:
                        frameresult.atoms.totxt(fmt.format(property='{property}', frame=frame+1))
                    if frameresult.domains is not None:
                        try:
                            frameresult.domains.totxt(fmt.format(property='domain_{property}', frame=frame+1))
                        except ValueError as e:
                            logger.warn(e.message)
                            logger.warn('The export of domain information could not be finished.')
                    if frameresult.surface_cavities is not None:
                        frameresult.surface_cavities.totxt(fmt.format(property='surface_cavities_{property}', frame=frame+1))
                    if frameresult.center_cavities is not None:
                        frameresult.center_cavities.totxt(fmt.format(property='center_cavities_{property}', frame=frame+1))
                    # TODO: try/except
                # gather results
                fileresults.append(frameresult)
            allresults.append(fileresults)
        return allresults


class CalculationCache(object):
    """
    Stores calculation results. Associates the input file with a
    'hdf5' file containing the calculated results.
    This is realized with a single directory containing hdf5 files that
    are named after the SHA256 value of the absolute path of the input file.
    Additionally, an index file "index.txt" is created, which contains
    the input file paths and the resulting hdf5 file names.
    """

    # TODO: replacement strategy
    def __init__(self, directory, max_cachefiles=0):
        """
        constructor

        Parameters
        ----------
        directory : string
            path of the directory in which the cahce files are stored.
        max_cachefiles : int, optional
            max files of chache. The default is 0.

        Returns
        -------
        None.

        """
        self.directory = directory
        self.max_cachefiles = max_cachefiles
        self.index = dict()
        self.indexfilepath = self.abspath("index.txt")
        if not os.path.isdir(directory):
            os.mkdir(directory)
        self.buildindex()
        self.writeindex()

    def __contains__(self, filepath):
        """
        **Parameters:**
            `filepath` :
                path to the input file

        **Returns:**
            If a cache file for the given input file exists.
        """
        sourcefilepath = core.file.get_abspath(filepath)
        cachefilepath = self.abspath(self.cachefile(sourcefilepath))
        return os.path.isfile(cachefilepath)

    def __getitem__(self, filepath):
        """
        Get the cache file for a given input file.

        **Parameters:**
            `filepath` :
                path to the input file

        **Returns:**
            A :class:`core.file.HDF5File` object.
            If no cache file exist for the input file, a new one is created.
        """
        sourcefilepath = core.file.get_abspath(filepath)
        cachefilepath = self.abspath(self.cachefile(sourcefilepath))
        if sourcefilepath not in self.index:
            self.index[sourcefilepath] = cachefilepath
            self.writeindex()
        cachefile = core.file.HDF5File(cachefilepath, sourcefilepath)
        # TODO: what if not sourcefilepath == cachefile.info.sourcefilepath
        return cachefile

    def abspath(self, filename):
        return os.path.abspath(os.path.join(core.file.get_abspath(self.directory), filename))

    def cachefile(self, filepath):
        return sha256(filepath.encode('utf-8')).hexdigest() + ".hdf5"

    def buildindex(self):
        self.index = dict()
        filenames = set(f for f in os.listdir(self.directory) if os.path.splitext(f)[1] == "hdf5")
        filenames.discard('discretization_cache.hdf5')
        for filename in filenames:
            cachepath = self.abspath(filename)
            try:
                cachefile = core.file.HDF5File(cachepath)
                filepath = cachefile.info.sourcefilepath
                if filepath is not None \
                        and filename == self.cachefile(filepath):
                    self.index[filepath] = filename
            except (IOError, AttributeError, RuntimeError):
                pass

    def cleanindex(self):
        if 0 < self.max_cachefiles < len(self.index):
            # get list of cachefiles and their last modification times (mtime)
            file_mtimes = []
            for filepath, cachefile in sorted(self.index.items()):
                cachepath = os.path.join(self.directory, cachefile)
                cache_mtime = os.path.getmtime(cachepath)
                file_mtimes.append((cache_mtime, cachepath))
            # delete older cachefiles so that only max_cachefiles are left
            file_mtimes.sort()
            for _, cachepath in file_mtimes[:-self.max_cachefiles]:
                os.remove(cachepath)
            # rebuild the index
            self.buildindex()

    def writeindex(self):
        self.cleanindex()
        with open(self.indexfilepath, "w") as f:
            for filepath, cachefile in sorted(self.index.items(),
                                              key=lambda x: x[0]):
                #print >>f, filepath + "; " + cachefile
                print(filepath + "; " + cachefile, file=f)

calculation = Calculation()

