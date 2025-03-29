from .. import data
from .io_utils import get_abspath

class InputFile(object):
    """Abstract class to manage a structure data
    
    Abstract class to manage a structure data file that contains atom data for one or more frames.
    Sub-classes need methods ``read_info()`` and ``read_atoms(frame)``.

    Attributes
    ----------
    path : string
        Absolute path to a structure file
    _info : core.data.FileInfo object
        Meta data including #frames, unit cell info (volume)
    has_info : bool
        If the structure infomation has been already read
    """

    def __init__(self, path):
        """Initialization

        Parameters
        ----------
        path : string
            absolute path to a structure file
        """
        self.path = get_abspath(path)
        self._info = data.FileInfo()        
        self.has_info = False

    @property
    def info(self):
        """ Getter for meta data of the structure file
        
        Returns
        -------
        _info : core.data.FileInfo object
            Meta data of the structure file
        """
        # read meta data if it has not been read yet 
        if not self.has_info:
            try:
                self.read_info()
            except IOError as e:
                # logger.error(str(e))
                pass

            #generate unit cell data from the atomic configuration
            if self._info.volume is None:
                self._info.volume_guessed = True
                minx, maxx = float('inf'), float('-inf')
                miny, maxy = float('inf'), float('-inf')
                minz, maxz = float('inf'), float('-inf')
                for frame in range(self._info.num_frames):
                    atoms = self.get_atoms(frame)
                    minx = min(minx, atoms.positions[:, 0].min())
                    maxx = max(maxx, atoms.positions[:, 0].max())
                    miny = min(miny, atoms.positions[:, 1].min())
                    maxy = max(maxy, atoms.positions[:, 1].max())
                    minz = min(minz, atoms.positions[:, 2].min())
                    maxz = max(maxz, atoms.positions[:, 2].max())
                self._info.volumestr = 'ORT %f %f %f' % (maxx-minx, 
                                                         maxy-miny, 
                                                         maxz-minz)
        return self._info

    def get_atoms(self, frame=0, *args):
        """ Read atomic data
        Read atomic data in the specified frame.
        
        Parameters
        ----------
        frame : int
            The number of frames
            
        Returns
        ----------
        atoms : core.data.Atoms object 
            Atoms object that contains structural data
        
        Raises
        ----------
        IndexError : 
            If the frame is not in the file
        FileError :
            If there are problems with the data in the file
        IOError : 
            If the file cannot be read
        """
        return self.read_atoms(frame,*args)

    def read_info(self):
        raise NotImplementedError

    def read_atoms(self, frame):
        raise NotImplementedError