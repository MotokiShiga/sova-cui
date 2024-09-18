# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 16:16:20 2022

@author: H. Morita
"""

import os
from .. import data
from .io_utils import get_abspath

class InputFile(object):
    """
    Abstract access to a file that contains atom data for one or more frames.
    Subclasses need to implement the ``readinfo()`` and ``readatoms(frame)`` methods.

    Attributes
    ----------
    id : int
        id.
    """
    def __init__(self, path):
        """
        Parameters
        ----------
        path : string
            absolute path to the file
        """
        self.path = get_abspath(path)
        self._info = data.FileInfo()        
        self.inforead = False

    @property
    def info(self):
        """
        :class:`core.data.FileInfo` object that contains metadata
        """        
        if not self.inforead:
            try:
                self.readinfo()
            except IOError as e:
                # logger.error(str(e))
                pass
            if self._info.volume is None:
                self._info.volume_guessed = True
                minx, maxx = float('inf'), float('-inf')
                miny, maxy = float('inf'), float('-inf')
                minz, maxz = float('inf'), float('-inf')
                for frame in range(self._info.num_frames):
                    atoms = self.getatoms(frame)
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

    def getatoms(self, frame=0,*args):
        """
        Read atom data for a specified frame.
        
        Parameters
        ----------
        frame : string
            the frame number
            
        Returns
        ----------
        an : object 
            class:`core.data.Atoms` object
        
        Raises
        ----------
        IndexError : 
            if the frame is not in the file
        FileError :
            if there are problems with the data in the file
        IOError : 
            if the file cannot be read
        """
        return self.readatoms(frame,*args)

    def readinfo(self):
        raise NotImplementedError

    def readatoms(self, frame):
        raise NotImplementedError