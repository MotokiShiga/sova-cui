# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 21:14:13 2022

@author: H. Morita
"""

import os

SEARCH_PATH = os.getcwd()

def get_abspath(path):
    """
    Return an absolute path for an input file. The difference to abspath from
    os.path is that this function uses SEARCH_PATH instead of the current
    working directory

    Parameters
    ----------
    path : string
        file path.

    Returns
    -------
    an absolute path : string
        an absolute path for an input file.
    """
    
    global SEARCH_PATH
    if not os.path.isabs(path) and SEARCH_PATH is not None:
        # relative paths use the SEARCH_PATH instead of the current working
        # directory, as that is modified at startup. The SEARCH_PATH is set
        # to the initial cwd in startBatch and startGUI.
        path = os.path.join(SEARCH_PATH, path)
    return os.path.abspath(path)


class FileError(Exception):
    """
    Exception which will be raised by the subclasses of :class:`InputFile`
    if errors occur when handling the data inside the file.
    """

    def __init__(self, message, cause=None):
        super(FileError, self).__init__(message, cause)
        self.message = message
        self.cause = cause

    def __str__(self):
        if self.cause is None:
            f = "{}: {}"
        else:
            f = "{}: {} (caused by: {})"
        return f.format(self.__class__.__name__, self.message, str(self.cause))
