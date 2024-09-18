# Author: Motoki Shiga and Hidetoshi Morita
# email: motoki.shiga.b4@tohoku.ac.jp
#
# Copyright (c) 2024 
# Please see the LICENSE file for further information.


from setuptools import setup, find_packages
import re

DESCRIPTION = "SOVA (Structural Order Visualization and Analysis) with Python"
NAME = 'sovapy'
AUTHOR = 'Motoki Shiga & Hidetoshi Morita'
MAINTAINER = 'Motoki Shiga'
AUTHOR_EMAIL = 'motoki.shiga.b4@tohoku.ac.jp'
URL = 'https://github.com/MotokiShiga/sova-cui'
LICENSE = 'MIT License'
DOWNLOAD_URL = 'https://github.com/MotokiShiga/sova-cui'
PYTHON_REQUIRES = ">=3.10"
# Get the current version number:
with open('sovapy/__init__.py') as fd:
    VERSION = re.search("__version__ = '(.*)'", fd.read()).group(1)

INSTALL_REQUIRES = [
        'ase>=3.22.1',
        'h5py',
        'igraph>=0.11.2',
        'matplotlib',
        'networkx>=3.1',
        'numpy>=1.23.1',
        'tqdm',
        'scipy>=1.8',
        'spglib>=2.0',
        'PyCifRW>=4.4.5',
    ]

CLASSIFIERS = [
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Information Analysis',
    'Topic :: Scientific/Engineering :: Physics',
    'Topic :: Scientific/Engineering :: Chemistry',
    'Topic :: Software Development :: Libraries :: Python Modules',
]

setup(
    name=NAME,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    maintainer=MAINTAINER,
    maintainer_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    url=URL,
    version=VERSION,
    download_url=DOWNLOAD_URL,
    python_requires=PYTHON_REQUIRES,
    install_requires=INSTALL_REQUIRES,
    classifiers=CLASSIFIERS,
    include_package_data=True,
    packages=find_packages(),   
)
