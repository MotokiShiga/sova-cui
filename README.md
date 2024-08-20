# sova-cui-dev



# Installation
1. Clone the SOVA repo

  ```sh
  git clone https://github.com/***/sova-cui.git
  ```

2. Install the pip package    
  ```sh
  # Option 1: install from PyPI
  pip install sova
  
  # Option 2: install from source
  # To enable verbose mode, use the -v option.
  pip install [-v] .
  ```

For the usage, see example codes in the directory 'examples'.

## Requirements

<!-- framework & version -->

| framework  | version |
| --------------------- | ---------- |
| Python                | 3.10.0     |
| numpy                | 1.23.5      |
| scipy | 1.8.1    |
| ase                 | 3.23.0     |
| h5py              | 3.7.0    |
| networkx                 | 3.1     |
| igraph                 | 0.11.3     |
| spglib               | 2.0.2     |
| PyCifRW             | 4.4.5      |
  
The versions of other packages can be found in the setup.py file.

## Acknowledgements

SOVA reuses source codes of the following package:

- [pyMolDyn](https://github.com/sciapp/pyMolDyn)
- [cif2cell](https://pypi.org/project/cif2cell/#description)

## Usage

Basic usage:

```python
from sova.core.file import File

path = "./data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(path)
atoms = f.getatoms()
```

## Examples
```sh
1_pdf_xyz       : PDF analysis from xyz file (amorphous SiO2)
2_pdf_cif       : PDF analysis from xyz file (beta-cristobalite)
3_pdf_cfg       : PDF analysis from cfg file generated by RMC++ (amorphous SiO2)
4_coordination  : Coordination number analysis (amorphous SiO2)
5_bond_angle    : Bond angle analysis (amorphous SiO2)
6_polyhedra     : Polyhedral symmetry analysis (q-value) (amorphous SiO2)
7_ring          : Ring analysis  (beta-cristobalite)
8_cavity        : Cavity analysis (amorphous SiO2)
9_save_result   : Save and load calculated results
``` 

Execute example1.py:

```sh
Histogram [O-O], [O-Si], [Si-Si]
```
<img src="docs/Figure_1_1.png" height=200 />

``` sh
Partial and Total g(r), S(Q), etc.
```
<img src="docs/Figure_1_2.png" height=400 />
