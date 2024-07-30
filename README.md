# sova-cui-dev



# Installation and Run
1. Clone the SOVA repo

  ```sh
  git clone https://github.com/***/sova-cui.git
  ```

2. Install the pip package    
  ```sh
  # Option 1: install from PyPI
  pip install sova
  
  # Option 2: install from source
  To enable verbose mode, use the -v option.
  pip install [-v] .
  ```
3. Run
  ```sh
  # Console
  > cd example
  > python example1.py
  
  # Jupyter notebook
  Please execute the files located in the 'example_jypyter' folder.
  ```    

## package

<!-- framework & version -->

| framework  | version |
| --------------------- | ---------- |
| Python                | 3.10.0     |
| numpy                | 1.23.5      |
| scipy | 1.8.1    |
| ase                 | 3.23.0     |
| h5py              | 3.7.0    |
| networkx                 | 3.1     |
| spglib               | 2.0.2     |
| PyCifRW             | 4.4.5      |
  
The versions of other packages can be found in the setup.py file.

## Usage

Basic usage:

```python
from sova.core.file import File

path = "./data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(path)
atoms = f.getatoms()
```

## Example
```sh
example1 : PDF analysis of the XYZ file with cell information.
example2 : PDF analysis of the XYZ file with cell information. sort [Si,O]
example3 : PDF analysis of the XYZ file with no cell information.
example4 : PDF analysis of the CIF file
example5 : PDF analysis of the cfg file (Reverse Monte Carlo file format)
example6 : Coordination number calculation
example7 : Angle distribution calculation (cfg file)
example8 : Angle distribution calculation (CIF file)
example9 : Polyhedra calculation (CIF file)
example10 : Polyhedra calculation XYZg file)
example11 : RINGs calculation
example12 : Cavity calculation
example13 : Data save and load.
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
