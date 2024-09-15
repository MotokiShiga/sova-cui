# sova-cui

This version of SOVA can be installed in macOS and Linux.    
We are preparing a binary version for Windows.

# Installation of python package **sova**
1. Clone the SOVA repo

  ```sh
  git clone https://github.com/***/sova-cui.git
  ```
  And move to the downloaded directory
  ```sh
  cd sova-cui
  bash 

2. To compile and to generate so or dll files, run  
for macos and linux
  ```sh
  bash run_install_mac_linux.sh
  ```
For windows, use "x64 Native Tools Command Prompt for Visual Studio 2022" to run 
  ```sh
  install_win.bat
  ```

3. To install SOVA, run    
  ```sh
  pip install .  
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
  
The versions of other packages can be found in requirements.txt.

You can make the virtual environment for sova by
```
python -m venv sova-cui
source sova-cui/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

## Acknowledgements

SOVA reuses source codes of the following package:

- [GR](https://github.com/sciapp/gr)
- [pyMolDyn](https://github.com/sciapp/pyMolDyn)
- [cif2cell](https://pypi.org/project/cif2cell/#description)


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

Execute example/1_pdf_xyz.py:

```sh
Histogram [O-O], [O-Si], [Si-Si]
```
<img src="docs/Figure_1_1.png" height=200 />

``` sh
Partial and Total g(r), S(Q), etc.
```
<img src="docs/Figure_1_2.png" height=400 />
