# SOVApy (Structural Order Viaualization and Analysis with Python)

SOVA can be installed in Windows, macOS and Linux.    
(The package name to be imported is "sovapy".)

# Install from PIPY 

```sh
pip install sovapy
```

# Build and Install sovapy
1. Clone the SOVA repo

  ```sh
  git clone https://github.com/MotokiShiga/sova-cui.git
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
If you encounter an issue with "vcruntime140.dll", 
please install "Microsoft Visual C++ 2015 Redistributable".

3. To install SOVA, run    
  ```sh
  pip install .  
  ```

For the usage, see example codes in the directory 'examples'.

## Environment

Major packages used for our development

<!-- framework & version -->
| Package  | Version |
| --------------------- | ---------- |
| Python                | 3.10.0     |
| ase                 | 3.22.1     |
| h5py              | 3.7.0    |
| igraph                 | 0.11.3     |
| matplotlib   | 3.6.3  |
| networkx                 | 3.1     |
| numpy                | 1.23.5      |
| PyCifRW             | 4.4.5      |
| scipy | 1.8.1    |
| spglib               | 2.0.2     |

  
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
