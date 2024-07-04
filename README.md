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
  
## Usage

Basic usage:

```python
from sova.core.file import File

path = "./data/amorphous_md/a_SiO2_speed1e11K_rand.xyz"
f = File.open(path)
atoms = f.getatoms()
```

## Example
example1 : PDF analysis of the XYZ file with cell information (a_SiO2_speed1e11K_rand.xyz: atomic order randomized).
example2 : セル情報ありxyzファイル（a_SiO2_speed1e11K_rand.xyz）[Si,O]にソートしたPDF解析
example3 : セル情報なしxyzファイル（a_Si.xyz）のPDF解析
example4 : cifファイル（sio2_beta_cristobalite222.cif）のPDF解析
example5 : cfgファイル（sio.cfg）のPDF解析
example6 : 配位数計算
example7 : 角度分布解析(sio.cfgファイル)
example8 : 角度分布解析(a_Si.cifファイル)
example9 : 多面体解析(sio2_beta_cristobalite222.cif)
example10 : 多面体解析(a_SiO2_speed1e11K_rand.xyz)
example11 : ring解析
example12 : cavity解析
example13 : データの保存、読み込み