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
 