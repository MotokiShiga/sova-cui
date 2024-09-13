# install external library



# Installation
1. Clone the GR repo

  ```sh
  git clone https://github.com/sciapp/gr.git -b v0.73.7 --depth 1 
  ```

2. Build the GR
```sh
cd gr
mkdir build
cd build
cmake ..
cmake --build .
```

3. Copy a source file 
```
cd ..
cp lib/gr3/gr3_mc.c ../../sova/computation/cavity_calculation/extension/
```

4. Build libalgorithm.so
```
make
```