# MAXWELL EQUATIONS SOLVER

This project is a Finite Elements Method (FEM) solver for Maxwell's equations.

<p align="center">
        <img src=equations.png />
</p>

## Requirements
To compile the library :
- CMake  
- A C++ compiler  
- The PyBind11 library  
- The Intel Threading Building Blocks (TBB) library

To generate the reference documentation :
- Doxygen  

To use it with Python :
- Python 3+  
- The Numpy and Matplotlib libraries  

NB: The TBB library is not included in this repository, but can be found here:  
https://github.com/oneapi-src/oneTBB

## How to install
1. Make sure that the TBB library is installed next this folder, as below:
```console
user@linux:~/path$ tree  
.  
+-- maxwell  
|   +-- doc
|   |   +-- ...
|   |
|   +-- python
|   |   +-- ...
|   |
|   +-- src
|   |   +-- ...
|   |
|   +-- equations.png  
|   +-- README.md  
|
+-- TBB
|   +-- build  
|   |   +-- ...
|   |
|   +-- cmake
|   |   +-- ...
|   ...  
```

If the library is installed elsewhere, the path to the library in the src/CMakeLists.txt file will have to be changed.  

2. Go to the code directory and execute the following commands to create a build directory to compile the library, a data directory to store the raw simulation outputs, and a results directory to save the postprocessed results.
```console
user@linux:~/path/maxwell$ mkdir build  
user@linux:~/path/maxwell$ mkdir data  
user@linux:~/path/maxwell$ mkdir results  
```

3. Execute the following commands to compile the C++ library and generate the reference documentation which will be available in the doc/html folder (file index.html).
```console
user@linux:~/path/maxwell$ cd build  
user@linux:~/path/maxwell/build$ cmake ../src -DCMAKE_INSTALL_PREFIX=$PWD/install    
user@linux:~/path/maxwell/build$ make install    
user@linux:~/path/maxwell/build$ make reference_doc
```

## How to use
```console
user@linux:~/path/maxwell$ cd python
user@linux:~/path/maxwell/python$ python simulation.py
```