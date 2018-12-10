# Requirements

- To run our Hilbert-Join your hardware needs to support AVX-512 instructions.
- GNU compiler version >= 5.1

To explicitly ensure, that CMake will use the GNU compiler use:

```{console, engine='sh'}
export CXX=g++
export CC=gcc
```

# Build with CMake

to build this project you need to type the following commands into your shell:

```{console, engine='sh'}
git clone https://gitlab.cs.univie.ac.at/martinp16cs/hilbertJoin.git
cd hilbertJoin
mkdir build
cd build
cmake ..
make -j
```