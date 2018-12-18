# Requirements

- To run our Hilbert-Join your hardware needs to support AVX-512 instructions.
- GNU compiler version >= 5.1
 
### AVX-512 instructions

- If we apply our Hilbert-curve, auto-vectorization will get eliminated, nevertheless we belive that future compilers will profit from the locality assumptions in the Hilbert curve.

To explicitly ensure, that CMake will use the GNU compiler use:

```{bash, engine='sh'}
export CXX=g++
export CC=gcc
```

# Build with CMake

to build this project you need to type the following commands into your shell:

```{bash, engine='sh'}
git clone https://gitlab.cs.univie.ac.at/martinp16cs/hilbertJoin.git
cd hilbertJoin
mkdir build
cd build
cmake ..
make -j
```

# Parameters

There are two different parameters KBLOCK and STRIPES, where we tried out different strategies. 
In our final experiments (see paper) we _always_ use the following setting:

- KBLOCK=2
- STRIPES=14

