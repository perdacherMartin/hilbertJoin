Accepted for [SIGMOD-2019](http://sigmod2019.org/sigmodcfp) in Amsterdam from 30th of June to 5th of July 2019.

# Requirements

- To run our Hilbert-Join your hardware needs to support AVX-512 instructions.
- GNU compiler version >= 5.1

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
export KBLOCK=4
export NUM_THREADS=_MAX_NUMBER_OF_CORES_
cmake ..
make -j
```

# Parameters

There are two different parameters KBLOCK and STRIPES, where we tried out different strategies. 
In our experiments (see paper) we _always_ use the following setting:

- KBLOCK=4
- STRIPES=14

For uniform data we suggest to use the following parameter settings:

- KBLOCK=16
- STRIPES=1


# Freqeuntly Asked Quesionts

[see FAQ](FAQ.md)