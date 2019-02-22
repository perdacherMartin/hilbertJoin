Accepted for [SIGMOD-2019](http://sigmod2019.org/sigmodcfp) in Amsterdam from 30th of June to 5th of July 2019.

# Requirements

- To run our Hilbert-Join your hardware needs to support AVX-512 instructions.
- GNU compiler version >= 5.1
- cmake version >= 3.7.0
- Linux package: *build-essential*, including *GNU make* version >= 4.1 

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

# Datasets used in our publication

#### Synthetic data (as comma seperated files)

- [Uniform_200K](https://ucloud.univie.ac.at/index.php/s/LaPLUmXQKsldvcO)
- [Uniform_select](https://ucloud.univie.ac.at/index.php/s/pUPFeZDXGtGNEpa)

#### Real data
- [BigCross](https://ucloud.univie.ac.at/index.php/s/ITlAQkZfIGFTvTD)
- [Activity recognition]() from UCI Data Repository
- [Higgs]() from UCI Data Repository
- [IoT botnet]() from UCI Data Repository

# Issues

Feel free to report [issues](https://gitlab.cs.univie.ac.at/martinp16cs/hilbertJoin/issues) about the code.

# Freqeuntly Asked Quesionts

[see FAQ](FAQ.md)