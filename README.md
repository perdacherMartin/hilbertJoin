Paper: ["Cache-oblivious High-performance Similarity Join" - Martin Perdacher (University of Vienna); Claudia Plant (University of Vienna); Christian Böhm (Ludwig-Maximilians-Universität)](http://dx.doi.org/10.1145/3299869.3319859)
at [SIGMOD-2019](http://sigmod2019.org/sigmodcfp) in Amsterdam from 30th of June to 5th of July 2019.

# Reproducability

We want to make our experiments transparent and comprehensible. Code and experimental data for almost all of our figures is available at our [cloud-repository](https://ucloud.univie.ac.at/index.php/s/09fYaJlKclk0Iq8).

# Requirements

- To run our Hilbert-Join your hardware needs to support AVX-512 instructions.
- GNU compiler version >= 5.1
- cmake version >= 3.7.0
- git version >= 1.8.3.1
- Linux package: *build-essential*, including *GNU make* version >= 4.1 

### Random number generators
- We use the random number generator provided by Intel&copy; MKL. Therefore, a working [Intel&copy; MKL](https://software.intel.com/en-us/mkl) environment should be installed. Ensure, that the environment variable `$MKLROOT` [is set correctly](https://software.intel.com/en-us/mkl-linux-developer-guide-scripts-to-set-environment-variables).


# Before compilation

To explicitly ensure, that CMake will use the GNU compiler use:

```{bash, engine='sh'}
export CXX=g++
export CC=gcc
```

Lookup the [compiler-flag](https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html) for your hardware. Change the `-march` flag in your `CMakeLists.txt` depending on the hardware.

Example configuration for Skylake processors:
```{bash, engine='sh'}
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -march=skylake -ffast-math -fassociative-math -O3 -fopenmp -lmkl_core -lmkl_intel_lp64 -lmkl_intel_thread -liomp5")
```

If your hardware does not support AVX-512, you could use [Intel&copy; Software Development Emulator (SDE)](https://software.intel.com/en-us/articles/intel-software-development-emulator) to emulate AVX-512 registers.

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

# Example calls

### Self-join

For a selfjoin with random generated uniform data [0.0, 1.0):
`./hilbertSelfJoinCardinality -n 200000 -e 0.2 -d 64 -t 64`

- `-n` are the number of objects in set A
- `-e` epsilon
- `-d` number of features (or dimensions)
- `-t` number of threads
 
For a selfjoin with a dataset from a file:
`./hilbertSelfJoinCardinality -n 200000 -e 0.2 -d 64 -t 64 -f uniform_200000x64.csv`

- `-f` filename 
    Each value is separated by a comma ',' and has _d_ objects in each line. The file has _n_ lines without a header.
    You could also use a binary format ".bin". 

### Join

Join between two sets `A` and `B` with random generated uniform data [0.0, 1.0):
`./hilbertJoinCardinality -n 200000 -m 200000 -e 0.2 -d 20 -t 64`

where 
- `-n` are the number of objects in set A
- `-m` are the number of objects in set B
 
and files could be specified with

- `-f` file for set A
- `-g` file for set B

# Datasets used in our publication

Note: use `.csv` files without header!

#### Synthetic data (as comma seperated files)

- [Uniform_200K](https://ucloud.univie.ac.at/index.php/s/LaPLUmXQKsldvcO)
- [Uniform_select](https://ucloud.univie.ac.at/index.php/s/pUPFeZDXGtGNEpa)

#### Real data
- [BigCross](https://ucloud.univie.ac.at/index.php/s/ITlAQkZfIGFTvTD) see [[2]](https://doi.org/10.1145/2133803.2184450)
- [Activity recognition](http://archive.ics.uci.edu/ml/datasets/heterogeneity+activity+recognition) see [[1]](https://doi.org/10.1145/2809695.2809718) from UCI Data Repository
- [Higgs](https://archive.ics.uci.edu/ml/datasets/HIGGS) see [[3]](https://www.nature.com/articles/ncomms5308) from UCI Data Repository 
- [IoT botnet](https://archive.ics.uci.edu/ml/datasets/detection_of_IoT_botnet_attacks_N_BaIoT) see [[4]](http://wp.internetsociety.org/ndss/wp-content/uploads/sites/25/2018/02/ndss2018_03A-3_Mirsky_paper.pdf) from UCI Data Repository 

# Issues

Feel free to report [issues](https://gitlab.cs.univie.ac.at/martinp16cs/hilbertJoin/issues) about the code.

# Freqeuntly Asked Quesionts

[see FAQ](FAQ.md)

# Comparison partners

- [BLAS-join](https://gitlab.cs.univie.ac.at/Google-TPU/BLAS-join/)
- [Super-EGO](https://www.ics.uci.edu/~dvk/code/SuperEGO.html). 
  <br/> List of changes:
  - Changed from float to double precision.
  - Adapted to a self-join problem.
  - We take care of the dimensions, within our compilation process (D must be known at compile time).
  - Removed storing of results, and added thread independent counting.

# References

- [1] Allan Stisen, Henrik Blunck, Sourav Bhattacharya, Thor Siiger Prentow, Mikkel Baun Kjærgaard, Anind K. Dey, Tobias Sonne, Mads Møller Jensen:
Smart Devices are Different: Assessing and MitigatingMobile Sensing Heterogeneities for Activity Recognition. SenSys 2015: 127-140
- [2] Marcel R. Ackermann, Marcus Märtens, Christoph Raupach, Kamil Swierkot, Christiane Lammersen, Christian Sohler:
StreamKM++: A clustering algorithm for data streams. ACM Journal of Experimental Algorithmics 17(1) (2012)
- [3] P. Baldi, P. Sadowski, and D. Whiteson. 2014. Searching for exotic particles in high-energy physics with deep learning. Nature Commu- nications 5 (02 Jul 2014), 4308 EP –. Article.
- [4] Yisroel Mirsky, Tomer Doitshman, Yuval Elovici, and Asaf Shabtai. 2018. Kitsune: An Ensemble of Autoencoders for Online Network Intrusion Detection. In 25th Annual Network and Distributed System Security Symposium, NDSS 2018, San Diego, California, USA, February 18-21, 2018.

