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

# Datasets used in our publication

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
- [Super-EGO](https://www.ics.uci.edu/~dvk/code/SuperEGO.html). <br/> Notes: we adapted the code to perform a self-join.
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

