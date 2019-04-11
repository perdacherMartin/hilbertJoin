### AVX-512 instructions

_What is the correct setting of parameters for uniform distributed data?_
There are two parameters, _KBLOCK_ and _stripes_. Throughout our experiments we have always used _KBLOCK=4_ and _stripes=14_ which are the default settings for any kind of distribution. For uniform distributed data _KBLOCK=16_ is faster, but the variation of _stripes_ has no effect. 

_Why do you use AVX-512 instructions?_
If we apply our Hilbert-curve in Intel or GNU compilers, auto-vectorization will get eliminated. Writing code with AVX instructions simulates the behaviour of having an implemented auto-vectorized approach. Nevertheless, we belive that future compilers will profit from the locality assumptions of the Hilbert curve.

_What are the Parameters KBLOCK and STRIPES?_

- KBLOCK: check after KBLOCK dimensions, whether the $`\varepsilon`$-distance is already exeeded. 
- STRIPES: How many EGO-Stripes are used. See Section _"2.2. Determining the bounds"_ in our paper. 

_How to set KBLOCK and STRIPES?_

KBLOCK should be smaller then the dimension of the dataset. Within our distance calculation, we check after KBLOCK dimensions whether we have exceeded $`\varepsilon^2`$ or not. 

Best fitting values for active dimesnions are $`0,1,2,3,4,5`$ which corresponds to $`1,2,5,14,41 (=((3^j)+1)/2)`$ stripes, (for more details see paper Section 3.1 "Determination of the Bounds")

In our experiments (see paper) we _always_ use the following setting:
- KBLOCK=4
- active dimensions=3, which are exactly 14 stripes

For uniform data we suggest to use the following parameter settings:
- KBLOCK=16
- STRIPES=1

_What does the output mean?_

Here an example output. 

N;D;JPPP;THREADS;EPSILON;STRIPES;KBLOCK;TIME;ALGTIME;SORTTIME;INDEXTIME;REORDERTIME;COUNTS;LOADPERCENT;WH
200000;64;0.000000;64;0.20000000000000;14;4;0.794607;0.579982;0.130889;0.514304;0.083736;0;0.061758;0.000000

- N ... number of objects
- D ... dimensionality (number of features)
- JPPP ... join-partners per point _nSelectivity_ (see Section 4.1.3).
- THREADS ... number of threads used
- STRIPES ... bounds (Section 3.1 in paper)
- KBLOCK ... check after each _KBLOCK_ objects, whether we have exceeded epsilon distance. 
- TIME ... time spent for the total algorithm
- ALGTIME ... time spent for join
- SORTTIME ... time spent for sorting
- INDEXTIME ... time spent for determining bounds (Section 3.1)
- REORDERTIME ... time spent for reordering the dimensions (proposed by Super-EGO)
- COUNTS ... cardinaities
- LOADPERCENT ... load in percent
- WH ... energy in watthours (currently turned off)

