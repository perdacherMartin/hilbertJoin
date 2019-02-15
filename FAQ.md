### AVX-512 instructions

_What is the correct setting of parameters for uniform distributed data?_
There are two parameters, _KBLOCK_ and _stripes_. You can vary and play with them, but throughout our experiments we have always used _KBLOCK=4_ and _stripes=14_ which are the default settings for any kind of distribution. For uniform distributed data _KBLOCK=16_ is faster, but the variation of _stripes_ has no effect. 

_Why do you use AVX-512 instructions?_
If we apply our Hilbert-curve in Intel or GNU compilers, auto-vectorization will get eliminated. Writing code with AVX instructions simulates the behaviour of having an implemented auto-vectorized approach. Nevertheless, we belive that future compilers will profit from the locality assumptions of the Hilbert curve.

_What are the Parameters KBLOCK and STRIPES?_

- KBLOCK: check after KBLOCK dimensions, whether the $$\varepsilon$$ distance is already exeeded. 
- STRIPES: How many EGO-Stripes are used. See Section 2.2 ``Determining the bounds''. 

_How to set KBLOCK and STRIPES?_

In our experiments (see paper) we _always_ use the following setting:
- KBLOCK=4
- STRIPES=14

For uniform data we suggest to use the following parameter settings:
- KBLOCK=16
- STRIPES=1