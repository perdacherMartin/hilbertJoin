### AVX-512 instructions

_What is the correct setting of parameters?_
There are two parameters, _KBLOCK_ and _stripes_. You can vary and play with them, but throughout our experiments we have always used _KBLOCK=4_ and _stripes=14_. These are the default settings for our approach. For uniform distributed data _KBLOCK=16_ is faster. 

_Why do you use AVX-512 instructions?_
If we apply our Hilbert-curve in Intel or GNU compilers, auto-vectorization will get eliminated. Writing code with AVX instructions simulates the behaviour of having an implemented auto-vectorized approach. Nevertheless, we belive that future compilers will profit from the locality assumptions of the Hilbert curve.
