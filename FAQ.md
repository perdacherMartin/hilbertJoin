### AVX-512 instructions

_Why do you use AVX-512 instructions?_
If we apply our Hilbert-curve in Intel or GNU compilers, auto-vectorization will get eliminated. Writing code with AVX instructions simulates the behaviour of having an implemented auto-vectorized approach. Nevertheless, we belive that future compilers will profit from the locality assumptions of the Hilbert curve.
