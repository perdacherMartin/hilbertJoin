### AVX-512 instructions

_Why do you use AVX-512 instructions?_
If we apply our Hilbert-curve in Intel or GNU compilers, auto-vectorization will get eliminated. Applying AVX instructions simulates the behaviour of having auto-vectorized approach already implemented. Nevertheless, we belive that future compilers will profit from the locality assumptions of the Hilbert curve.
