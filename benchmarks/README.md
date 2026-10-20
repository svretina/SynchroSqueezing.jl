# WSST benchmarks

Use the isolated benchmark environment so BenchmarkTools is not a production
dependency:

```sh
julia --project=benchmarks -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=benchmarks benchmarks/wsst.jl
```

The script reports warmed runtimes and allocations for phase estimation,
nearest and linear WSST after CWT, the complete CWT+WSST pipeline, standard
whole-matrix synthesis, arbitrary-mask standard synthesis, and the separate
FFT-based finite-bank frame inverse.
