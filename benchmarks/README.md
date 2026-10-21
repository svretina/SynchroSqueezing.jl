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

Frequency-constrained ridge solvers can be benchmarked independently with:

```sh
julia --project=benchmarks benchmarks/ridges.jl
```

The optional first argument sets the number of frequency bins. The script
reports maximum-energy, dense-DP, and physically banded-DP costs separately.

Compare first-order, exact second-order, and prediction-window-pruned costs
with:

```sh
julia --project=benchmarks benchmarks/predictive_ridges.jl
```

Optional arguments set the number of frequency bins and time samples. Exact
second-order storage grows quadratically with the number of bins inside the
search range; use representative constrained ranges when benchmarking.
