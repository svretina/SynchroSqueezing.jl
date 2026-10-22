# SynchroSqueezing.jl

`SynchroSqueezing.jl` is a pure-Julia package under development for
synchrosqueezed time-frequency transforms, frequency-ridge extraction, and
oscillatory component reconstruction.

Its architecture reuses mature Julia transforms where they are mathematically
suitable, while keeping synchrosqueezing, ridge tracking, and reconstruction
as explicit package layers.

`ContinuousWavelets.jl` is the selected CWT backend. The package includes a
periodic adapter with derivative support, metadata validation, first-order
wavelet synchrosqueezing, residual accounting, and direct WSST inversion.

## Status

Milestones 03 through 04C are complete. The public API provides a periodic CWT, its exact
spectral time derivative, stable phase-frequency estimates, nearest and linear
first-order WSST kernels, and reconstruction for real and analytic complex
signals. Milestone 03A adds standard scalar-admissibility inversion and
arbitrary masked component reconstruction, frequency-constrained global ridge
tracking, and conservative multi-component ownership.
Milestone 04B adds exact constant-velocity second-order dynamic programming
and an explicitly heuristic prediction-window-pruned strategy.
Milestone 04C adds analytic ridge-component envelopes, weighted log-linear
and nonlinear exponential damping fits, local damping rates, explicit sample
selection, local WSST-prominence diagnostics, persistent operational mode
termination, and separate fit-versus-reconstruction sensitivity diagnostics.
Milestone 05A is in progress: the repository now includes a frozen
full-pipeline scientific regression, reproducible stage-by-stage profiling,
and a memory-bounded ridge-DP optimization. On the recorded 1024-sample
workload, the first optimized ridge kernel was about 26.6 times faster while
preserving the dense reference path exactly.

See the [roadmap](roadmap/README.md), [numerical conventions](docs/src/conventions.md),
[dependency audit](docs/src/dependency_audit.md), and [provisional backend
decision](docs/src/backend_decision.md), and [CWT adapter
contract](docs/src/backend_adapter.md).
Performance methodology and current measurements are documented in
[`docs/src/performance.md`](docs/src/performance.md).

## Development

Run the Julia test suite with:

```sh
julia --project -e 'using Pkg; Pkg.test()'
```

Python is never a production dependency. The pinned reference environment is
described in [`reference/python/README.md`](reference/python/README.md).

## License

The project license has not been selected yet.
