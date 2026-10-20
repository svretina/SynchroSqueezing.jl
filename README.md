# SynchroSqueezing.jl

`SynchroSqueezing.jl` is a pure-Julia package under development for
synchrosqueezed time-frequency transforms, frequency-ridge extraction, and
oscillatory component reconstruction.

Its architecture reuses mature Julia transforms where they are mathematically
suitable, while keeping synchrosqueezing, reconstruction, and future ridge
tracking as explicit package layers.

`ContinuousWavelets.jl` is the selected CWT backend. The package includes a
periodic adapter with derivative support, metadata validation, first-order
wavelet synchrosqueezing, residual accounting, and direct WSST inversion.

## Status

Milestone 03 is complete. The public API provides a periodic CWT, its exact
spectral time derivative, stable phase-frequency estimates, nearest and linear
first-order WSST kernels, and reconstruction for real and analytic complex
signals. Milestone 03A adds standard scalar-admissibility inversion and
arbitrary masked component reconstruction. Ridge extraction begins in
Milestone 04.

See the [roadmap](roadmap/README.md), [numerical conventions](docs/src/conventions.md),
[dependency audit](docs/src/dependency_audit.md), and [provisional backend
decision](docs/src/backend_decision.md), and [CWT adapter
contract](docs/src/backend_adapter.md).

## Development

Run the Julia test suite with:

```sh
julia --project -e 'using Pkg; Pkg.test()'
```

Python is never a production dependency. The pinned reference environment is
described in [`reference/python/README.md`](reference/python/README.md).

## License

The project license has not been selected yet.
