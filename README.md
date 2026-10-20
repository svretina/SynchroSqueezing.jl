# SynchroSqueezing.jl

`SynchroSqueezing.jl` is a planned pure-Julia package for synchrosqueezed
time-frequency transforms, frequency-ridge extraction, and oscillatory
component reconstruction.

The package is in its foundation phase. Its architecture will reuse mature
Julia implementations of continuous wavelet transforms and STFTs where they
are mathematically suitable. The package’s new numerical work will be the
synchrosqueezing, ridge-tracking, and reconstruction layers.

`ContinuousWavelets.jl` is the selected CWT backend. The package includes a
small periodic adapter with derivative support, metadata validation, and CWT
inversion; synchrosqueezing itself comes in later milestones.

## Status

Milestone 02 is complete. The current public API provides a periodic CWT,
its exact spectral time derivative, a phase-based frequency estimate, and
full CWT inversion. It does not yet provide an SST or ridge API.

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
