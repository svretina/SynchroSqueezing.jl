# SynchroSqueezing.jl

`SynchroSqueezing.jl` is a planned pure-Julia package for synchrosqueezed
time-frequency transforms, frequency-ridge extraction, and oscillatory
component reconstruction.

The package is in its foundation phase. Its architecture will reuse mature
Julia implementations of continuous wavelet transforms and STFTs where they
are mathematically suitable. The package’s new numerical work will be the
synchrosqueezing, ridge-tracking, and reconstruction layers.

`ContinuousWavelets.jl` is the selected CWT backend. The package will add its
SST-specific adapter, derivative support, and metadata validation in the next
milestone.

## Status

Milestone 00 is in progress. The current source only establishes package,
testing, documentation, and reproducible-reference infrastructure; it does
not yet provide a public transform API.

See the [roadmap](roadmap/README.md), [numerical conventions](docs/src/conventions.md),
[dependency audit](docs/src/dependency_audit.md), and [provisional backend
decision](docs/src/backend_decision.md).

## Development

Run the Julia test suite with:

```sh
julia --project -e 'using Pkg; Pkg.test()'
```

Python is never a production dependency. The pinned reference environment is
described in [`reference/python/README.md`](reference/python/README.md).

## License

The project license has not been selected yet.
