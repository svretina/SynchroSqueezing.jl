# Milestone 00 — Foundations and dependency audit

## Objective

Create the package, documentation, tests, and reference environment while deciding which existing Julia libraries can safely provide CWT and STFT infrastructure. Do not implement a custom CWT in this milestone.

## Deliverables

- Package skeleton, project metadata, CI, documentation scaffolding, deterministic test infrastructure, and an optional Python reference environment.
- `docs/src/conventions.md` defining FFT sign convention, angular versus cyclic frequency, sample interval, coefficient layout, positive-frequency ordering, padding, normalization, and edge-reliability conventions.
- `docs/src/dependency_audit.md` covering `ContinuousWavelets.jl`, `ContinuousWavelet.jl`, `Wavelets.jl`, `DSP.jl`, `FFTW.jl`, `WaveletsExt.jl`, and any maintained Julia reassignment, ridge-tracking, or quadratic-distance-transform package found during the audit.
- A pinned `ssqueezepy` reference environment; it must not become a Julia runtime dependency.

## Dependency classification

| Package | Intended role | Initial status |
| --- | --- | --- |
| `ContinuousWavelets.jl` | Primary CWT/inverse-CWT candidate | Provisional backend candidate |
| `ContinuousWavelet.jl` | Planned CWT comparison candidate | Evaluation-only until selected |
| `DSP.jl` | STFT, windows, signal utilities for eventual FSST | Candidate dependency/extension |
| `FFTW.jl` | FFT planning and possible derivative support | Core candidate |
| `Wavelets.jl` | Discrete-wavelet denoising/preprocessing | Optional only |
| `NPZ.jl` | Python fixture interchange | Test-only candidate |

For every audited package, record its public API, registry/compatibility status, maintenance/testing evidence, complex-signal support, coefficient layout and normalization, frequency/scale accessibility, derivative support, inverse quality, and plans/workspaces. State whether it is required, optional, development-only, or rejected.

## Acceptance gate

The package loads, CI and deterministic test plumbing run, conventions are documented, and the audit gives a defensible shortlist. No backend is selected based solely on speed.

## Explicit non-goals

- No production CWT, wavelet family, general STFT, or custom FFT implementation.
- No WSST, ridges, SST2, GPU path, or performance claim.
