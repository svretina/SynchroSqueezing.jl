# Dependency audit

Audit date: 2026-10-20. Julia’s local General registry was available during
this audit. Package behavior below is based on documented public interfaces and
must be rechecked against the exact resolved versions before a package becomes
a production dependency.

## Decision summary

No external numerical library is a hard runtime dependency in Milestone 00.
Milestone 01 selected `ContinuousWavelets.jl` v1.2.2 as the provisional CWT
backend, subject to the derivative and metadata checks documented in the
[backend decision](backend_decision.md). `ContinuousWavelet.jl` was evaluated
from a Git-pinned source because it is not registered, then rejected for the
production path. `DSP.jl` remains the provisional STFT backend for a later
FSST track. `Wavelets.jl` remains optional and is not a substitute for an
analytic continuous transform.

| Package | Publicly documented capability | CWT/SST relevance | Decision |
| --- | --- | --- | --- |
| [ContinuousWavelets.jl](https://github.com/UCD4IDS/ContinuousWavelets.jl) | 1-D CWT, Fourier-domain wavelets, `cwt`, `icwt`, `computeWavelets`, `getMeanFreq`, and multiple inverse-frame choices | Supports tested complex and Float32 inputs; returns time × scale matrices and exposes analyticity warnings | Selected core CWT dependency (`1.2`); adapter/derivative validation remains M02 work |
| [ContinuousWavelet.jl](https://hmatuschek.github.io/ContinuousWavelet.jl/stable/) | Planned CWT object, Morlet and Cauchy/Paul wavelets, transformed-array object, and synthesis | Source-level test confirms time × voice output, Float64-real-only transform dispatch, and poor nominal-scale reconstruction | Rejected for the production path; retained only as a pinned experiment comparison |
| [DSP.jl](https://docs.juliadsp.org/dev/periodograms/) | Complex STFT, windows, periodogram/spectrogram utilities | Appropriate FSST infrastructure, but standard STFT hop/centering/phase and inverse requirements need a dedicated adapter | Candidate for later optional/core FSST dependency; not used for WSST |
| [FFTW.jl](https://github.com/JuliaMath/FFTW.jl) | FFT interface, plans, FFTW controls | Needed only if a selected backend cannot publicly supply a spectral CWT derivative | Defer; use for a narrow derivative bridge only after public/backend-extension options are exhausted |
| [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl) | DWT/WPT/MODWT, thresholding, denoising, noise estimation, GPU array support | Discrete transforms and preprocessing are not the analytic CWT required for standard WSST | Optional preprocessing extension only |
| [WaveletsExt.jl](https://github.com/UCD4IDS/WaveletsExt.jl) | Companion extensions/utilities to the Julia wavelet ecosystem | No verified public primitive required by the WSST core | Inspect only if the CWT prototype identifies a concrete missing facility |
| [NPZ.jl](https://github.com/JuliaIO/NPZ.jl) | NumPy `.npy`/`.npz` interchange | Allows Julia tests to consume small Python fixtures without Python at runtime | Test-only candidate |
| [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) | Repeatable Julia microbenchmarking | Required for warm, allocation-aware performance evidence | Development-only candidate |
| [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl) | Julia documentation builds | Required to publish conventions, theory, validation, and audit decisions | Documentation-only dependency |

## Backend findings to verify experimentally

### ContinuousWavelets.jl

Its README documents Fourier-domain storage of wavelets and transforms computed
using an FFT/rFFT, pointwise multiplication, and an inverse transform. It
documents `cwt`, `icwt`, `computeWavelets`, `getMeanFreq`, analytic and real
families, and inverse-frame options including `NaiveDelta`, `PenroseDelta`, and
`DualFrames`. A shown CWT result has shape samples × scales, so an adapter is
required for the package’s scale/frequency × time contract. The documentation
also warns that selected wavelets can have undesirable low/high-frequency
support, making an explicit SST-suitability check mandatory.

Open questions for Milestones 01–02:

1. Exact scale/frequency units when a physical sample rate is supplied.
2. Real and complex input behavior, father-wavelet contribution, and boundary
   semantics.
3. Whether public APIs expose frequency-domain kernels sufficient to calculate
   `∂tW` with the backend’s exact normalization.
4. Allocation behavior, plans/workspaces, and Float32 behavior.
5. Whether `icwt` accuracy and its selected inverse frame are compatible with
   a derived inverse WSST.

### ContinuousWavelet.jl

Its documentation describes `ContinuousWaveletTransform` as a reusable planned
transform with precomputed scale kernels, a potentially attractive design for
repeated waveform analyses. It documents Morlet and Cauchy wavelets,
time × voice coefficient storage, and `synthesis` for reconstruction. Its
documented `transform` signature accepts `AbstractArray{Float64,1}`; complex
input support and precision generality therefore require direct verification.
The documentation available for this audit is substantially older than the
current `ContinuousWavelets.jl` release information, so maintenance and Julia
compatibility are explicit selection risks.

### DSP.jl

`DSP.Periodograms.stft` returns complex DFT coefficients in a frequency × time
matrix; the documentation states that each column is one windowed segment. This
is suitable raw material for FSST, unlike power-only spectrogram output. It does
not by itself establish the derivative-window, moment-transform, phase,
centering, hop-size, or synthesis conventions needed for an invertible FSST.
Those form a separate later milestone.

## Missing algorithms and reuse position

This audit found no verified maintained Julia package whose public interface
directly supplies all of the needed phase reassignment, conservative complex
deposition, SST inversion, or constrained continuous-ridge dynamic programming.
Milestone 04 will conduct a narrower solver audit before implementing the
specialized Viterbi-style recurrence. Avoid adopting a broad optimization
framework unless it demonstrably provides the desired cost, constraints, and
validated complexity.

## Dependency policy

- Production runtime: Julia plus only dependencies that survive the backend
  selection and expose the metadata needed for correct SST.
- Extensions: optional preprocessing and GPU/plotting facilities.
- Development/test: `BenchmarkTools`, `NPZ`, and generated Python fixtures.
- Documentation: `Documenter` only.
- Python: external comparison environment only; never imported by the package.

## Sources

- [ContinuousWavelets.jl README and API examples](https://github.com/UCD4IDS/ContinuousWavelets.jl)
- [ContinuousWavelet.jl transform documentation](https://hmatuschek.github.io/ContinuousWavelet.jl/stable/transform/)
- [ContinuousWavelet.jl synthesis documentation](https://hmatuschek.github.io/ContinuousWavelet.jl/stable/synthesis/)
- [DSP.jl STFT documentation](https://docs.juliadsp.org/dev/periodograms/)
- [Wavelets.jl README](https://github.com/JuliaDSP/Wavelets.jl)
