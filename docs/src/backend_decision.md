# CWT backend decision

## Decision

**Select `ContinuousWavelets.jl` v1.2.2 as the package's CWT backend.** It is
now a declared core dependency, constrained to the 1.2 release series. Its
derivative support, normalization metadata, and inversion compatibility must
still pass the next milestone’s adapter tests before an SST API is exposed.

`ContinuousWavelet.jl` v0.1.0 is rejected as a production backend and retained
only in the isolated experiment environment. Its planned-transform design is
interesting, but it is unregistered, must be Git-pinned, pulls `Plots.jl` and
graphics dependencies into its runtime, accepts only `Vector{Float64}` in its
documented transform method, and could not reconstruct the nominally matched
test tone accurately.

## Experiment environment

Reproduce the tests with:

```sh
julia --project=experiments -e 'using Pkg; Pkg.instantiate()'
julia --project=experiments experiments/continuouswavelets_baseline.jl
julia --project=experiments benchmarks/backend_comparison.jl
```

The environment locks `ContinuousWavelets.jl` at v1.2.2 and pins the
unregistered alternative to commit `23c49cdafee0ae21aeba2755924d6ef59b4e5637`.
It is intentionally separate from the package runtime environment.

## Baseline result: ContinuousWavelets.jl

Configuration: `Morlet(π)`, `NoAve()`, `β=2`, `fsample=1024 Hz`, `N=2048`.
The CWT returns a `2048 × 30` `ComplexF64` matrix in **time × scale** layout.
`getMeanFreq` returned an ascending 2.340215–261.948159 Hz grid, and the
scale-averaged peak for a 64 Hz test tone was 62.058680 Hz. Using
`PenroseDelta()` yielded relative full-reconstruction errors of `6.03e-3` for
the tone and `1.66e-2` for the linear chirp.

Complex input was accepted (`ComplexF64` output); Float32 input was accepted
(`ComplexF32` output). The selected configuration emits an analyticity warning
for its lowest-frequency wavelet. The SST adapter must expose a reliability or
admissibility policy that removes/marks such rows rather than treating every
returned coefficient as SST-safe.

## Nominal comparison result

For comparison only, the alternative used `MorletWavelet(1.0)` and scales
`fs ./ frequencies` from the selected backend. Its documentation identifies a
unit centre frequency, so this is a nominal coordinate match—not equality of
normalization, wavelet shape, boundary treatment, or inverse convention.

| N | CWT scales | ContinuousWavelets median CWT | ContinuousWavelet median transform | CW inverse error | Alternative synthesis error |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2,048 | 30 | 5.891 ms | 5.509 ms | 6.03e-3 | 9.54e-1 |
| 8,192 | 38 | 18.510 ms | 26.693 ms | 2.91e-3 | 9.54e-1 |
| 16,384 | 42 | 48.004 ms | 58.569 ms | 2.17e-3 | 9.54e-1 |
| 65,536 | 50 | 522.163 ms | 362.797 ms | 1.29e-3 | 9.54e-1 |

Each timing is a same-machine warm median of ten single evaluations and
includes the backend’s per-transform allocations. The returned transforms are
not mathematically equivalent enough for the last column to establish a
speed ranking. In particular, the alternate backend’s lower allocation count
does not offset its missing complex/Float32 transform support or failed
reconstruction under the tested setup.

At `N=2048`, `ContinuousWavelets.jl` allocated 17.49 MB across 1,449
allocations per CWT; the alternative allocated 5.27 MB across 301 allocations.
This identifies plan/workspace reuse as a later integration concern, not a
reason to duplicate CWT implementation now.

For the `N=2048` comparison, `ContinuousWavelets.jl` configuration construction
took 106.814 ms and its documented wavelet-kernel construction took 1.515 ms;
the alternative’s planned-kernel construction took 2.100 ms. These are not
equivalent initialization stages—the former backend makes the Fourier daughters
separate from configuration creation—so they are recorded as integration data,
not compared as a winner/loser result.

## Integration contract and constraints

The forthcoming adapter must:

1. Convert the documented time × scale coefficients to the package’s internal
   frequency/scale × time view without an unconditional full transpose.
2. Preserve `fsample`, wavelet parameters, `NoAve()` choice, scale order,
   physical frequency coordinates, inverse-frame selection, and boundary
   behavior in result metadata.
3. Reject or mark rows whose analytic-support diagnostic is unsuitable for
   phase reassignment.
4. Derive SST quadrature and inverse weights from this backend convention; do
   not reuse the backend’s `icwt` constants without verification.

## Open Milestone 02 questions

- Can public `computeWavelets` output and `cwt` overloads supply a spectral
  `∂tW` using the backend’s exact scales and normalization?
- Are reusable daughters and FFT plans available through stable public methods,
  and can they eliminate the observed per-transform allocation cost?
- Which inverse frame and restricted reliable-frequency band give the required
  WSST reconstruction behavior?
- What padding/cropping semantics and cone-of-influence equivalent does the
  backend expose?
