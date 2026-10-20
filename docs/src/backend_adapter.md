# CWT backend adapter

Milestone 02 exposes a small adapter around `ContinuousWavelets.jl`. It owns
the backend-specific Fourier details while the rest of the package sees one
layout and one set of physical units.

```julia
using SynchroSqueezing

fs = 1024.0
t = (0:1023) ./ fs
signal = cos.(2π * 64.0 .* t)

plan = plan_cwt(length(signal); fs)
result = cwt(plan, signal)
estimate = instantaneous_frequency(result; threshold = 1e-8)
recovered = inverse_cwt(result)
```

## Result contract

`result.coefficients` and `result.derivative` have shape
`frequency × time`; positive frequencies are strictly ascending and measured
in Hz. They are `PermutedDimsArray` views of the backend's native
`time × scale` matrices, so changing the orientation does not copy either
full matrix.

`result.frequency_weights` contains midpoint-bin widths in Hz.
`result.squeezing_weights` contains the positive log-frequency quadrature used
by WSST, and `result.synthesis_response` records the response of the weighted
sum of the actual backend daughters. Standard-compatible results also store
the independently computed direct-synthesis `admissibility` constant and the
measured 1%-error `reliable_frequency_range`.

`result.metadata` records the sample rate and interval, source and padded
lengths, crop range, boundary choice, backend normalization label, inverse
style, and exact backend configuration. The current adapter is periodic, so
the padded length equals the source length and `crop == 1:signal_length`. The
source mean is stored separately as `dc_component` and excluded from phase
reassignment.

`result.valid_rows` flags wavelets whose relative DC response is below the
requested `dc_tolerance` (default `0.01`). `result.reliability` expands that
row-level diagnostic over time. This is deliberately conservative: invalid
coefficients are not silently treated as valid SST inputs.

## Derivative and frequency convention

The adapter builds its derivative from the public
`ContinuousWavelets.computeWavelets` daughters used for the CWT itself. For a
periodic signal it evaluates

```math
\partial_t W = \operatorname{IFFT}\left(i\omega\,\widehat{x}_+\,\widehat{\psi}\right),
\qquad \omega = 2\pi f.
```

The FFT sign and scaling therefore agree with the selected backend's
one-sided transform convention. `instantaneous_frequency` returns
`imag(∂ₜW * conj(W)) / (2π * abs2(W))` in Hz, marks unreliable or
below-threshold entries invalid, and writes `NaN` rather than a misleading
zero for those entries.

The derivative experiment can be reproduced with:

```sh
julia --project experiments/cwt_derivative.jl
```

For the 64 Hz periodic tone, the test suite checks the peak phase estimate to
machine precision and checks the spectral derivative against centered
differences. The latter is only a validation baseline; it is not used in the
production calculation.

## Current limits

Only `boundary=:periodic` is supported. Reflective and zero-padded modes need
their own derivative and edge-reliability semantics before they can be exposed.
Plans are also precision-specific: a Float32 plan accepts Float32 or
ComplexF32 signals, and likewise for Float64. These explicit constraints keep
units, layout, and validity behavior predictable. The standard path currently
supports generalized Morse wavelets; other backend wavelets remain available
only through the explicitly selected legacy frame path.
