# WSST component reconstruction

Component reconstruction uses a known time-dependent frequency region without
performing ridge tracking:

```julia
curve = 40 .+ 10 .* time
mask = abs.(reshape(result.frequencies, :, 1) .- reshape(curve, 1, :)) .<= 4
component = reconstruct_region(result, mask)
```

The mask selects only stored complex SST bins. Synthesis sums those bins and
applies the same scalar admissibility factor as the full oscillatory inverse.
It does not access the original signal, CWT, derivative, excluded-coefficient
residual, or frame-response filter.

Deterministic tests cover widely separated real and analytic-complex tones,
three unequal-amplitude tones, parallel linear chirps, a modulated quadratic
chirp paired with an exponential chirp, and mask-width sensitivity. The easy
two-tone case reconstructs both components below 1% relative L² error and
checks amplitude, phase, and cross-component contamination independently.

The reconstruction of all unselected SST bins is also calculated. This
SST-domain complement is distinct from both the retained weighted-CWT residual
and the input-space difference between the original signal and reconstructed
components.

Components whose instantaneous frequencies approach within a wavelet's
resolution violate the separation assumptions. The validation suite keeps
such a case as an explicit limitation: output remains finite, but component
error is expected to become large. Crossing-aware ownership remains a later
ridge-tracking problem.
