# WSST backend validation

## Decision

Three options were evaluated in the required order.

1. The unmodified ContinuousWavelets.jl defaults were rejected for standard
   synthesis because adaptive spacing, width variation, frame rescaling, and
   filter-wise peak normalization prevent one scalar admissibility identity.
2. A thin adapter succeeds and is selected. It asks the backend for a fixed
   logarithmic Morse bank, disables global frame normalization, and restores
   the analytical Fourier amplitudes. FFT analysis and derivative machinery
   remain backend-owned.
3. A dedicated FFTW CWT is therefore unnecessary and was not implemented.

The default wavelet is generalized Morse with `gamma=3`, `beta=20`, and unit
peak-frequency parameter. Its response is zero at DC and supported only on
the positive Fourier grid. The nominal physical center is obtained directly
from the analytical peak rather than an unreliable discrete mean at the first
few FFT bins.

## Finite response

Validation computes, without correcting it,

```math
R[m]=\sum_k q_kH_k[m]
```

and compares `R/Cψ` with one. For `N=1024`, `fs=1024 Hz`, and 64 Hz, the
absolute response errors are approximately:

| Voices/octave | `abs(R/Cψ - 1)` |
| ---: | ---: |
| 4 | `4.42e-4` |
| 8 | `1.51e-10` |
| 16 | `1.44e-15` |

With 16 voices/octave, the contiguous region satisfying a 1% response-error
criterion is approximately 2–334 Hz. The exact upper limit depends on signal
length and sample rate; every plan stores its measured
`reliable_frequency_range`. Accuracy degrades near DC and Nyquist because the
finite scale bank truncates the continuous integral.

At fixed `fs=1024 Hz`, increasing the periodic record from 128 through 1024
samples extends the measured lower reliable boundary from 16 to 2 Hz (16, 8,
4, and 2 Hz respectively). This is the expected added low-frequency scale
coverage, rather than an improvement manufactured by response correction.

The alternative `method=:frame` divides by `R[m]` and remains useful for
whole-signal diagnostics. It is deliberately unavailable in
`reconstruct_region`, since filtering a restricted region by the global frame
response is not the standard component formula.

## ssqueezepy comparison

The pinned ssqueezepy 0.6.6 generalized-Morse L¹ response with the same
`gamma=3`, `beta=20` parameters reports `Cψ=0.6553435757496617`; Julia's
analytical value is `0.6553434149791622` (relative difference about `2.5e-7`,
from ssqueezepy's numerical integration). The reference generator records raw
CWT, derivative, phase, SST, standard inverse, and known-band components.
Reflect padding and scale-bound selection differ, so boundary samples and raw
matrix extents are compared through analytical physical quantities rather
than asserted elementwise equal.

The measured values are also stored in the machine-readable repository report
`validation/wsst_03a.json`.
