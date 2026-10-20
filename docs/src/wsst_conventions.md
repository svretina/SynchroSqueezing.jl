# WSST conventions

## Backend contract

WSST reuses `ContinuousWavelets.jl`; it does not implement another CWT. The
adapter records the generated Fourier daughters, dimensionless scale
coordinates, physical mean frequencies, sample rate, boundary mode, inverse
style, and source length. Public matrices are frequency/scale × time, while
the backend's time × scale CWT is exposed through a non-copying view.

The standard path requests a fixed logarithmic generalized-Morse bank. The
backend initially normalizes each `p=Inf` Fourier daughter by its sampled peak;
the adapter restores the analytical Morse amplitude and disables global frame
rescaling. The recorded normalization is `:analytic_l1`. Backend inverse scale
`s` is proportional to physical center frequency and equals reciprocal
dilation up to the fixed wavelet-frequency factor.

## Derivative and phase

For periodic input, the adapter evaluates

```math
D_k = \operatorname{IFFT}(i\omega\widehat x_+\widehat\psi_k),
\qquad \omega=2\pi f,
```

using the same daughters as the CWT. Thus `D_k` differentiates with respect to
seconds. Tests at two sample rates verify the units and Fourier sign.

## Quadrature and synthesis response

The scale quadrature is constructed in the actual log center-frequency
coordinate. Endpoint weights use adjacent spacing and interior weights use
midpoint widths:

```math
q_1=\tfrac12\Delta\log f_1,\quad
q_k=\tfrac12(\Delta\log f_{k-1}+\Delta\log f_k),\quad
q_K=\tfrac12\Delta\log f_{K-1}.
```

These positive weights give stable squeezed coefficients and use the
trapezoidal endpoint rule. The adapter also measures the finite discrete
synthesis response from the actual filters,

```math
H[m]=\sum_k q_k\widehat\psi_k[m].
```

This response is stored for validation and the explicitly selected
`method=:frame` alternative. Standard inversion uses the independently
calculated scalar admissibility constant and never divides by this response.

## Signals, DC, and boundaries

Real signals and positive-frequency analytic complex signals are supported.
Two-sided complex input is rejected because a one-sided CWT cannot recover its
negative-frequency content. The sample mean is removed before analysis and
stored as a separate DC channel; it is never assigned a zero-frequency ridge.

Only `boundary=:periodic` is supported. There is therefore no finite-support
edge region: the row reliability mask records daughter DC leakage,
independently of phase validity. Zero and reflective padding remain
unsupported until their time-dependent reliability semantics are defined.
