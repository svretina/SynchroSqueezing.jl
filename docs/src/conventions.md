# Numerical conventions

These are the package-level conventions to be preserved by every selected
backend and every SST algorithm. A backend that cannot expose enough metadata
to satisfy them is not suitable for the production path.

## Sampling and frequency units

- A uniformly sampled signal has sample rate `fs` in Hz and `dt = 1 / fs` in
  seconds.
- Test grids use `t = (0:N-1) / fs`; no sample is duplicated at the endpoint.
- Cyclic frequency is denoted `f` and measured in Hz. Angular frequency is
  `ω = 2πf` and measured in rad/s. Dimensionless wavelet scale is never
  presented as Hz without an explicit backend-specific conversion.
- Positive-frequency analyses require `0 < fmin < fmax < fs / 2`. DC and the
  Nyquist bin are represented explicitly rather than inferred from row number.

## Fourier and derivative convention

The production transform convention is recorded with each backend adapter.
For a continuous-time convention using `exp(-iωt)` in the forward transform,
the physical-time derivative is evaluated spectrally as multiplication by
`iω` before the inverse transform. The first-order reassignment frequency is

```math
\widehat f(a,t) = \frac{1}{2π}\operatorname{Im}\left(
    \frac{\partial_t W(a,t)}{W(a,t)}
\right).
```

The implementation must verify its backend’s FFT sign and scaling before using
this identity. Invalid low-magnitude coefficients are tracked by a mask; they
are never interpreted as zero-frequency estimates.

## Layout and ordering

The public result contract will use `coefficients[frequency_or_scale, time]`.
Positive-frequency coordinates are ascending. Backends with time × scale
storage are adapted at the boundary; numerical kernels must not silently
reorder data, and full-matrix copies/transposes require a measured rationale.

## Normalization, weights, and inversion

Every result stores the analysis-wavelet identity and parameters, transform
normalization, physical grid, bin edges, reassignment quadrature weights,
padding/crop metadata, and reconstruction normalization. The inverse SST
weight is derived from the selected backend convention; it is never fitted to
match a reference output.

## Boundaries and reliability

Padding mode, original and padded lengths, crop indices, and any
cone-of-influence/reliability mask are metadata. Tests report boundary errors
separately from reliable-interior errors.

## Precision and determinism

Float32 and Float64 paths are tested independently. A single-thread execution
mode is the numerical regression reference. Parallel paths must preserve the
documented result within established floating-point tolerances.
