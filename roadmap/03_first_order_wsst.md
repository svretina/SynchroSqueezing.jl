# Milestone 03 — First-order WSST

## Objective

Implement the new SST-specific machinery on top of the selected backend: stable instantaneous-frequency estimation, complex-coefficient reassignment, and mathematically consistent inversion.

## Scope

- Phase estimator based on `Im(dW * conj(W)) / (2π * abs2(W))` above a documented validity threshold.
- Absolute and relative magnitude thresholds plus a validity mask; invalid estimates are never silently converted to zero frequency.
- Ascending physical-frequency grid with explicit bin edges.
- Nearest-bin reassignment and conservative linear-bin deposition. Both accumulate complex coefficients and account for out-of-range/nonfinite contributions in optional diagnostics.
- Quadrature/reconstruction weights derived from the backend’s actual CWT normalization.
- Full inverse WSST and a minimal public result type retaining only necessary metadata; optional retention of CWT, derivative, and phase diagnostics.

## Validation

Build deterministic analytical tests for zero input, complex and real tones, linear/quadratic/exponential chirps, amplitude modulation, noise, odd/even lengths, and near-DC/Nyquist behavior. Verify CWT derivative behavior, tone-frequency recovery, coefficient-conservation identities, Float32/Float64 paths, and interior-region full reconstruction. Start with scalar reference kernels; fuse phase and deposition only after their independent results agree.

Compare with a pinned `ssqueezepy` release in stages: wavelet response, scales/grid, raw CWT, derivative, instantaneous frequency, SST matrix, and inverse. Supply identical scales, thresholds, wavelet parameters, padding, and frequency bins wherever possible. Agreement with Python is evidence, not the sole correctness criterion.

## Acceptance gate

First-order WSST localizes isolated tones at the correct frequency, conserves accepted weighted coefficients within floating-point error, reconstructs well-resolved noiseless Float64 signals to a justified interior error target (initially below `1e-3`), and has reproducible Python comparison fixtures.

## Explicit non-goals

- Ridge extraction, multi-component ownership, GPU execution, FSST, and SST2.
