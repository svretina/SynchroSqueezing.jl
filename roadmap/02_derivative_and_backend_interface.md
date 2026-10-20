# Milestone 02 — Derivative support and backend interface

## Objective

Make the chosen CWT backend usable by WSST without leaking its storage conventions into the public API. The central technical question is how to obtain the exact physical-time derivative of the CWT.

## Deliverables

- `experiments/cwt_derivative.jl` and a written numerical recommendation.
- `src/backends/interface.jl` (or equivalent) defining an internal, multiple-dispatch interface for:
  - complex coefficients;
  - physical scales/frequencies and coefficient quadrature weights;
  - time derivative of the transform;
  - reconstruction normalization and full inverse;
  - padding, crop, and edge-reliability metadata.
- An adapter for the selected backend that exposes a consistent internal frequency-or-scale × time view without copying/transposing complete matrices unless profiling demonstrates that it is necessary.
- SST-wavelet validation, including analyticity, frequency coverage, and inverse compatibility.

## Derivative investigation

Prefer, in this order:

1. A documented public backend facility that can produce the spectral derivative.
2. A small upstream enhancement that exposes derivative-CWT or in-place support.
3. A narrowly scoped FFTW implementation of `IFFT(iω * FFT(x) * wavelet_response)` using the backend’s exact wavelet kernels, scale grid, and normalization.

Compare the selected derivative with analytical pure-tone expectations and centered finite differences. A finite-difference derivative is a validation baseline, not the production default. Do not use undocumented backend internals unless they are isolated and version-tested.

## Acceptance gate

The derivative is correct on analytical interior regions, units are explicitly physical Hz after division by `2π`, all adapter metadata required for inversion is available, and tests show that layouts and valid masks are unambiguous.

## Explicit non-goals

- No reimplementation of a complete CWT backend.
- No SST2, FSST, or performance threading yet.
