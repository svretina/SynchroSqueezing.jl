# WSST numerics

## Thresholds and validity

`NoThreshold()` accepts every nonzero finite coefficient.
`AbsoluteThreshold(a)` uses coefficient magnitude `abs(W) > a`.
`RelativeThreshold(r)` uses `r * maximum(abs, W)` over finite values in the
complete transform. Relative thresholding currently supports only
`reference=:global`; unsupported references fail at construction.

The phase ratio rescales each coefficient/derivative pair before forming its
numerator and denominator. This avoids avoidable underflow and overflow for
subnormal and very large coefficients. A point is phase-valid only when the
coefficient, derivative, denominator, and result are finite and its magnitude
passes the threshold. Negative and zero frequencies can be valid estimates,
but a positive output grid classifies them as out of range.

## Grids and deposition

Linear, logarithmic, and arbitrary increasing positive centers are supported.
Generated ranges denote centers. Edges are arithmetic midpoints for linear or
arbitrary grids and geometric midpoints for log grids. Explicit edges must be
finite, positive, increasing, and contain their corresponding centers.

Nearest reassignment selects the physically closest center. Linear
reassignment preserves both the complex contribution and its first frequency
moment between neighboring centers. At an outer half-bin it deposits into the
edge center; beyond the outer edge it records an exclusion without clamping.

The optimized kernel evaluates phase and deposits in one pass, without a full
phase matrix. The separately callable reference path first constructs the
phase matrix. Randomized tests require exact agreement of their deposited
matrices under the deterministic scalar execution model.

## Precision and storage

Plans are precision-specific. Float32 input produces ComplexF32 coefficients;
Float64 input produces ComplexF64. `keep_cwt`, `keep_derivative`, and
`keep_phase` default to false. The reconstruction residual is optional through
`preserve_residual`; disabling it reduces storage but disables full inversion.

The periodic algorithm has no cone-of-influence crop. Tests near DC and
Nyquist are robustness checks, not claims that edge-of-band resolution matches
a well-resolved interior tone.
