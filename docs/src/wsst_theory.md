# WSST theory

The first-order wavelet synchrosqueezing transform (WSST) sharpens an analytic
CWT by moving each complex coefficient from its analysis-wavelet center to the
coefficient's local phase frequency. It changes frequency, not time.

For the package CWT coefficient ``W_k(t)`` and its physical-time derivative
``D_k(t)``, the estimated cyclic frequency is

```math
\widehat f_k(t) =
\frac{\operatorname{Im}(D_k(t)\overline{W_k(t)})}
     {2\pi |W_k(t)|^2}.
```

The derivative is in units per second, so the result is in Hz. Undefined or
numerically rejected estimates remain `NaN` and have a separate validity mask.

## Bin-integrated coefficients

Let ``q_k`` be the positive scale quadrature weight. Nearest reassignment adds
``q_k W_k(t)`` to the closest output center. Linear reassignment splits that
same complex number between adjacent centers using their physical-frequency
distance. Therefore, for every time column,

```math
\sum_j T_j(t) = \sum_{k\in\mathcal A(t)} q_k W_k(t)
```

up to floating-point accumulation error. `T` stores bin-integrated complex
coefficients, not a density and not power. Inversion sums bins directly; it
does not multiply by output-bin widths. This identity is not energy
conservation because complex values can interfere after deposition.

`freqrange=(fmin, fmax)` gives the first and last output-bin centers. Explicit
centers and edges are also supported. Out-of-range estimates are retained in
the optional reconstruction residual rather than clamped into an edge bin.

Milestone 03A adds reconstruction from any explicitly supplied frequency–time
mask. Automatic ridge extraction and construction of those masks remain part
of Milestone 04.
