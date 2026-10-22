# Second-order wavelet synchrosqueezing

Second-order WSST sharpens coefficients from chirped modes by correcting the
first-order phase frequency with a local chirp-rate estimate. Select it with
the `order` keyword:

```julia
second = wsst(
    signal;
    fs = 1024.0,
    order = 2,
    cwt_options = (; voices_per_octave = 32),
    freqrange = (100.0, 450.0),
    nfrequencies = 700,
    threshold = RelativeThreshold(1e-5),
)
```

`order` denotes synchrosqueezing order; it is unrelated to generalized Morse
wavelet order.

## Operators and units

Let ``W(t,s)`` be the package's L1-normalized analytic CWT, ``W_t`` and
``W_{tt}`` its physical-time derivatives, and ``M`` the transform obtained by
differentiating the sampled Fourier daughter with respect to physical cyclic
frequency and dividing by ``i2\pi``. Thus ``M/W`` has units of seconds. The
complex first-order frequency is

```math
\widetilde f = \frac{1}{i2\pi}\frac{W_t}{W}.
```

The complex time-reassignment operator and its derivative are

```math
\widetilde\tau = t + \frac{M}{W}, \qquad
\partial_t\widetilde\tau = 1 + \frac{M_t}{W} -
    \frac{M W_t}{W^2}.
```

The local modulation operator is

```math
q = \frac{(i2\pi)^{-1}\left(W_{tt}/W-(W_t/W)^2\right)}
         {\partial_t\widetilde\tau},
```

in Hz/s. The corrected frequency deposited by WSST2 is

```math
\widehat f^{[2]} = \operatorname{Re}\left(
    \widetilde f - q\frac{M}{W}
\right).
```

These equations are the physical-unit form of Pham and Meignen's
time-derivative WSST2 formulation. All transforms are computed spectrally from
the same sampled Morse daughter bank. No numerical derivative of a ridge is
used.

With the package Fourier convention
``X(f)=\int x(t)e^{-i2\pi f t}\,dt``, physical-time differentiation uses
the multiplier ``i2\pi f``. For a sampled Morse daughter ``D_s(f)`` the
implemented physical-frequency derivative is

```math
\partial_f D_s(f) = \frac{D_s(f)}{f}
\left(\beta - \gamma x^\gamma\right),
```

where ``x`` is the backend's normalized positive frequency at scale ``s``.
The zero-frequency derivative is set to zero. The moment transform is then
``M=\mathcal{F}^{-1}\{X\,\partial_fD_s\}/(i2\pi)``; this factor and the
spectral time-derivative multipliers fix the signs and the seconds, Hz, and
Hz/s units in the equations above.

## Diagnostics and fallback

The transparent operator API exposes every scientific intermediate:

```julia
plan = plan_cwt(length(signal); fs = 1024.0, voices_per_octave = 32)
operators = second_order_operators(plan, signal)
estimate = second_order_frequency(
    operators;
    threshold = RelativeThreshold(1e-5),
    denominator_threshold = sqrt(eps(Float64)),
    fallback = :first_order,
    max_abs_chirp_rate = 2_000.0,
)
```

`estimate` contains corrected and first-order frequencies, complex time
reassignment, its derivative, chirp rate, a second-order validity mask, and a
fallback mask. Coefficients are second-order-valid only when all inputs are
finite, the CWT magnitude and backend reliability gates pass,
``abs(∂ₜτ̃)`` exceeds `denominator_threshold`, the frequency is positive, and
the optional chirp-rate bound passes.

With `fallback=:first_order`, an invalid second-order estimate uses the valid
first-order phase frequency. With `fallback=:exclude`, it is retained in the
optional reconstruction residual. Invalid values never become frequency-bin
indices.

Every `WSSTResult` records the order, estimator, implementation, denominator
threshold, fallback policy, chirp-rate limit, observed second-order-valid
fraction, fallback fraction, excluded fraction, and backend wavelet configuration in
`reassignment_metadata`.

The three fractions use all reliable, above-threshold coefficients with a
finite first-order frequency as their common denominator. When at least one
candidate exists, they sum to one for WSST2: a candidate is corrected,
assigned by fallback, or excluded.

## Reconstruction and limitations

WSST2 moves the original complex, quadrature-weighted CWT coefficients. It
does not replace them with frequency or chirp-rate values. The same standard
inverse, arbitrary masks, ridge reconstruction, analytic envelope, damping,
activity, and ownership APIs therefore apply to first- and second-order
results.

The initial implementation deliberately supports the standard analytic Morse
backend and the scalar reference deposition path. It is periodic, so boundary
errors must be reported separately. A sharper representation is not proof of
better component identification: compare frequency error, concentration,
reconstruction, damping, invalid fractions, and parameter sensitivity.

The source audit and pinned external commit are recorded in
`reference/wsst2/README.md`.
