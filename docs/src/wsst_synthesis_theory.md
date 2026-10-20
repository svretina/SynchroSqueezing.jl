# Standard WSST synthesis theory

The standard backend represents the analytic L¹ CWT

```math
W_x(a,b)=\frac1a\int x(t)\overline{\psi((t-b)/a)}\,dt.
```

ContinuousWavelets.jl parameterizes filters by inverse dilation `s=1/a`.
Hence `dlog(a)=-dlog(s)`, and the positive quadrature weights are the absolute
trapezoidal widths in `log(s)`.

For the real generalized Morse response used by the default adapter,

```math
\widehat\psi(u)=A(u/f)^{\beta}\exp(-(u/f)^\gamma),\qquad
A=2(e\gamma/\beta)^{\beta/\gamma},
```

the direct-synthesis admissibility constant is

```math
C_\psi=\int_0^\infty\frac{\widehat\psi(u)}u\,du
=\frac{A}{\gamma}\Gamma(\beta/\gamma).
```

`admissibility_constant` implements this expression and an independent
log-frequency quadrature. For the default `gamma=3`, `beta=20` wavelet both
give approximately `0.65534341497916`.

The scale sum is

```math
S(t_n)\approx\sum_k q_k W_k(t_n).
```

For real input, standard synthesis is `2real(S/Cψ)`. For input that is
already analytic complex, it is `S/Cψ`: applying another factor of two would
double the amplitude incorrectly. DC is outside the strictly analytic CWT and
is included only by the explicitly requested full-input mode.

Because the SST stores bin-integrated values, summing all frequency rows gives
the same `S(t)` as summing the accepted weighted CWT coefficients. Output-bin
widths are not applied again.

```julia
oscillatory = inverse_wsst(result) # method=:standard by default
component = reconstruct_region(result, mask)
full = inverse_wsst(result; mode=:full)
```

`reconstruct_region` accepts any Boolean frequency × time mask and requires
neither input samples nor retained CWT/derivative matrices.
