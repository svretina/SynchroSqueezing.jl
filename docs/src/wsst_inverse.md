# Inverse WSST

Summing the bin-integrated matrix produces the discrete scale integral
``S(t)``. The default standard inverse divides by the scalar admissibility
constant `Cψ` and applies the signal convention:

- real input: twice the real part of `S/Cψ`;
- analytic complex input: retain the complex analytic reconstruction;
- two-sided complex input: reject explicitly.

DC is recovered from the separately stored mean only for `mode=:full` or
`:residual`. `mode=:oscillatory` uses the deposited SST matrix alone.

```julia
result = wsst(signal; fs=1024.0, preserve_residual=true)
xosc = inverse_wsst(result; mode=:oscillatory)
xres = inverse_wsst(result; mode=:residual)
xfull = inverse_wsst(result; mode=:full)
```

`:oscillatory` uses deposited bins only. `:residual` uses finite weighted CWT
coefficients excluded by reliability, threshold, invalid phase, or output
range, and adds DC. `:full` synthesizes deposited plus residual contributions
and adds DC. Full or residual inversion throws when the residual was not
retained or nonfinite coefficients made exact accounting impossible.

The residual is computed from excluded weighted CWT coefficients. It is never
formed by subtracting an SST reconstruction from the input. Thus a successful
full inverse still tests quadrature and synthesis.

`inverse_wsst(...; method=:frame)` retains the Milestone 03 diagonal finite-bank
correction for whole-result diagnostics. That method performs an FFT and
divides by the measured response. It is explicitly distinct from standard
WSST synthesis and cannot be requested for masked component reconstruction.

`reconstruct_region(result, mask)` performs standard synthesis from any
frequency × time subset. It never includes DC or an excluded-coefficient
residual and needs no retained CWT or input waveform.
