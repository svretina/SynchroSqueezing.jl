# Pham-Meignen WSST2 reference

Milestone 06 was audited against:

- D.-H. Pham and S. Meignen, *Second-order Synchrosqueezing Transform:
  The Wavelet Case, Comparisons and Applications*.
- Reference implementation: <https://github.com/phamduonghung/WSST2>
- Audited Git commit: `2721d07e943abcc6faf7bfcaf974b5dfad597882`
- Primary source file: `Wsst2_new.m`

The external MATLAB repository is a validation dependency only. It is not
vendored and is not used at package runtime. Reproduce the source audit with:

```sh
git clone https://github.com/phamduonghung/WSST2 /tmp/wsst2-reference
git -C /tmp/wsst2-reference checkout 2721d07e943abcc6faf7bfcaf974b5dfad597882
```

## Convention mapping

The MATLAB code computes `WT`, `omega`, `tau`, `phipp`, and `omega2`. Their
counterparts in this package are:

| MATLAB | SynchroSqueezing.jl | Units |
|---|---|---|
| `WT` | `operators.cwt.coefficients` | coefficient |
| `omega` | `estimate.first_order` | Hz |
| `tau` | `estimate.time_reassignment - t` | seconds |
| `phipp` | `estimate.chirp_rate` | Hz/s |
| `omega2` | `estimate.values` | Hz |

The reference uses normalized sample frequency and descending dilation
scales. This package uses physical seconds/Hz, ascending center frequency,
periodic boundaries, L1-normalized generalized Morse daughters, and
trapezoidal finite-bank log-scale quadrature. Consequently, raw matrices are
not compared without converting units, scale order, boundary handling, and
normalization.

The implemented estimator follows equations (21)-(27) of the paper. It uses
the reference time-derivative formulation, including the complex modulation
operator. It does not use `Wsst2_omega_old.m`, which separately takes real
parts before the correction and is retained upstream only as an older variant.

