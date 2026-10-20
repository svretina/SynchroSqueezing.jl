# Reconstruction audit

This audit predates the Milestone 03A correction and records why the original
Milestone 03 inverse could not be used for component extraction.

## Findings

1. **CWT normalization.** The original adapter used ContinuousWavelets.jl's
   `p=Inf` filters, adaptive scale spacing, changing scale-dependent wavelet
   width, and global `frameBound=1` normalization. This was not one fixed
   L¹-normalized mother wavelet sampled on a logarithmic dilation grid.
2. **Scale quadrature.** It used midpoint widths in log mean-frequency, with
   full endpoint weights. The corrected standard path uses fixed logarithmic
   inverse scales and trapezoidal endpoint weights.
3. **SST storage.** Reassignment already stored bin-integrated complex
   contributions. This part was correct and remains unchanged.
4. **Old synthesis correction.** The old inverse calculated
   ``R[m]=\sum_k q_k H_k[m]`` and divided the Fourier transform of the summed
   SST coefficients by `R[m]`. This is a discrete diagonal frame inverse.
5. **Retained data.** The diagonal inverse did not require input samples, but
   `mode=:full` also used a retained weighted-CWT residual and DC. The new
   standard oscillatory and masked inverses use only the SST matrix and scalar
   reconstruction metadata.
6. **Meaning of the reported machine precision.** The roughly `1e-12` tone
   errors primarily measured cancellation by `1/R[m]`; narrow output ranges
   additionally relied on excluded-coefficient residuals. They did not verify
   the direct scale-integral formula.
7. **Masked regions.** No arbitrary time-frequency mask API existed. Applying
   the global `1/R` filter to a selected band would not implement standard
   component synthesis.
8. **Signal classes.** An FFT negative-frequency energy test distinguishes
   real, positive-frequency analytic complex, and unsupported two-sided
   complex input.
9. **Analyticity.** ContinuousWavelets.jl stores analytic families on the
   positive rFFT grid. Its default Morlet configuration had measurable DC
   leakage, while the selected generalized Morse wavelet has exactly zero DC
   response for positive beta. Backend warnings at the upper scale describe
   Nyquist truncation, not negative-frequency leakage.
10. **Standard compatibility.** The backend can support standard WSST after a
    thin adapter: fixed log spacing, no global frame rescaling, one fixed Morse
    wavelet, and restoration of the analytical Fourier amplitude removed by
    the backend's per-filter peak normalization.

## Baseline and regression policy

Before this correction, 133 Julia tests passed; the periodic DC-plus-tone
frame-corrected reconstruction error was about `8.6e-13`, and the example
chirp ridge RMSE was about `0.072 Hz`. These remain useful implementation
regressions, but frame accuracy is now explicitly separate from standard SST
synthesis. Existing threshold, phase, grids, deposition, diagnostics, and
residual code was retained.

The production default is now the standard-compatible path. The former
behavior remains available with `standard=false` at planning time and
`method=:frame` during inversion.
