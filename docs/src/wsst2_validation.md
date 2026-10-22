# WSST2 validation

Milestone 06 validates the operator before treating sharper plots as evidence.
The deterministic suite covers Float32 and Float64, positive and negative
linear chirp rates, a constant tone, a fast real chirp, zero input,
ill-conditioned time reassignment, explicit first-order fallback, chirp-rate
limits, coefficient conservation, inversion, ridge extraction, and analytic
component features.

For the paper's linear-FM analytical case, the real part of the local
modulation operator equals the chirp rate and the corrected frequency is
exact. The numerical tests evaluate those quantities at the
dominant CWT scale in the reliable interior. Float64 errors are near machine
precision; Float32 uses a precision-appropriate tolerance.

The quantitative comparison uses identical CWT, threshold, output grid, and
deposition settings for WSST1 and WSST2. On the recorded 512 Hz, 1024-sample,
80 Hz/s chirp, WSST2 reduces the median energy-weighted frequency width while
maintaining sub-bin dominant-frequency error. The zero-chirp case is also
tested so second order is not credited for an improvement where none is
expected.

`validation/wsst2_06.json` records the exact configuration, interior and
boundary frequency errors, concentration widths, chirp-rate error, validity
fractions, reconstruction errors, and runtime/allocation measurements. The
recorded 1024-sample comparison is:

| Chirp rate (Hz/s) | SST1 interior RMSE (Hz) | SST2 interior RMSE (Hz) | SST1 width (Hz) | SST2 width (Hz) |
|---:|---:|---:|---:|---:|
| 0 | 0.0945 | 0.0945 | 0.0945 | 0.0945 |
| 20 | 0.0645 | 0.0642 | 0.0578 | 0.0556 |
| 80 | 0.0830 | 0.0642 | 0.0780 | 0.0556 |

The zero-rate case is unchanged, the 20 Hz/s improvement is negligible at
this grid resolution, and the 80 Hz/s case is the first recorded regime with
a clear reduction in both ridge-frequency error and energy width. This is a
numerical regime statement, not a universal crossover rate. Boundary errors
remain much larger: at 80 Hz/s their RMSE is 9.39 Hz for SST1 and 5.31 Hz for
SST2. The periodic boundary must therefore remain outside scientific fits.

At 80 Hz/s, 81.9% of threshold-eligible coefficients receive a second-order
estimate and 18.1% use the configured first-order fallback. Full finite-frame
reconstruction errors are below ``6\times10^{-16}`` for both orders because
reassignment moves the original weighted coefficients rather than replacing
them.

The BNS-like example in `examples/wsst2_bns.jl` is a controlled synthetic
stress test with two damped, nonlinearly evolving post-merger-like modes. Its
frequencies and damping constants are declared inputs, not inferred universal
BNS relations. It reports ridge, reconstruction, damping, support, boundary,
and sensitivity metrics for both WSST orders in `validation/wsst2_bns.json`.

For the rapidly rising mode, the interior ridge RMSE changes from 0.546 Hz to
0.492 Hz; for the second mode it changes from 3.453 Hz to 2.251 Hz. The second
mode's damping error improves from -3.71% to -1.80%, while the first mode's
damping error worsens from -4.70% to -6.82%. Component reconstruction also
worsens slightly for both modes. This is exactly why frequency concentration,
waveform reconstruction, and damping are reported as separate outcomes.

The Soultanis example in `examples/soultanis_waveforms.jl` additionally applies
the implementation to the eight Z4c MPA1 equal-mass numerical-relativity
strains distributed with a sibling `post-merger-bilby` checkout. The waveform
files are not vendored here: their directory may be supplied explicitly or via
`SOULTANIS_DIR`. Each strain is aligned at peak amplitude, resampled to 16,384
Hz, and restricted to -2–30 ms around merger. The resulting CWT, WSST1, and
WSST2 comparison is stored in
`examples/soultanis_cwt_wsst_comparison.png`; preprocessing and reassignment
fractions are recorded in `validation/soultanis_wsst2.json`.

The associated Fourier comparison uses the post-merger `h+` strain, a Tukey
window with ``\alpha=0.1``, and zero-padding to 16,384 samples. Padding reduces
the plotted FFT bin spacing to 1 Hz but does not improve the approximately
33 Hz physical resolution of a 30 ms record. The transform comparison uses the
Tukey window to suppress periodic-boundary artifacts but deliberately does not
pad to 16,384 samples; doing so would mostly increase CWT/WSST memory and
runtime. This is an application example, not a claim about physical mode
classification or catalogue-wide validation.

The separate `examples/soultanis_figure1.jl` workflow reproduces the visual
conventions of Figure 1 in arXiv:2111.08353v2 before comparing CWT, WSST1, and
WSST2: the 1.25 + 1.25 solar-mass reference case, 0–15 ms after merger,
1.35–4.12 kHz, linear per-panel coefficient normalization, and the `turbo`
color map. The mass-dependent white tracks come directly from the paper's
Eqs. (10)–(12). The paper does not specify its wavelet family and parameters,
scale grid, boundary treatment, or exact coefficient normalization, so the
example explicitly records its analytic Morse substitute and does not claim
an exact numerical reproduction. See `validation/soultanis_figure1.json` for
the complete separation between stated, figure-read, unavailable, and chosen
parameters.

Two further Soultanis experiments probe weak lower-frequency structure. The
focused experiment compares WSST1 and WSST2 over 1.35–2.50 kHz with a
64-voice CWT and a color scale normalized inside that band. The subtraction
experiment reconstructs the dominant WSST2 component from a ±150 Hz complex
ridge band, subtracts it in the time domain, projects only numerical
negative-frequency leakage out of the complex residual, and recomputes both
WSST orders. The plots are
`examples/soultanis_focused_secondary_wsst.png` and
`examples/soultanis_dominant_subtraction.png`; quantitative diagnostics are in
`validation/soultanis_secondary_modes.json`.

These experiments deliberately do not label candidate ridges as
``f_\mathrm{spiral}`` or ``f_{2-0}``. Several candidates change frequency or
supported duration after subtraction, and the 0–15 ms residual norm remains
roughly 47–68% of the analytic-input norm across the sequence. Thus,
subtraction reveals additional morphology but has not yet established that it
improves physical mode recovery. A synthetic injection/oracle-subtraction
study is still required for that conclusion.

The follow-up `examples/soultanis_figure5_bands.jl` replaces the single broad
secondary search with independent Figure 5-guided searches for `fspiral`,
`f2-0`, `f2+0`, and, for the high-mass cases where it is shown, `fspiral-0`.
Each interval is reassigned separately with both WSST orders and then drawn on
common time-frequency axes using the paper's Figure 2 color convention. A
Fourier-spectrum column exposes the exact prior applied to each waveform. See
`examples/soultanis_figure5_secondary_composite.png` and
`validation/soultanis_figure5_bands.json`.

The Fourier column uses the paper's
``h_{\mathrm{eff},+}=f|\widetilde h_+|`` definition, converts the supplied
``r h_+`` data to strain at 40 Mpc, uses the post-merger interval identified by
the Figure 2 caption, and uses a logarithmic physical-amplitude axis. A
right-hand cosine/Tukey taper rolls off the final 10% without suppressing the
physical merger onset. This is a preprocessing choice, not a paper parameter.
The Tukey window discussed by the paper belongs to its analytic waveform model,
while the separately mentioned unspecified window belongs to the
lapse-function spectrum. The secondary dashed lines reproduce
positions digitized from Figure 5; local spectral maxima inside the shaded
bounds are retained separately as diagnostics. The cyan band is evaluated
from the paper's mass-dependent analytic drift relation and its line marks the
corresponding spectral maximum.

The paper separately states that peak extraction uses the full signal including
inspiral. Accordingly, the secondary dashed positions are reproduced from
Figure 5 rather than redefined from the displayed post-merger spectrum.

The secondary search ranges and shaded bounds remain approximate readings of
Figure 5 rather than author-supplied table values. The unshaded `fspiral` line
receives a documented ±120 Hz search window. Each WSST band is independently
normalized in the composite, so its color answers “is there coherent structure
within this prior?” rather than comparing physical amplitudes across modes.
Overlapping high-mass priors can also follow the same coefficient structure
and must not be interpreted as independently identified modes without an
injection study. The paper's detector-sensitivity curves are omitted.

## Performance cost

The committed one-thread measurements used Julia 1.13.0-rc1, FFTW with one
thread, and the same 1024-sample input for both orders:

| Stage | Median time | Allocated memory |
|---|---:|---:|
| CWT and first derivative | 4.82 ms | 31.62 MiB |
| WSST1 reassignment after CWT | 22.54 ms | 14.21 MiB |
| Extra WSST2 operator bundle, including CWT | 20.97 ms | 69.93 MiB |
| Second-order frequency from prepared operators | 29.77 ms | 15.86 MiB |
| WSST2 reassignment after prepared operators | 50.79 ms | 30.07 MiB |
| WSST1 end to end | 29.82 ms | 45.83 MiB |
| WSST2 end to end | 83.84 ms | 100.00 MiB |

These are reproducible local measurements rather than portable performance
promises. The WSST2 implementation is intentionally scalar and single-threaded
until the reference path has accumulated more independent validation.

Regenerate the comparison with:

```bash
julia --project=benchmarks benchmarks/wsst2.jl
julia --project=examples examples/wsst2_bns.jl
```
