# Amplitude and damping estimation

Milestone 04C estimates damping from a reconstructed complex component, not
from WSST ridge brightness. For a ridge reconstruction ``z(t)``, the analysis
uses

```math
A(t)=|z(t)|,\qquad \phi(t)=\operatorname{unwrap}(\arg z(t)),\qquad
f(t)=\frac{1}{2\pi}\frac{d\phi}{dt}.
```

For real input, request `analytic=true` during reconstruction. The resulting
complex signal retains the positive-frequency component; its real part equals
the ordinary real reconstruction and its magnitude is the envelope.

```julia
component = reconstruct_component(
    result, ridge; bandwidth_hz = 6.0, analytic = true,
)
features = component_features(
    component.signal, times;
    valid_mask = ridge.valid .& ridge.active,
    amplitude_floor = 1e-3,
)
fit = fit_exponential_damping(
    features;
    window = (0.25, 2.0),
    method = :loglinear,
)
```

## Constant damping fit

The model is

```math
A(t)=A_r\exp[-\gamma(t-t_r)],\qquad \tau=\gamma^{-1}.
```

The default `:loglinear` method performs weighted least squares on
`log(A)`. The reference time is the weighted mean time, which improves
conditioning. `:nonlinear` is an independent Gauss–Newton fit in amplitude
space and can be useful when additive envelope noise makes the logarithm
strongly heteroscedastic.

`DampingFit` contains the damping time and rate, fitted envelope, amplitude
and log-amplitude residuals, exact selection mask, effective time window,
RMSE, ``R^2``, convergence state, and linearized standard errors. A
nonpositive decay rate is flagged as `status == :nondecaying`; when positive
damping is required its `damping_time` is `NaN`.

Weights, an explicit fitting window, a validity mask, and an amplitude floor
control sample inclusion. This matters when a mode becomes inactive or fades
into the noise. Excluded samples remain visible in `fit.selected`; they are
never silently reintroduced.

## Local rate

`local_damping_rate(features)` estimates

```math
\gamma(t)=-\frac{d\log A(t)}{dt}.
```

This is a time-dependent diagnostic, not another estimate of a constant
``\tau``. Differentiation amplifies noise, so an odd `smoothing_window` may be
used, and the result should be interpreted only where the envelope is safely
above its reliability floor.

## Separate uncertainty sources

`fit.damping_time_standard_error` is the regression uncertainty conditional
on the selected ridge, reconstruction band, and fit interval. It does not
measure ridge or reconstruction uncertainty.

Use `damping_sensitivity` with alternate ridges, bandwidths, and windows to
measure that separate source:

```julia
sensitivity = damping_sensitivity(
    result,
    [baseline_ridge, alternate_ridge];
    bandwidths_hz = [5.0, 6.0, 7.0],
    windows = [(0.25, 1.8), (0.35, 2.0)],
)
```

Each case retains its complete fit. The reported minimum, maximum, and spread
describe between-case sensitivity. Nearby-mode leakage, amplitude modulation,
low SNR, boundary transients, or a poor ridge can produce a visibly
nonexponential envelope and residual structure even when the formal fit
standard error is small.

## Two damped oscillators

Run the complete CairoMakie validation example with:

```bash
julia --project=examples examples/damped_oscillators.jl
```

It generates two damped chirps with different damping times, extracts both
ridges, reconstructs analytic components, performs damping fits, evaluates
bandwidth/window sensitivity, and writes `examples/damped_oscillators.png`.
The plot and terminal report include ridge RMSE, damping-time error,
reconstruction similarity, relative ``L^2`` error, and fit uncertainty.

## Fade-out stress test

When one mode has a much shorter damping time, its ridge tracker still returns
a mathematical path after the physical component has fallen into the noise.
That path must not be interpreted as continuing evidence for the mode. Use an
explicit activity or envelope floor, report the last reliable time, and fit
only the supported interval.

```bash
julia --project=examples examples/damping_fadeout_stress.jl
```

This example uses ``\tau_1=0.32`` s and ``\tau_2=2.60`` s with additive white
noise. It compares a reliability-aware fit against a fit whose end time is
progressively extended beyond the short mode's known synthetic detectability
limit. The latter exposes the expected noise-floor bias rather than hiding it
inside a single damping estimate.

## Operational mode termination

`ridge_prominence` compares the mean squared WSST magnitude in a narrow band
around a ridge with two neighboring sidebands. This is local transform
contrast—not detector SNR and not a calibrated statistical significance.
Both sidebands must exist, so ridges too near a frequency boundary are marked
invalid independently.

```julia
prominence = ridge_prominence(
    result,
    ridge;
    ridge_bandwidth_hz = 1.0,
    background_band_hz = (2.5, 7.0),
)
```

`mode_support` combines the reconstructed component's amplitude relative to
its own reference with an optional prominence threshold:

```julia
support = mode_support(
    features;
    prominence,
    relative_amplitude_threshold = 0.05,
    prominence_threshold = 100.0,
    persistence_seconds = 0.08,
    recovery_factor = 1.25,
    reference_mask = early_reliable_samples,
    boundary_mask = wavelet_reliable_samples,
)
fit = fit_exponential_damping(features; valid_mask = support.fit_mask)
```

The persistence interval prevents a single bad sample from terminating a
mode. `recovery_factor` adds hysteresis while termination remains pending.
Once termination is confirmed, `support.supported` stays false. The original
ridge is not altered and remains available for debugging. `stop_reason`
distinguishes relative-amplitude, insufficient-prominence, combined,
invalid-component, boundary, and end-of-data stops. A mode supported through
the record has `stop_reason == :end_of_data` and `stop_time === nothing`.
`fit_mask` is stricter than the state mask: it excludes transient failed
samples as well as everything after termination.

`supported_reconstruction` packages the untouched raw reconstruction together
with a separately masked signal. `mode_support_sensitivity` repeats support
selection and damping fitting over several amplitude thresholds; its damping
spread measures decision sensitivity, not conditional regression uncertainty.

The thresholds are operational analysis choices and are retained in
`support.config`. The amplitude reference is component-specific; it is never
the maximum of the full multi-component waveform.

Run the two-component comparison with:

```bash
julia --project=examples examples/mode_termination_two_components.jl
```

It compares 2%, 5%, and 10% amplitude cutoffs and their damping fits, tests
prominence-only thresholds independently, applies an absorbing
5%-plus-prominence cutoff, preserves the raw DP path and reconstruction, and
restricts the usable component and damping estimate to the supported prefix.

The companion `examples/mode_termination_validation.jl` verifies a
nonexponential envelope and sweeps the separation between a weak short-lived
mode and a 3.4-times stronger long-lived mode. Close modes cause premature
termination and biased or invalid damping estimates in this deliberately hard
case. This limitation is reported directly rather than hidden by selecting a
threshold against synthetic truth.
