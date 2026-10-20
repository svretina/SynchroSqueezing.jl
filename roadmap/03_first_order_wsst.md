# Milestone 03 — First-order WSST

## Objective

Implement the new SST-specific machinery on top of the selected backend: stable instantaneous-frequency estimation, complex-coefficient reassignment, and mathematically consistent inversion.

## Scope

- Phase estimator based on `Im(dW * conj(W)) / (2π * abs2(W))` above a documented validity threshold.
- Absolute and relative magnitude thresholds plus a validity mask; invalid estimates are never silently converted to zero frequency.
- Ascending physical-frequency grid with explicit bin edges.
- Nearest-bin reassignment and conservative linear-bin deposition. Both accumulate complex coefficients and account for out-of-range/nonfinite contributions in optional diagnostics.
- Quadrature/reconstruction weights derived from the backend’s actual CWT normalization.
- Full inverse WSST and a minimal public result type retaining only necessary metadata; optional retention of CWT, derivative, and phase diagnostics.

## Validation

Build deterministic analytical tests for zero input, complex and real tones, linear/quadratic/exponential chirps, amplitude modulation, noise, odd/even lengths, and near-DC/Nyquist behavior. Verify CWT derivative behavior, tone-frequency recovery, coefficient-conservation identities, Float32/Float64 paths, and interior-region full reconstruction. Start with scalar reference kernels; fuse phase and deposition only after their independent results agree.

Compare with a pinned `ssqueezepy` release in stages: wavelet response, scales/grid, raw CWT, derivative, instantaneous frequency, SST matrix, and inverse. Supply identical scales, thresholds, wavelet parameters, padding, and frequency bins wherever possible. Agreement with Python is evidence, not the sole correctness criterion.

## Acceptance gate

First-order WSST localizes isolated tones at the correct frequency, conserves accepted weighted coefficients within floating-point error, reconstructs well-resolved noiseless Float64 signals to a justified interior error target (initially below `1e-3`), and has reproducible Python comparison fixtures.

## Explicit non-goals

- Ridge extraction, multi-component ownership, GPU execution, FSST, and SST2.

# Review of Milestone 03 — First-order WSST

Verdict: the milestone is well scoped but not yet implementation-complete. It identifies the correct algorithms and validation goals, but several mathematical and architectural decisions are still underspecified. An agent following it literally could produce a numerically plausible SST that fails to invert correctly, mishandles low frequencies, or cannot be integrated cleanly with `ContinuousWavelets.jl`.

Readiness assessment

Needs revision

Scope definition

## Strong

Mathematical specification

## Incomplete

Backend integration

## Incomplete

Testing strategy

## Good foundation

I would preserve the existing scope, but strengthen it in nine places.

## 1. Critical missing requirements

| Priority | Issue                             | Required correction                                                                                                    |
| -------- | --------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| P0       | Inversion convention              | Derive and independently verify the reconstruction constant and discrete quadrature weights for the actual backend     |
| P0       | Real versus complex signals       | Define positive-frequency analytic reconstruction and complex-signal support explicitly                                |
| P0       | DC and low-frequency residual     | Specify what happens to coefficients not represented in the SST frequency range                                        |
| P0       | Backend compatibility             | Establish exact shape, normalization, scale ordering, derivative convention, and father-wavelet treatment              |
| P0       | Bin deposition                    | Distinguish preservation of weighted coefficient sums from preservation of physical frequency integrals                |
| P0       | Numerical stability               | Specify threshold semantics, zero denominators, nonfinite values, overflow, and scale-dependent coefficient magnitudes |
| P1       | Boundary effects                  | Define a cone-of-influence or reliability mask, including derivative reliability                                       |
| P1       | API contracts                     | Define return types, array orientation, allocation options, and reconstruction behavior                                |
| P1       | Convergence and metamorphic tests | Add tests that detect discretization errors even when Python gives similar results                                     |

### A. The largest gap is the definition of inversion

The milestone says:

> Quadrature/reconstruction weights derived from the backend’s actual CWT normalization.

This is necessary but insufficient.

For an \\(L^2\\)-normalized continuous wavelet transform,

\\[ W_x(a,b)=\int x(t)a^{-1/2}\psi^\*((t-b)/a)\\,dt, \\]

the scale-integrated expression contains

\\[ \int_0^\infty W_x(a,b)a^{-3/2}\\,da. \\]

For an \\(L^1\\)-normalized transform,

\\[ W_x^{(1)}(a,b)=\int x(t)a^{-1}\psi^\*((t-b)/a)\\,dt, \\]

the corresponding expression is

\\[ \int_0^\infty W_x^{(1)}(a,b)\frac{da}{a}. \\]

The inverse constant depends on the Fourier convention and the wavelet. For a real signal and positive-frequency analytic wavelet, a typical reconstruction has the form

\\[ x(b)=\frac{2}{C\_\psi} \operatorname{Re}\left[ \int_0^\infty W_x^{(1)}(a,b)\frac{da}{a} \right], \\]

with \\(C\_\psi\\) defined consistently with that convention.

The Python `ssqueezepy` implementation explicitly applies a wavelet admissibility coefficient and a factor of two for analytic-wavelet real-signal reconstruction. Its CWT uses an \\(L^1\\)-normalization option.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub



Required improvement: introduce a tested `ReconstructionConvention` object and explicitly compute the constant and discrete scale quadrature weights.

Do not hardcode a normalization factor based solely on a successful sinusoid test.

### B. ContinuousWavelets.jl is not necessarily normalized like ssqueezepy

`ContinuousWavelets.jl` exposes a wavelet scaling normalization and frame-bound configuration, as well as multiple inverse-frame choices. It is therefore unsafe to assume its coefficients are interchangeable with `ssqueezepy` coefficients.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub

+1



The milestone must require verification of the selected backend's actual generated filters and normalization.

Also, a backend's inverse CWT is not automatically a mathematical validation of the proposed inverse WSST. Inverse CWT may use a dual-frame operator while SST reconstruction relies on integrating reassigned coefficients.

These are different reconstruction mechanisms and must be tested independently.

### C. Full reconstruction must define treatment of omitted frequencies

The current statement that the SST should reconstruct the full signal is ambiguous.

Consider a real-valued signal

\\[ x(t)=c+\cos(2\pi f_0t). \\]

An SST containing only positive-frequency analytic oscillations will generally not reconstruct the constant offset \\(c\\) unless its low-frequency contribution is represented separately.

Likewise, coefficients discarded by thresholding or reassignment outside the requested frequency range cannot simply be recovered from the squeezed matrix.

I would introduce explicit reconstruction modes:

- `:oscillatory` — reconstruct the represented positive-frequency oscillatory content.
- `:full` — reconstruct the full input when the result contains the necessary low-frequency/residual information.
- `:residual` — recover the unsqueezed or excluded contribution, where it has been retained.

If full reconstruction cannot be guaranteed for a particular configuration, the API must report that rather than silently return an incomplete signal.

### D. Complex-signal reconstruction requires an explicit convention

A complex analytic signal cannot generally be reconstructed with the same real-part operation used for a real signal.

For analytic signals, reconstruction should retain the complex integrated coefficients with the corresponding one-sided normalization.

For arbitrary complex signals containing both positive- and negative-frequency content, positive-frequency-only WSST is insufficient.

Therefore, Milestone 03 should support real and positive-frequency analytic complex signals initially, and either reject general two-sided complex input or require an explicitly documented two-sided transform.

This is better than promising unrestricted complex-signal support without defining how negative frequencies are handled.

### E. Conservative deposition needs a precise definition

Suppose scale \\(a_k\\) contributes

\\[ q\_{k,n}=w_kW(a_k,t_n). \\]

A nearest-bin kernel deposits the full contribution into one bin. A linear-bin kernel deposits fractions \\(\alpha q\_{k,n}\\) and \\((1-\alpha)q\_{k,n}\\) into neighboring bins.

This conserves the complex coefficient sum:

\\[ \sum_j T\_{j,n} = \sum\_{k\in\mathrm{accepted}}q\_{k,n}. \\]

But it does not necessarily conserve energy,

\\[ \sum_j |T\_{j,n}|^2, \\]

because complex coefficients can interfere.

It also does not imply conservation of a frequency integral unless the stored values are defined consistently as bin-integrated coefficients or frequency densities.

Recommendation: store the SST matrix as bin-integrated complex contributions. This makes the inverse a sum rather than an additional frequency quadrature, simplifying correctness and interoperability.

### F. The validity threshold needs more precision

The phase estimator

\\[ \widehat f= \frac{\operatorname{Im}(\partial_tW\overline W)} {2\pi |W|^2} \\]

is correct for the chosen time-phase convention, but the numerical contract is incomplete.

Specify:

- Whether the relative threshold is global, per scale, or per time column.
- Whether thresholding uses amplitudes or squared amplitudes.
- How the denominator is checked before evaluation.
- How NaNs and infinities propagate into diagnostics.
- Whether invalid estimates are represented by a mask or NaN.
- Whether coefficients whose frequencies lie outside the output grid are invalid or merely excluded from deposition.

A coefficient can have a perfectly valid phase estimate yet fall outside the requested output band. These are different statuses and should be tracked separately.

### G. The derivative must have known physical units

The transform derivative should be with respect to physical time in seconds:

\\[ \partial_tW. \\]

If the backend computes derivatives in sample-index coordinates,

\\[ \partial_tW=f_s\\,\partial_nW. \\]

Failing to apply this conversion changes the frequency estimates by a factor of the sampling rate.

The milestone should explicitly verify this using two equivalent discretizations of the same physical sinusoid.

### H. Error thresholds need a convergence protocol

An unconditional inverse-error threshold of \\(10^{-3}\\) across all listed signals is not appropriate.

Near-DC signals, strongly chirping components, signals close to Nyquist, and sharp transients can have significantly different reconstruction errors.

Separate:

1. Algebraic identities, which should agree to floating-point precision.
2. Numerical quadrature convergence, whose error decreases with increasing scale resolution.
3. Reconstruction accuracy, evaluated on well-resolved signals in valid interior regions.
4. Known limitation cases, which should test robustness rather than require high reconstruction fidelity.

### I. Correctness and optimizations need separate interfaces

The milestone mentions scalar kernels and subsequent fusion, which is good.

But specify that both implementations remain testable:

```
synchrosqueeze_reference!(...)
synchrosqueeze!(...)
```

The optimized implementation must agree with the reference on randomized coefficient matrices, not merely produce similar ridge plots.

This catches deposit indexing errors, race conditions, and threshold inconsistencies that signal-level tests may not expose.

## 2. Recommended implementation contract

The most important architectural decision is to define the SST as a transformation of an existing backend result, rather than tie the reassignment code to a specific CWT library.

ContinuousWavelets.jl backend

CWT coefficients, scales, frequency mapping, normalization

Validated CWT adapter

`W`, `dW`, scales, quadrature weights, boundary metadata, reconstruction constant

Phase estimator

frequency + validity

Reassignment

complex deposition

SSTResult

Bin-integrated coefficients + reconstruction and validity metadata

Inverse WSST

Oscillatory reconstruction / full reconstruction with retained residual

### Minimal proposed API

```

using Synchrosqueezing

# The CWT configuration is provided by the
# previously implemented backend adapter.
result = wsst(
    x,
    backend;
    fs=8192.0,
    freqrange=(100.0, 3500.0),
    nfreq=512,
    frequency_spacing=:log,
    reassignment=:nearest,
    threshold=RelativeThreshold(1e-8),
    keep_cwt=false,
    keep_derivative=false,
    keep_phase=false,
    preserve_residual=true
)

# Complex, bin-integrated SST coefficients
Tx = result.coefficients

# Ascending frequencies in Hz
freqs = result.frequencies
edges = result.frequency_edges

# Reconstruct available oscillatory content
xosc = inverse_wsst(result; mode=:oscillatory)

# Reconstruct full signal if residual retained
xfull = inverse_wsst(result; mode=:full)

# Optional standalone phase estimate
fhat, valid = instantaneous_frequency(W, dW; fs=8192.0)

```

A `SSTResult` should expose input sampling information, wavelet/backend metadata, output frequency grid, reconstruction convention, retained residual, and optional diagnostics.

The exact Julia types should be designed in this milestone, but public functions should not depend on undocumented internals of `ContinuousWavelets.jl`.

## 3. Fully revised Milestone 03

The following is the replacement specification I would give the development agent.

# Milestone 03 — First-order Wavelet Synchrosqueezing Transform (WSST)

## Objective

Implement a mathematically consistent, numerically stable, and invertible first-order wavelet synchrosqueezing transform in Julia, using the CWT backend selected in the preceding milestones.

The implementation must reuse `ContinuousWavelets.jl` or the selected equivalent backend for wavelet analysis. Do not introduce a separate general-purpose CWT implementation.

Milestone 03 must produce a complete and independently validated first-order WSST pipeline with correct physical frequency estimation, complex coefficient reassignment, explicit treatment of omitted signal content, and well-defined inverse reconstruction.

Ridge extraction and reconstruction of individual ridge-localized modes belong to later milestones.

---

## 0. Prerequisites and integration audit

Before implementation, verify that the preceding milestones expose:

- CWT coefficients for real and supported complex signals.
- Physical sample spacing or sampling frequency.
- Wavelet scale centers and ordering.
- Mother-wavelet Fourier convention.
- Coefficient normalization.
- Scale-frequency mapping.
- Time derivative of CWT coefficients.
- Original and padded signal lengths.
- Padding/cropping metadata.
- Supported wavelet families.
- Boundary reliability information, where available.
- Inverse CWT functionality.
- Low-frequency or father-wavelet treatment.

Write a backend-convention document describing these properties.

Specifically inspect whether `ContinuousWavelets.jl` uses a tight-frame or other discrete normalization, whether scale weights are implicit in its generated filters, and whether coefficient magnitudes vary with scale because of the chosen normalization.

Do not apply textbook CWT reconstruction weights blindly to backend coefficients.

If the backend interface cannot expose sufficient information for mathematically correct inversion, add the smallest necessary adapter functionality or improve the backend integration before proceeding.

### Deliverable

`docs/src/wsst_conventions.md`

It must specify all normalization, coordinate, and reconstruction conventions unambiguously.

---

## 1. Mathematical definitions

### 1.1 Continuous wavelet transform

Use the CWT definition and normalization selected in the preceding milestone.

For reference, a conventional L1-normalized analytic CWT is

[\
W_x(a,b)=\int x(t)a^{-1}\psi^\*((t-b)/a),dt.\
]

Its Fourier representation follows from the backend's Fourier convention.

No independent CWT implementation is required in Milestone 03.

### 1.2 Time derivative

The derivative must be

[\
D_x(a,b)=\partial_b W_x(a,b)\
]

with respect to physical time.

When computing derivatives spectrally, account for the backend's FFT convention and sampling interval.

Verify that the derivative of a pure complex exponential satisfies

[\
D_x(a,b)=i2\pi f_0W_x(a,b)\
]

away from boundary effects, for a positive-frequency analytic tone within the wavelet passband.

### 1.3 Instantaneous frequency

For each valid coefficient:

[\
\widehat f_x(a,b)=\
\frac{\operatorname{Im}\
[D_x(a,b)\overline{W_x(a,b)}]}\
{2\pi|W_x(a,b)|^2}.\
]

The estimate is measured in Hz.

This must be checked against the backend's complex phase convention. If the convention conjugates the time-phase dependence, correct the sign at the adapter level and document why.

### 1.4 SST coefficient construction

For the chosen CWT normalization, derive scale quadrature weights (w_k) such that the relevant continuous reconstruction integral is approximated by

[\
Q(b)=\sum_k w_kW_x(a_k,b).\
]

Each accepted coefficient contributes to the SST frequency bin corresponding to its estimated instantaneous frequency.

Define

\sum\_{k\in\mathcal A(j,n)}\
w_kW\_{k,n}\
]

for nearest-bin reassignment.

Define corresponding weighted deposition into two neighboring bins for linear reassignment.

The output matrix must store **bin-integrated complex coefficients**, not an unspecified spectral density.

### 1.5 Reconstruction

Derive the reconstruction relation

\mathcal R\left(\
C^{-1}\sum_jT\_{j,n}\
\right),\
]

where:

- (C) is the reconstruction constant appropriate to the backend.
- (\mathcal R) handles the analytic/real-signal convention.
- The summation covers represented positive-frequency oscillatory content.

For a real signal represented by positive-frequency analytic wavelets, (\mathcal R) typically includes twice the real part.

For a positive-frequency analytic complex signal, the reconstruction must retain its complex value using the appropriate one-sided normalization.

For arbitrary complex signals containing negative-frequency content, either provide a correct two-sided implementation or raise a clear unsupported-input error.

The distinction between an inverse CWT and inverse WSST must be documented and independently tested.

---

## 2. Threshold and validity infrastructure

Implement threshold policies as explicit, validated types.

Suggested types:

- `AbsoluteThreshold`
- `RelativeThreshold`
- `NoThreshold`

Relative thresholding must specify its reference: global transform maximum, scale-wise maximum, or another documented statistic.

Do not implement ambiguous threshold behavior.

Thresholds apply to coefficient magnitude (|W|), not implicitly to squared magnitude.

A point is valid for reassignment only if:

1. The transform coefficient is finite.
2. The derivative coefficient is finite.
3. The magnitude exceeds the selected threshold.
4. The denominator is numerically safe.
5. The resulting instantaneous frequency is finite.

A valid estimate falling outside the SST frequency grid must be classified separately as out of range.

Invalid instantaneous frequencies must never silently become zero.

Provide a standalone phase-estimation implementation returning the estimated frequencies and a validity mask.

For memory-sensitive execution, permit direct phase evaluation inside the deposition kernel without allocating the full phase matrix.

### Required edge cases

Test:

- Exact zero coefficient.
- Subnormal coefficient magnitude.
- Very small but nonzero coefficients.
- Zero derivative.
- Nonfinite coefficient or derivative.
- Large magnitudes approaching overflow risk.
- Negative-frequency estimate.
- Zero-frequency estimate.
- Frequencies outside the requested output interval.

Use stable numerical arithmetic and avoid unnecessary complex division.

---

## 3. Frequency-grid infrastructure

Provide construction and validation for the physical output frequency grid.

Support:

- Linearly spaced frequency centers.
- Logarithmically spaced frequency centers.
- Explicit user-supplied frequency centers.
- Explicit bin edges.
- Ascending frequency ordering.

Require strictly increasing and finite bin edges, with positive frequencies for the default one-sided WSST.

For logarithmic grids, all frequency centers must be positive.

Specify whether the user-provided frequency range describes bin centers or outer bin edges. The selected interpretation must be consistent throughout the API.

Implement fast frequency-to-bin mappings for uniform linear and logarithmic grids.

For arbitrary nonuniform grids, implement a correct search-based mapping.

Nearest-bin reassignment must use an explicitly documented nearest-center or containing-bin policy; do not conflate these definitions.

For linear deposition, calculate interpolation weights using neighboring physical frequency centers, and handle extreme bins consistently.

Clearly define behavior at exact boundaries and bin centers.

Do not clamp out-of-range frequencies silently into the first or last bin.

---

## 4. Reference reassignment algorithms

Implement two independent, readable reference routines.

### 4.1 Nearest-bin deposition

For each time index and wavelet scale:

1. Read the CWT coefficient.
2. Evaluate its validity.
3. Calculate the instantaneous frequency.
4. Locate the destination frequency bin.
5. Calculate the scale-weighted complex contribution.
6. Add the contribution to the selected output bin.
7. Record exclusions when diagnostics are enabled.

### 4.2 Linear-bin deposition

For each accepted contribution:

1. Find the neighboring frequency centers.
2. Calculate their interpolation weights.
3. Deposit fractions of the same complex contribution into the two bins.
4. Ensure the two weights sum to one.
5. Handle grid boundaries according to a documented policy.

Use the same threshold logic, reconstruction weights, and out-of-range behavior in both implementations.

### 4.3 Coefficient conservation

For every time column, verify

\sum\_{k\in\mathrm{deposited}}w_kW\_{k,n}\
]

within floating-point accumulation error.

This is a complex coefficient-sum identity. Do not interpret it as conservation of squared magnitude or physical energy.

A frequency-weighted integral should be tested only if the representation is explicitly stored as a frequency density.

---

## 5. Excluded coefficients and residual bookkeeping

Categorize coefficients into:

- Deposited successfully.
- Rejected by magnitude threshold.
- Rejected because of nonfinite input.
- Rejected because of unstable phase estimation.
- Excluded because their instantaneous frequency lies outside the output grid.

Track these categories in optional diagnostics.

Define a reconstruction residual that accounts for finite, valid weighted coefficients that were not deposited.

Do not attempt to silently reconstruct information from NaN or infinite coefficients.

Support:

- `preserve_residual=true`: retain the necessary reconstruction contribution.
- `preserve_residual=false`: permit lower-memory operation without promising full reconstruction.

Low-frequency or father-wavelet information requires special handling.

When such information is provided by the backend and is necessary for reconstructing the original input, preserve it separately from the oscillatory SST coefficients.

Do not assign father-wavelet coefficients an artificial zero-frequency ridge.

Define clear behavior for a signal with a nonzero constant offset.

---

## 6. Inverse WSST

Implement:

- `inverse_wsst(result; mode=:oscillatory)`
- `inverse_wsst(result; mode=:full)`

The oscillatory inverse reconstructs the content represented in the squeezed matrix.

The full inverse includes the retained low-frequency contribution and valid excluded-coefficient residual, when available.

If the result lacks the information required for full inversion, throw an informative error.

### Reconstruction requirements

- Correct handling of real signals.
- Correct handling of supported analytic complex signals.
- Explicit rejection or correct treatment of two-sided complex signals.
- Correct normalization and quadrature.
- Correct array dimensions and element types.
- No dependence on ridge information.
- No reapplication of frequency-bin widths to bin-integrated SST coefficients.

Validate inversion independently of the backend's inverse CWT.

In addition, compare the full inverse WSST against both the input signal and the backend's inverse CWT reconstruction.

Discrepancies must be analyzed in terms of normalization, scale truncation, frame effects, thresholding, frequency exclusion, and boundaries.

---

## 7. Public API and result types

Design a minimal but extensible API.

Suggested entry points:

- `wsst(x, backend; kwargs...)`
- `instantaneous_frequency(W, dW; kwargs...)`
- `synchrosqueeze(W, dW, metadata; kwargs...)`
- `inverse_wsst(result; mode=...)`

Provide an in-place reassignment method and a reference implementation separately.

### Result structure

Define a concrete, parametric `SSTResult`.

It must retain or expose:

- Bin-integrated complex SST coefficients.
- Ascending physical frequency centers.
- Frequency-bin edges.
- Sampling frequency and interval.
- Original sample count.
- Wavelet and backend configuration.
- CWT scale coordinates.
- Scale quadrature and reconstruction convention.
- Threshold and reassignment configuration.
- Padding and crop metadata.
- Signal type and analytic/real reconstruction convention.
- Low-frequency reconstruction information, if retained.
- Residual reconstruction information, if retained.
- Optional phase, validity, CWT and derivative arrays.
- Boundary reliability information, where available.

Avoid unnecessary duplication of metadata already held in immutable backend configuration objects.

### Storage orientation

Use frequency × time as the public SST matrix orientation.

Backend CWT arrays may use time × scale orientation.

Do not force a full matrix copy merely to reconcile the two representations; handle layouts through adapters, views, or specialized kernels.

### Error handling

Validate and report:

- Unsupported wavelet configuration.
- Invalid sampling frequency.
- Nonuniform input sampling without a supported resampling stage.
- Empty or invalid frequency grid.
- Inconsistent transform and derivative dimensions.
- Unsupported input signal type.
- Missing reconstruction metadata.
- Incompatible normalization.
- Invalid threshold configuration.

---

## 8. Boundary effects and numerical reliability

A finite signal does not provide equally reliable wavelet coefficients at every time.

Use backend wavelet-support and padding information to calculate or approximate a frequency-dependent reliability mask.

Where appropriate, expose a cone-of-influence diagnostic.

If the backend does not supply sufficient support metadata, document and implement a conservative support estimate.

The reliability mask must be separate from the phase-validity mask:

- Phase validity concerns numerical evaluation of the instantaneous frequency.
- Boundary reliability concerns contamination from finite support, padding, or spectral differentiation.

Neither mask should silently alter results unless a filtering option explicitly requests it.

Test periodic, zero, and reflect padding where supported.

Verify that numerical differentiation does not introduce unexpected endpoint artifacts into otherwise reliable interior samples.

---

## 9. Analytical testing

Use deterministic fixtures with fixed random seeds.

### 9.1 Elementary inputs

Test:

- Zero signal.
- Constant nonzero signal.
- Real sinusoid.
- Complex positive-frequency exponential.
- Unsupported or explicitly two-sided negative-frequency complex exponential.
- Two separated sinusoidal components.
- Amplitude-modulated sinusoid.

### 9.2 Chirps

Test:

- Linear chirp.
- Quadratic instantaneous-frequency chirp.
- Exponential chirp.
- Chirps with changing amplitude.
- Chirps near the supported frequency boundaries.

Use exact analytical phase and instantaneous-frequency functions.

### 9.3 Numerical conditions

Test:

- Odd and even signal lengths.
- Power-of-two and non-power-of-two lengths.
- Float32 and Float64.
- Very small signal amplitudes.
- Signal amplitude scaling.
- Nonfinite transform coefficients.
- Low signal-to-noise ratios.
- Near-DC and near-Nyquist tones.
- Different sampling rates describing the same continuous signal.

### 9.4 Algebraic invariants

Verify:

1. The phase estimator returns the correct frequency for ideal complex tones.
2. The validity mask rejects ill-conditioned coefficients.
3. Deposited complex coefficients satisfy the conservation identity.
4. Nearest-bin and linear-bin deposition agree on exact bin-center frequencies.
5. Linear interpolation weights sum to one.
6. Linear deposition preserves the first frequency moment when both neighboring centers are available.
7. The full inverse is insensitive to the location of deposited coefficients when all their contributions remain included.
8. The fused optimized kernel agrees with the scalar reference.
9. Transform results are deterministic under the supported execution model.
10. Output arrays and frequency grids have documented shapes and ordering.

### 9.5 Convergence studies

Perform systematic refinements over:

- Number of wavelet scales.
- Voices per octave.
- Output frequency-bin resolution.
- Signal sampling rate.
- Padding size or boundary support.
- Float32 versus Float64 precision.

Report instantaneous-frequency error and reconstruction error against resolution.

Distinguish errors caused by CWT discretization from those caused by synchrosqueezing.

Do not claim convergence without observing the expected numerical trend in an appropriate smooth test case.

---

## 10. Python cross-validation

Pin the Python version and commit or release used for `ssqueezepy`.

Generate reference fixtures including:

- Input signal.
- Wavelet parameters and Fourier response.
- Scale coordinates.
- Raw CWT.
- CWT derivative.
- Instantaneous-frequency estimates.
- Threshold and validity information.
- SST frequency coordinates.
- Complex SST coefficients.
- Full inverse reconstruction.

Compare in stages:

1. Wavelet response.
2. CWT coefficients.
3. Derivative coefficients.
4. Instantaneous-frequency estimates.
5. Reassignment destination bins.
6. SST coefficient matrix.
7. Coefficient-conservation identities.
8. Reconstructed signal.

Align wavelet normalization, CWT scale grids, padding, thresholding, and output frequency bins whenever possible.

Do not demand direct elementwise equality between implementations using different wavelet frames or discretizations.

For incompatible conventions, compare both implementations to independently known analytical quantities.

Include one synthetic fixed-coefficient fixture bypassing the CWT backend, to test only the phase and reassignment algorithms.

Document discrepancies with actual numerical measurements rather than qualitative visual statements.

---

## 11. Performance requirements

First implement a straightforward scalar reference algorithm.

Then introduce an optimized kernel with:

- In-place output accumulation.
- Reusable scratch buffers.
- No unnecessary phase matrix allocation.
- Avoidance of per-coefficient heap allocations.
- Specialized bin indexing for regular grids.
- Layout-aware iteration.
- Type-stable code.
- Optional CWT/derivative retention.

Keep the reference implementation available for testing.

Measure separately:

- Phase estimation.
- Frequency-to-bin mapping.
- Nearest-bin deposition.
- Linear-bin deposition.
- Complete WSST after the CWT.
- End-to-end WSST including the backend CWT.

Record Julia allocations and warmed runtimes using BenchmarkTools.jl.

Do not implement GPU kernels, distributed execution, or complex thread synchronization in this milestone.

Threading can be considered only if it is already part of the established package infrastructure and introduces no additional validation risk.

---

## 12. Required deliverables

### Source code

- Phase estimator and validity handling.
- Threshold policies.
- Frequency-grid implementation.
- Scale quadrature and reconstruction convention.
- Nearest-bin deposition.
- Linear-bin deposition.
- Reference reassignment kernel.
- Optimized reassignment kernel.
- SST result type and metadata.
- Low-frequency and residual accounting.
- Inverse WSST.

### Tests

- Standalone phase-estimator tests.
- Frequency-grid and interpolation tests.
- Reassignment conservation tests.
- Reconstruction tests.
- Numerical precision tests.
- Analytical signal tests.
- Convergence studies.
- Python reference comparisons.

### Documentation

Create:

- `docs/src/wsst_theory.md`
- `docs/src/wsst_conventions.md`
- `docs/src/wsst_inverse.md`
- `docs/src/wsst_numerics.md`
- `docs/src/wsst_validation.md`

Explain all mathematical assumptions, supported signal classes, normalization conventions, finite-frequency coverage limitations, and inverse behavior.

### Examples

Create:

- `examples/wsst_tone.jl`
- `examples/wsst_chirp.jl`
- `examples/wsst_reconstruction.jl`

Examples must be runnable without Python.

---

## 13. Acceptance gate

Milestone 03 is complete only when the following criteria are satisfied.

### Mathematical correctness

- The phase estimator recovers the frequency of a well-resolved complex tone in Hz.
- Real-signal reconstruction uses the mathematically correct analytic-wavelet convention.
- Complex-signal reconstruction is correct for the supported signal class.
- CWT normalization and quadrature weights are documented and numerically verified.
- Reassignment preserves the sum of accepted complex weighted coefficients.
- Full inverse reconstruction correctly handles retained low-frequency and excluded finite contributions.

### Numerical verification

- Analytical tests pass.
- Nearest-bin and linear-bin kernels pass all algebraic consistency tests.
- Scalar and optimized implementations agree within justified floating-point tolerances.
- Tests pass for Float32 and Float64 with precision-appropriate thresholds.
- Well-resolved noiseless Float64 signals achieve a justified interior reconstruction target, provisionally (E_x < 10^{-3}).
- Reconstruction accuracy and instantaneous-frequency accuracy exhibit appropriate convergence under refinement.
- Boundary error is measured separately.

### Interoperability

- Python comparison fixtures are reproducible.
- Differences between Julia and Python are explained by documented implementation conventions or measured numerical errors.
- The Julia package does not require Python at runtime.

### Software quality

- Public APIs and data types are documented.
- Invalid inputs and unsupported modes fail explicitly.
- No undefined frequency estimates are silently converted to physical zero frequency.
- The primary implementation does not duplicate the selected CWT backend.
- Automated tests pass in CI.

---

## Explicit non-goals

- Automatic ridge extraction.
- Multiridge tracking.
- Ridge-band component reconstruction.
- Mode ownership or classification.
- Second-order synchrosqueezing.
- STFT-based synchrosqueezing.
- GPU acceleration.
- General-purpose wavelet implementations.
- Advanced multithreaded optimization.
- Scientific fitting of physical BNS mode-frequency relations.

## Final execution requirement

Implement the milestone incrementally.

Begin by deriving the discrete reconstruction formula from the actual CWT backend and confirming its consistency with analytical test signals.

Then implement the scalar phase and reassignment kernels, followed by inversion and residual accounting.

Only after the complete reference pipeline is validated should the optimized reassignment kernel be introduced.

Do not mark the milestone complete based solely on attractive spectrograms or Python agreement.

Completion requires numerical evidence of mathematical consistency, reconstruction accuracy, and reproducibility.

## 4. Two final recommendations

### Do not let residual bookkeeping conceal a broken inverse

There is a subtle but important testing danger in the revised milestone: if `preserve_residual=true` stores enough information to reproduce the original signal, then full reconstruction could appear perfect even when the SST inversion is incorrectly normalized.

Require three separate reconstruction tests:

| Test                                                                         | What it validates                                       |
| ---------------------------------------------------------------------------- | ------------------------------------------------------- |
| Inverse of all deposited SST coefficients                                    | SST normalization, quadrature, and complex accumulation |
| SST reconstruction plus independently computed excluded-coefficient residual | Complete represented CWT reconstruction                 |
| Reconstruction including the backend's low-frequency channel                 | Recovery of the original signal, including DC           |

The residual should be calculated from the weighted CWT coefficients excluded from deposition, not from the difference between the original signal and the SST reconstruction.

This ensures reconstruction tests actually validate the numerical algorithms.

### Add an inverse-consistency test that is independent of frequency assignment

For a given CWT and fixed set of accepted coefficients, compute both nearest-bin and linear-bin synchrosqueezing.

Because both preserve the same complex weighted sum, their full oscillatory inverses should agree to rounding error, even though their time-frequency representations differ.

This is an especially powerful test for deposition bugs because it simultaneously checks indexing, interpolation, accumulation, and inversion.

It also provides a foundation for later ridge reconstruction: once the global inverse is correct, the next milestone can investigate which subsets of frequency bins correspond to individual physical components.

Overall: with these changes, Milestone 03 becomes a complete specification for a scientifically reliable first-order WSST, rather than just an implementation checklist. I would keep ridge extraction outside this milestone and require rigorous inversion and cross-backend validation before moving to Milestone 04.