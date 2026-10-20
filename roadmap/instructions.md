# Development specification: Synchrosqueezing.jl

Standalone Julia package · Scientific computing · Time-frequency analysis · Ridge extraction and signal reconstruction

## 1. Project objective

Build a high-performance, pure-Julia package named `Synchrosqueezing.jl` implementing synchrosqueezed time-frequency transforms, instantaneous frequency estimation, automatic ridge tracking, and reconstruction of individual oscillatory components.

The package must be scientifically validated, numerically accurate, modular, and substantially optimized for repeated analysis of large datasets.

It should support arbitrary signals, but its architecture should accommodate applications to gravitational-wave signals from binary neutron star (BNS) post-merger simulations, where several simultaneous oscillatory modes exhibit rapidly evolving frequencies.

The package must remain independent of Python at runtime. Python will be used exclusively to generate reference results, perform cross-implementation comparisons, and measure performance.

### Principal reference implementation

The most useful Python baseline is [ssqueezepy](https://github.com/OverLordGoldDragon/ssqueezepy?utm_source=chatgpt.com). It provides CWT, STFT, their synchrosqueezed versions, inverse transformations, generalized Morse wavelets, and automatic ridge extraction. Its numerical kernels already use Numba, parallel execution, and optional GPU acceleration. Consequently, outperforming it should be treated as an experimentally testable objective, not an assumed advantage of Julia.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub

+1

The proposed package should implement the following.

| Capability                                    | Priority |
| --------------------------------------------- | -------- |
| FFT-based continuous wavelet transform (CWT)  | P0       |
| First-order synchrosqueezed CWT (WSST)        | P0       |
| Instantaneous frequency estimation            | P0       |
| Dynamic-programming ridge extraction          | P0       |
| Component and full-signal reconstruction      | P0       |
| Synthetic analytical verification             | P0       |
| Python `ssqueezepy` reference comparison      | P0       |
| Short-time Fourier transform SST (FSST)       | P1       |
| Multiple simultaneous and crossing ridges     | P1       |
| Ridge uncertainty and confidence diagnostics  | P1       |
| Second-order synchrosqueezing (SST2)          | P2       |
| Iterative ridge refinement and reconstruction | P2       |
| Threading, GPU support, batch processing      | P2–P3    |

Here P0 defines the first scientifically usable release, not merely a software demonstration.

## 2. Scientific foundations

### 2.1 Signal model

Consider a signal decomposed into amplitude- and frequency-modulated components:

\\[ x(t)=\sum\_{k=1}^{K} A_k(t)e^{i\phi_k(t)}+\epsilon(t), \\]

where

\\[ f_k(t)=\frac{\phi_k'(t)}{2\pi} \\]

is the instantaneous frequency of mode \\(k\\), \\(A_k(t)\\) is its amplitude, and \\(\epsilon(t)\\) represents noise or unresolved structure.

For real signals, the physically observed signal may be the real part of this expression.

The purpose of synchrosqueezing is to concentrate a conventional time-frequency representation around the curves \\(f_k(t)\\), making it easier to estimate the frequencies and reconstruct their associated components.

Two important assumptions underpin standard first-order SST:

- Each mode is sufficiently well represented locally by a slowly varying amplitude and phase.
- Distinct components are sufficiently separated in the time-frequency representation for reliable localization.

These conditions can break down during mode crossings, strong transients, abrupt frequency changes, and severe overlap. The package must explicitly test these limitations rather than imply that SST automatically resolves all overlapping signals.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://www.sciencedirect.com\&sz=32)

ScienceDirect

+1

### 2.2 Continuous wavelet transform

Use the convention

\\[ W_x(a,b)=\int\_{-\infty}^{\infty} x(t)\frac{1}{\sqrt{a}} \psi^\*\left(\frac{t-b}{a}\right)\\,dt, \\]

where \\(a>0\\) is scale and \\(b\\) is time.

The package must also support an \\(L^1\\)-normalized convention, with normalization \\(1/a\\), because this facilitates comparison with `ssqueezepy`.

The Fourier-domain implementation is

\\[ W_x(a,b)=\mathcal F^{-1} \left[ \widehat{x}(\omega)\sqrt{a}\\, \widehat{\psi}^{\\,\*}(a\omega) \right]\(b). \\]

For the \\(L^1\\) convention, the scale factor changes accordingly.

The implementation should initially support a complex analytic Morlet wavelet, followed by generalized Morse wavelets.

For Morlet, a starting definition is

\\[ \widehat{\psi}(\omega)= C e^{-(\omega-\omega_0)^2/2}\Theta(\omega), \\]

with explicit normalization and handling of the zero-frequency response. A strictly admissible implementation must address the zero-mean condition rather than silently treating every Gaussian frequency response as admissible.

Engineering requirement: wavelet normalization, angular versus cyclic frequency, physical sampling interval, and FFT sign conventions must all be documented and independently tested.

### 2.3 Instantaneous frequency estimation

For sufficiently large wavelet coefficients, estimate

\\[ \boxed{ \widehat f_x(a,b) =\frac{1}{2\pi} \operatorname{Im} \left[ \frac{\partial_b W_x(a,b)} {W_x(a,b)} \right] } \\]

in Hz when the derivative is taken with respect to physical time.

The preferred derivative is computed spectrally:

\\[ \partial_b W_x= \mathcal F^{-1} \left[ i\omega\\,\widehat{x}(\omega) \sqrt{a}\widehat{\psi}^{\\,\*}(a\omega) \right]. \\]

This avoids finite-difference errors in the phase derivative and permits computing the transform and derivative using related FFT kernels.

Near zeros of \\(W_x\\), the ratio becomes numerically unstable. Introduce a threshold \\(\gamma\\), and invalidate frequency estimates wherever

\\[ |W_x(a,b)|\leq\gamma. \\]

Do not substitute a frequency of zero for invalid estimates. Prefer a validity mask or NaN in diagnostic arrays.

### 2.4 Synchrosqueezing

The wavelet coefficients are reassigned from their original scale locations into frequency bins centered around the estimated instantaneous frequency.

In continuous form, schematically,

\\[ T_x(f,b)= \int_0^\infty W_x(a,b)a^{-3/2} \delta(f-\widehat f_x(a,b))\\,da. \\]

For discrete logarithmically spaced scales, an implementation can take the form

\\[ T_x(f_j,b) \mathrel{+}= W_x(a_k,b)\\,\Delta a_k\\,a_k^{-3/2}, \\]

subject to \\(\widehat f_x(a_k,b)\\) falling within frequency bin \\(j\\).

For the \\(L^1\\)-normalized transform, the corresponding integration weight changes; on logarithmic scales it is naturally expressed using \\(d\log a\\). These formulations must not be mixed.

The package should implement two distinct reassignment strategies:

1. Nearest-bin reassignment: fast and simple, used as the initial reference.
2. Conservative linear-bin deposition: distributes each contribution between adjacent output bins while preserving the sum of assigned complex coefficients.

Complex coefficients must be accumulated, not merely their squared magnitudes. This is essential for inverse reconstruction.

### 2.5 Ridge extraction

Given a complex synchrosqueezed representation \\(T(f_j,t_n)\\), define a ridge path

\\[ r=(j_1,\ldots,j_N) \\]

as the solution of

\\[ r^\star= \underset{j_1,\ldots,j_N}{\operatorname{argmin}} \left[ \sum\_{n=1}^{N}D(j_n,n) + \lambda \sum\_{n=2}^{N} \left(f\_{j_n}-f\_{j\_{n-1}}\right)^2 \right], \\]

where \\(D\\) penalizes low signal energy, and \\(\lambda\\) penalizes rapid changes in frequency.

A useful data cost is

\\[ D(j,n)=- \log\left( \frac{|T(f_j,t_n)|^2} {\max\_\ell |T(f\_\ell,t_n)|^2+\varepsilon} +\varepsilon \right). \\]

For logarithmically spaced frequency bins, support penalization in log-frequency instead of physical frequency.

This formulation permits dynamic programming and forward-backward ridge recovery. The Python reference uses a related algorithm with squared jump penalties and exclusion bands around previously extracted ridges.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub

+1

### 2.6 Reconstruction

For a ridge \\(r_k(t)\\), reconstruct the associated oscillatory component by integrating synchrosqueezed coefficients inside a band surrounding it:

\\[ \widetilde x_k(t)= \frac{1}{C\_\psi} \int\_{|f-f_k(t)|<\Delta f_k(t)} T_x(f,t)\\,df, \\]

where \\(C\_\psi\\) is the wavelet-dependent reconstruction constant consistent with the chosen CWT convention.

For real signals represented using positive-frequency analytic wavelets, the physical reconstruction involves the appropriate real part and factor of two.

The package must provide full reconstruction, ridge-band reconstruction, and residual reconstruction. It must also demonstrate numerically that

\\[ x(t)\approx\sum_k x_k(t)+x\_{\rm residual}(t) \\]

within the discretization error and the limits imposed by the frequency coverage and edge effects.

Critical: correct reconstruction weights and normalizations must be derived and tested for each transform, rather than fitted to match Python output.

### 2.7 Second-order synchrosqueezing

First-order SST is exact in its instantaneous-frequency estimate for ideal locally monochromatic signals, but a rapidly varying chirp can produce frequency smearing.

Second-order SST attempts to correct this by incorporating local phase curvature, using derivatives of the STFT or CWT with respect to time and frequency/scale.

For the STFT, let

\\[ V_x(t,\eta)=\int x(u)g(u-t)e^{-i2\pi\eta(u-t)}\\,du. \\]

Its first-order frequency estimate is

\\[ \omega_x(t,\eta) =\frac{1}{2\pi} \operatorname{Im}\left(\frac{\partial_tV_x}{V_x}\right), \\]

for this particular demodulated-window convention.

Introduce the complex quantities

\\[ \widetilde\omega_x= \frac{1}{2\pi i}\frac{\partial_tV_x}{V_x}, \\]

\\[ \widetilde\tau_x= t-\frac{1}{2\pi i} \frac{\partial\_\eta V_x}{V_x}. \\]

The local chirp-rate correction can be expressed schematically through

\\[ q_x= \frac{\partial_t\widetilde\omega_x} {\partial_t\widetilde\tau_x}, \\]

giving

\\[ \omega_x^{(2)} = \operatorname{Re} \left[ \widetilde\omega_x+ q_x(t-\widetilde\tau_x) \right]. \\]

The complete discretized formulation needs to be derived carefully from one consistent STFT convention, including the derivatives, window moments, denominator thresholds, and treatment of degenerate points.

Do not implement SST2 by simply adding a numerical second derivative to the first-order instantaneous-frequency estimate. It requires the correct reassignment operators and stability conditions.

For fast-evolving BNS post-merger frequencies, second-order SST is a particularly worthwhile research extension, but first-order SST and inversion must be validated before introducing it.

## 3. Repository architecture

The implementation must separate mathematical definitions, numerical algorithms, data management, and high-performance execution kernels.

Use the following repository layout.

Synchrosqueezing.jl/

Package root and environment

├── Project.toml ├── README.md ├── LICENSE ├── CITATION.cff ├── CHANGELOG.md ├── .github/workflows/ │   ├── CI.yml │   └── Benchmark.yml │ ├── src/ │   ├── Synchrosqueezing.jl │   ├── types.jl │   ├── validation.jl │   ├── grids.jl │   ├── padding.jl │   ├── normalization.jl │   │ │   ├── wavelets/ │   │   ├── abstract.jl │   │   ├── morlet.jl │   │   └── morse.jl │   │ │   ├── transforms/ │   │   ├── cwt.jl │   │   ├── stft.jl │   │   └── derivatives.jl │   │ │   ├── squeezing/ │   │   ├── phase.jl │   │   ├── reassignment.jl │   │   ├── ws st.jl │   │   ├── fsst.jl │   │   └── second_order.jl │   │ │   ├── ridges/ │   │   ├── costs.jl │   │   ├── dynamic_programming.jl │   │   ├── multiridge.jl │   │   ├── refinement.jl │   │   └── diagnostics.jl │   │ │   ├── reconstruction/ │   │   ├── inverse_cwt.jl │   │   ├── inverse_sst.jl │   │   ├── components.jl │   │   └── residual.jl │   │ │   └── performance/ │       ├── plans.jl │       ├── workspaces.jl │       ├── kernels.jl │       └── threading.jl │ ├── test/ │   ├── runtests.jl │   ├── test_wavelets.jl │   ├── test_cwt.jl │   ├── test_phase.jl │   ├── test_sst.jl │   ├── test_ridges.jl │   ├── test_inverse.jl │   ├── test_multicomponent.jl │   ├── test_second_order.jl │   └── test_python_reference.jl │ ├── reference/ │   └── python/ │       ├── requirements.txt │       ├── generate_reference.py │       ├── validate_reference.py │       └── README.md │ ├── benchmarks/ │   ├── julia/ │   ├── python/ │   ├── run_comparison.py │   └── results/ │ ├── examples/ │   ├── basic_chirp.jl │   ├── two_components.jl │   ├── crossing_chirps.jl │   ├── noisy_signal.jl │   └── bns_synthetic.jl │ └── docs/     ├── make.jl     └── src/         ├── theory.md         ├── conventions.md         ├── api.md         ├── validation.md         ├── performance.md         └── limitations.md

In the implementation, use `wsst.jl` without a space. Optional Python tooling and large generated datasets must not be required to install or load the Julia package.

### Recommended dependencies

| Dependency                                                                               | Responsibility                                         |
| ---------------------------------------------------------------------------------------- | ------------------------------------------------------ |
| [FFTW.jl](https://github.com/JuliaMath/FFTW.jl?utm_source=chatgpt.com)                   | FFTs, inverse FFTs, reusable FFT plans                 |
| AbstractFFTs.jl                                                                          | FFT interface where useful                             |
| LinearAlgebra                                                                            | Array and linear-algebra operations                    |
| Statistics                                                                               | Basic statistical utilities                            |
| [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl?utm_source=chatgpt.com) | Repeatable Julia benchmarks, development-only          |
| Test, Random                                                                             | Unit tests and reproducible synthetic signals          |
| Documenter.jl                                                                            | API and mathematical documentation, documentation-only |

Implement FFT-based CWT directly rather than using another Julia CWT package as the package's fundamental computational dependency. [ContinuousWavelets.jl](https://github.com/UCD4IDS/ContinuousWavelets.jl?utm_source=chatgpt.com) is nevertheless worth inspecting and optionally comparing against: it already implements Fourier-domain CWTs, including Morlet-type wavelets.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub

+1

Use package extensions for optional GPU or plotting dependencies. Do not make plotting libraries mandatory.

## 4. Public API and data structures

The API should be convenient for scientific users while exposing enough configuration for reproducibility.

### 4.1 Proposed usage

The following is the target API to implement, not code from an existing library.

```

using Synchrosqueezing

fs = 8192.0
t = (0:8191) ./ fs

# Evolving-frequency synthetic signal
phase = 2π .* (100 .* t .+ 200 .* t.^2)
x = cos.(phase)

# Configure an analytic wavelet
wavelet = MorletWavelet(ω0=6.0)

# Obtain SST and CWT coefficients
result = wsst(
    x;
    fs=fs,
    wavelet=wavelet,
    voices_per_octave=32,
    freqrange=(50.0, 1200.0),
    threshold=:adaptive,
    padding=:reflect,
    order=1
)

# Recover the strongest continuous ridge
ridges = extract_ridges(
    result;
    n_ridges=1,
    penalty=2.0,
    bandwidth=5,
    method=:dynamic_programming
)

# Recover instantaneous frequency
f = ridges.frequencies[:, 1]

# Reconstruct the corresponding component
component = reconstruct(
    result,
    ridges;
    bandwidth=8,
    mode=:real
)

# Reconstruct the full input signal
full_signal = inverse_wsst(result)

```

Additional core entry points:

```
cwt(x; kwargs...)
stft(x; kwargs...)
wsst(x; kwargs...)
fsst(x; kwargs...)

instantaneous_frequency(transform; kwargs...)
extract_ridges(transform; kwargs...)
reconstruct(transform, ridge; kwargs...)

inverse_cwt(transform)
inverse_wsst(transform)
inverse_fsst(transform)
```

Naming choices can evolve, but choose one consistent scheme and avoid competing synonyms such as `ssq_cwt`, `sst_cwt`, and `wsst` throughout the public API.

### 4.2 Types

Define concrete, parametric types for key numerical structures.

```
abstract type AbstractWavelet end
abstract type AbstractTFResult end
abstract type AbstractSSTResult <: AbstractTFResult end

struct MorletWavelet{T} <: AbstractWavelet
    ω0::T
end

struct MorseWavelet{T} <: AbstractWavelet
    β::T
    γ::T
end

struct SSTConfig{T,W}
    fs::T
    wavelet::W
    nvoices::Int
    freqrange::Tuple{T,T}
    threshold::T
end
```

These are indicative type definitions, not the final complete implementation.

Define additional types:

- `CWTResult`: coefficients, scales, time grid, metadata.
- `SSTResult`: squeezed coefficients, original transform, frequency grid, validity metadata, normalization.
- `RidgeResult`: integer indices, frequencies in Hz, ridge amplitudes, validity mask, confidence diagnostics.
- `SSTPlan`: reusable immutable or safely shareable transform configuration.
- `SSTWorkspace`: mutable scratch storage, designed for reuse.

Keep arrays in a consistent shape:

\\[ \mathrm{coefficients}[j,n] = T(f_j,t_n). \\]

That is, frequency or scale on the first dimension, time on the second.

Specify positive-frequency ordering as ascending by default. A row-order conversion should be applied explicitly when interfacing with Python, never hidden within numerical routines.

### 4.3 Metadata requirements

Every transform result should retain sufficient information to reconstruct its input without needing users to remember the original configuration.

Store or make retrievable:

- Sampling frequency and sampling interval.
- Original signal length and padded length.
- Original real or complex signal type.
- Wavelet family and parameters.
- Frequency or scale grid, including bin edges.
- Transform normalization.
- Reassignment method and threshold.
- Boundary padding mode and crop indices.
- Reconstruction coefficient or quadrature weights.
- Cone-of-influence or reliability mask where applicable.

Do not store large optional intermediates unconditionally. Provide `keep_cwt`, `keep_derivative`, and `keep_phase` options.

## 5. Implementation instructions by module

### Module A — Grid generation, FFT conventions, and padding

Implement `grids.jl`, `padding.jl`, and `normalization.jl` before writing the actual transform.

Required functions include:

```
angular_fftfreq(n, dt)
linear_frequency_grid(fmin, fmax, nbins)
log_frequency_grid(fmin, fmax, voices_per_octave)
scale_grid(wavelet, fs, n, freqrange, nv)
pad_signal(x, strategy, target_length)
crop_original(x, padding_metadata)
```

Requirements:

1. Work for odd and even array lengths.
2. Correctly construct negative FFT frequencies, DC, and the Nyquist bin.
3. Never confuse rad/s, Hz, and dimensionless wavelet scale.
4. Validate \\(0\<f\_{\min}\<f\_{\max}\<f_s/2\\) for positive-frequency analyses.
5. Support `:none`, `:zero`, `:reflect`, and `:periodic` padding.
6. Preserve metadata enabling exact removal of padding.
7. Support non-power-of-two signal lengths.
8. Generate frequency-bin edges consistently with the centers.
9. Add tests confirming monotonicity, spacing, bin coverage, and FFT sign conventions.

The agent must create a documentation page specifying these conventions before attempting Python numerical parity.

### Module B — Analytic wavelets

Implement analytic wavelets in Fourier space using functions such as:

```
wavelet_fourier(wavelet, ω)
wavelet_peak_frequency(wavelet)
wavelet_admissibility(wavelet)
wavelet_support(wavelet; tolerance=...)
```

Start with Morlet and then introduce Morse wavelets:

\\[ \widehat\psi\_{\beta,\gamma}(\omega) = C\_{\beta,\gamma} \omega^\beta e^{-\omega^\gamma} \Theta(\omega). \\]

For generalized Morse wavelets, do not confuse the wavelet shape parameter \\(\gamma\\) with the coefficient magnitude threshold used by SST.

Tasks:

1. Implement analytical Fourier-domain evaluation.
2. Ensure zero response at nonpositive frequencies under the selected analytic convention.
3. Implement normalization and admissibility integrals.
4. Calculate effective peak frequency and scale-frequency mapping.
5. Test wavelet localization on monochromatic inputs.
6. Compare generated wavelet responses with Python for matching parameters.

### Module C — FFT-based CWT

The default CWT kernel should:

1. Pad the input once.
2. Compute its FFT once.
3. Iterate over wavelet scales.
4. Evaluate wavelet responses on the FFT grid.
5. Multiply the signal spectrum by each wavelet response.
6. Compute the inverse FFT for that scale.
7. Store coefficients in the output matrix.
8. Crop to the original time domain.

Its complexity is approximately

\\[ \mathcal O(N_a N\log N), \\]

where \\(N_a\\) is the number of scales.

Support in-place execution:

```
cwt!(output, workspace, plan, signal)
```

Avoid allocating a wavelet spectrum, product array, and output array inside every scale iteration.

An important optimization is to calculate the transform derivative during the same pass:

```
cwt_with_derivative!(W, dW, workspace, plan, x)
```

This permits FFT-spectrum reuse and provides both arrays required for SST.

Tests must verify the transform independently before introducing synchrosqueezing.

### Module D — Phase transform

Implement:

```
phase_transform!(frequency, W, dW, config)
```

Mathematical requirements:

\\[ f=\frac{\operatorname{Im}(dW/W)}{2\pi}. \\]

Avoid performing a complex division when it can be implemented using a stable ratio:

\\[ \operatorname{Im}\left(\frac{dW}{W}\right) = \frac{\operatorname{Im}(dW\overline W)}{|W|^2}. \\]

For values above the validity threshold, this is an efficient computational form.

Define threshold policies:

- `:absolute`: fixed coefficient magnitude threshold.
- `:relative`: fraction of a reference amplitude.
- `:adaptive`: threshold estimated from a documented robust noise estimator.

The first release should implement absolute and relative thresholds; adaptive thresholding can follow after testing.

Return invalid frequencies through a separate mask, or avoid allocating a frequency matrix altogether when reassignment can be fused.

### Module E — Synchrosqueezing kernel

Implement:

```
synchrosqueeze!(T, W, dW, config)
```

Use the following baseline algorithm:

```

For each time index n:
    Initialize T[:, n] to zero

    For each wavelet scale k:
        Read W[k, n]
        Check coefficient threshold

        If valid:
            Compute instantaneous frequency
            Determine destination frequency bin
            Calculate quadrature contribution
            Add the complex contribution to T[bin, n]

Return squeezed transform

```

Important numerical details:

- Derive the quadrature weight for the implemented wavelet normalization.
- Ignore values outside the requested output frequency interval, and record lost contributions when diagnostics are enabled.
- Never use integer rounding on logarithmic frequency axes without using the actual frequency-grid mapping.
- Handle out-of-range, nonfinite, and near-zero coefficient values.
- Support nearest-bin and conservative linear-bin accumulation.
- Test conservation of reassigned coefficient sums for all accepted coefficients.
- Keep DC and negative-frequency handling explicit.

Performance requirement: the first implementation must have a correct scalar reference kernel. Only optimize after an independently verified comparison with that kernel.

Julia arrays are column-major. With the frequency dimension first, each time column is contiguous, so iterating over frequency bins inside time loops can be cache-friendly.

For parallel execution, assign separate time columns to threads. Each thread writes exclusively to its own output columns, avoiding atomics and data races.

### Module F — STFT and Fourier SST

Introduce STFT with a Gaussian analysis window and explicit control over:

```
nfft
window_length
hop_length
window_type
padding
```

Define a window-centered transform convention and test its sinusoidal phase response before implementing the associated reassignment formula.

Implement the frequency derivative using the appropriate derivative-window STFT, and then implement `fsst`.

The first FSST release can restrict hop size to one sample where necessary for faithful reconstruction. General hop lengths should only be enabled once the inversion and sampling implications are properly resolved.

Treat STFT-based reconstruction separately from wavelet-based reconstruction. STFT reconstruction may depend on window normalization and overlap-add conditions, while WSST uses wavelet admissibility and scale quadrature.

### Module G — Dynamic-programming ridge tracking

Implement at least three algorithms:

| Method                     | Purpose                                            |
| -------------------------- | -------------------------------------------------- |
| Maximum energy per time    | Baseline; no continuity                            |
| Dynamic programming        | Accurate global continuous ridge                   |
| Banded dynamic programming | Faster ridge tracking with bounded frequency jumps |

For a cost function \\(D(j,n)\\), the dynamic-programming recurrence is

\\[ F_n(j)=D(j,n)+ \min_i \left[ F\_{n-1}(i)+\lambda d(i,j)^2 \right]. \\]

Store a predecessor index for each state, select the optimal final state, and backtrack.

A naive implementation costs \\(\mathcal O(N_t N_f^2)\\) and is unacceptable for large frequency grids.

Implement the following progression:

G1. Reference algorithm

Use the direct recurrence to verify exact correctness against exhaustive solutions on very small matrices.

G2. Banded algorithm

Restrict candidate transitions to

\\[ |i-j|\leq B. \\]

Complexity reduces to approximately

\\[ \mathcal O(N_tN_fB). \\]

Make \\(B\\) configurable. Clearly document that this restricts allowable trajectories.

G3. Squared-distance transform

For uniformly spaced transformed frequency coordinates, implement an exact quadratic distance-transform minimum, reducing the transition step to linear complexity in the frequency dimension under the appropriate assumptions.

For arbitrary nonuniform grids, use the reference algorithm or a rigorously validated alternative.

G4. Frequency-informed tracking

Support a prior frequency interval, maximum jump rate, and optional prior curve:

```
extract_ridges(
    result;
    freqrange=(1500.0, 4000.0),
    max_rate=10000.0,
    prior=nothing
)
```

Use physical frequency units when applying physical slope constraints.

### Module H — Multiple ridge extraction

The initial approach should extract the strongest ridge, suppress its neighborhood, then extract subsequent ridges.

Implement:

```
extract_ridges(T; n_ridges=3, bandwidth=5)
```

However, sequential extraction may fail at crossings or when components approach each other.

Introduce an optional joint multi-ridge algorithm with:

- Exclusion or separation penalties.
- Explicit allowance for crossing trajectories.
- Optional mode-specific frequency priors.
- Ridge birth and termination states.
- Missing-observation penalties.
- A maximum allowable frequency-change rate.

These features must not be allowed to compromise the performance or clarity of the initial single-ridge algorithm.

### Module I — Ridge reconstruction and residuals

Implement:

```
reconstruct(result, ridge; bandwidth=...)
inverse_wsst(result)
reconstruct_components(result, ridges; bandwidths=...)
```

Support fixed bandwidth in frequency bins and physical Hz. Add adaptive bandwidth based on local spectral energy later.

For multiple components, define an explicit ownership policy for overlapping reconstruction bands. Otherwise the same SST coefficient can be assigned to two components and energy can be double-counted.

A sensible initial policy is non-overlapping nearest-ridge assignment within each ridge's permitted reconstruction region.

Provide residual diagnostics, including reconstruction error and the proportion of coefficients excluded from reconstruction.

### Module J — Second-order SST

Only begin after modules A–I are stable.

Derive the higher-order expressions in `docs/src/theory.md`, implement the necessary STFT window-moment transforms, and test the result against signals whose phase is exactly quadratic.

The strongest validation is that a linear chirp should exhibit substantially improved frequency concentration relative to first-order SST when its assumptions hold.

Use independent theoretical and literature-based references. Do not claim correctness merely because an implementation produces a visually thinner ridge.

## 6. Testing strategy: analytical correctness before Python parity

The package should have three independent levels of validation:

1. Analytical tests: known solutions with exact phases, frequencies, and amplitudes.
2. Cross-implementation tests: identical inputs and configurations in Julia and `ssqueezepy`.
3. Regression and performance tests: numerical consistency, runtime, memory use, and scalability.

Python comparison is essential, but it must not replace analytical testing. If both implementations share a normalization error or an inappropriate threshold, agreement alone does not establish correctness.

### 6.1 Analytical signal suite

Implement a deterministic signal-generation module under `test/signals.jl`.

| Test signal                   | Definition                                     | Required validation                  |
| ----------------------------- | ---------------------------------------------- | ------------------------------------ |
| Constant tone                 | \\(\cos(2\pi f_0t)\\)                          | Frequency, amplitude, reconstruction |
| Linear chirp                  | \\(\cos[2\pi(f_0t+\frac12kt^2)]\\)             | \\(f(t)=f_0+kt\\)                    |
| Quadratic-frequency chirp     | \\(\cos[2\pi(f_0t+\frac12kt^2+\frac13qt^3)]\\) | \\(f(t)=f_0+kt+qt^2\\)               |
| Exponential chirp             | \\(f(t)=f_0e^{\alpha t}\\)                     | Nonlinear instantaneous frequency    |
| AM tone                       | \\((1+\epsilon\cos\Omega t)\cos(2\pi f_0t)\\)  | Envelope reconstruction              |
| Two tones                     | \\(\cos(2\pi f_1t)+\cos(2\pi f_2t)\\)          | Two-ridge separation                 |
| Two moving chirps             | \\(f_1(t),f_2(t)\\)                            | Multiridge tracking                  |
| Crossing chirps               | \\(f_1(t_c)=f_2(t_c)\\)                        | Crossing ambiguity                   |
| Chirp with Gaussian noise     | Known chirp + seeded noise                     | Accuracy versus SNR                  |
| Abrupt phase/frequency change | Piecewise-defined phase                        | Transient behavior                   |
| Complex analytic tone         | \\(e^{i2\pi f_0t}\\)                           | Complex phase and amplitude          |
| Zero signal                   | \\(x(t)=0\\)                                   | Numerical robustness                 |

Use sample grids that avoid an endpoint duplication:

```
t = (0:N-1) ./ fs
```

Explicitly test odd and even lengths and frequencies near DC and the Nyquist frequency.

### 6.2 Frequency and reconstruction metrics

Implement reusable test utilities.

Instantaneous frequency RMSE:

\\[ E_f= \sqrt{ \frac{1}{N\_{\rm valid}} \sum\_{n\in\mathcal V} \left[ f\_{\rm estimated}(t_n)-f\_{\rm exact}(t_n) \right]^2 }. \\]

Relative reconstruction error:

\\[ E_x= \frac{\\|\widetilde x-x\\|\_2}{\\|x\\|\_2}. \\]

Signal-to-reconstruction error ratio:

\\[ \mathrm{SRER} = 10\log\_{10} \frac{\sum_n|x(t_n)|^2} {\sum_n|x(t_n)-\widetilde x(t_n)|^2}. \\]

Relative complex transform error:

\\[ E_T= \frac{\\|T\_{\rm Julia}-T\_{\rm Python}\\|\_F} {\\|T\_{\rm Python}\\|\_F}. \\]

Also implement phase error, amplitude error, ridge coverage, concentration statistics, and the fraction of points excluded by edge-effect masks.

All metrics must use clearly defined valid regions. Errors at heavily padded boundaries must be reported separately rather than silently discarded.

### 6.3 Required analytical unit tests

Validation checklist

0/13 completed

FFT and frequency-grid tests

Correct positive/negative frequency indexing, physical units, odd/even lengths and Nyquist handling

Wavelet tests

Correct Fourier response, normalization, admissibility and scale-to-frequency mapping

CWT linearity

W(ax+by) ≈ aW(x)+bW(y)

CWT derivative

Compare spectral derivative with the exact derivative for monochromatic signals

Pure-tone phase

Recover the exact tone frequency in regions of significant coefficient amplitude

SST coefficient conservation

Sum of deposited contributions equals the sum of valid weighted CWT coefficients, allowing floating-point rounding

Ridge optimization

Match brute-force enumeration on small frequency-time matrices

Ridge jump bounds

Every recovered ridge respects configured admissible transitions

Full inverse

Recover pure tones, chirps and synthetic mode mixtures

Component reconstruction

Recover known modes individually and check the residual

Numerical precision

Validate Float32 and Float64 paths independently

Memory and thread safety

No shared-state races; deterministic results under independently processed time columns

SST2 chirp test

Compare second-order instantaneous frequency with the analytical chirp law

&#x20;Copy test checklist

Initial numerical acceptance thresholds should be established from observed convergence and precision rather than chosen to guarantee success.

As provisional targets, use relative full-reconstruction errors below \\(10^{-3}\\) for well-resolved noiseless Float64 test signals in the reliable interior, and substantially tighter tolerances for isolated arithmetic identities. These are engineering targets, not mathematical guarantees; adjust them only with documented evidence.

## 7. Python comparison harness

The external reference suite should use `ssqueezepy`, with exact version pinning.

The Python implementation supports `ssq_cwt`, `ssq_stft`, `extract_ridges`, and `issq_cwt`, making it particularly suitable for testing the initial Julia implementation.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub

+2

### 7.1 Reference generation

Create `reference/python/generate_reference.py`.

It should:

1. Read test parameters from a JSON manifest.
2. Generate or load the deterministic test signal.
3. Run the Python CWT.
4. Run the Python synchrosqueezed CWT.
5. Extract one or more ridges.
6. Reconstruct the full signal and individual components.
7. Save numerical arrays with metadata.
8. Record the Python package version, commit identifier where available, array types, and configuration.

An illustrative Python reference script:

```
    wavelet="morlet",    fs=fs,    nv=32,    padtype="reflect",    astensor=False,    flipud=False,)ridge_indices = extract_ridges(    Tx,    freqs,    n_ridges=1,    penalty=2.0,    bw=5,    transform="stft",)recovered = issq_cwt(    Tx,    wavelet="morlet",)np.savez_compressed(    "reference_linear_chirp.npz",    signal=x,    time=t,    cwt=Wx,    sst=Tx,    frequencies=freqs,    scales=scales,    ridge_indices=ridge_indices,    reconstruction=recovered,)with open("reference_metadata.json", "w") as f:    json.dump(        {            "package": "ssqueezepy",            "version": ssqueezepy.__version__,            "fs": fs,            "N": N,            "wavelet": "morlet",            "nv": 32,            "padding": "reflect",        },        f,        indent=2,    )
```

This is an illustrative starting harness, not a substitute for checking the installed Python release. In particular, the use of `transform="stft"` above intentionally requests a linear-frequency ridge penalty for an already synchrosqueezed frequency grid. That choice should be verified against the intended algorithm; the string accepted by the Python API is not itself a guarantee that the resulting penalty is physically appropriate.

Before producing canonical reference files, verify the actual function signatures and defaults against the pinned package version.

### 7.2 Cross-language comparison methodology

A crucial distinction:

Algorithmic agreement and exact numerical agreement are different requirements.

`ssqueezepy` may use different wavelet normalization, scale generation, padding, phase thresholds, quadrature, reassignment binning, or frequency orientation.

Therefore, compare in stages:

| Stage | Quantity                 | What to establish                                |
| ----- | ------------------------ | ------------------------------------------------ |
| 1     | Wavelet response         | Same normalization and Fourier convention        |
| 2     | Scales and grids         | Identical physical coordinates                   |
| 3     | Raw CWT                  | Coefficient agreement                            |
| 4     | CWT derivative           | Phase-derivative agreement                       |
| 5     | Instantaneous frequency  | Consistent reassignment estimates                |
| 6     | Synchrosqueezed matrix   | Equivalent binning and complex values            |
| 7     | Ridge indices            | Identical objective and optimization constraints |
| 8     | Ridge frequencies        | Agreement in Hz                                  |
| 9     | Reconstructed components | Comparable amplitude and phase                   |
| 10    | Complete reconstruction  | Correct inverse normalization                    |

When exact transforms are not directly comparable, repeat tests using identical explicitly supplied scales, grids, thresholds, and wavelet parameters.

A very useful additional test is to compare the Julia ridge tracker directly against Python on an identical synthetic time-frequency matrix, bypassing the transforms entirely. That isolates differences in the ridge algorithms.

Likewise, compare the reassignment kernels with identical precomputed \\(W\\), \\(\partial_tW\\), and output grids.

### 7.3 Python references must remain outside the Julia runtime

Do not introduce PythonCall.jl or PyCall.jl as runtime dependencies.

For automated reference testing, either read generated fixtures or run an optional subprocess through the test harness. A normal user should be able to install, import, and use the Julia package without Python.

Version the small canonical reference cases in the repository. Large test matrices can be generated on demand and excluded from Git.

## 8. Performance design

The package is intended to outperform Python where possible through efficient kernels, FFT reuse, low allocation rates, memory layout, and parallel execution.

Do not optimize on the assumption that Julia is automatically faster.

### 8.1 Performance-critical operations

The most expensive parts are expected to include:

- Forward and inverse FFTs for every CWT scale.
- Creation and multiplication of wavelet frequency responses.
- Calculation of the CWT derivative.
- Synchrosqueezing deposition into the frequency-time matrix.
- Ridge tracking on large time-frequency matrices.
- Repeated full-transform computations during parameter studies.

Profile these operations separately.

### 8.2 Reusable plans and workspaces

Design reusable plans so repeated transforms do not reconstruct FFT plans and scratch buffers.

Target API:

```

plan = plan_wsst(
    8192;
    fs=8192.0,
    wavelet=MorletWavelet(ω0=6.0),
    voices_per_octave=32,
    freqrange=(100.0, 3000.0)
)

workspace = allocate_workspace(plan)

T = Matrix{ComplexF64}(
    undef, plan.nfreq, plan.noutput
)

wsst!(T, workspace, plan, x)

# Reuse the plan and workspace for another signal
wsst!(T, workspace, plan, x2)

```

The plan should contain immutable configuration, cached frequency grids and wavelet responses, and FFT planning information.

The workspace should contain mutable scratch buffers. Explicitly document whether a plan may be shared across threads and require separate workspaces whenever necessary.

### 8.3 Memory strategy

For \\(N_f\\) frequency bins and \\(N_t\\) samples, one `ComplexF64` matrix requires

\\[ M=16N_fN_t\quad\text{bytes}. \\]

Interactive memory estimator

Frequency bins256

Time samples65,536

Number of complex coefficient matrices

1

3

5

Estimated matrix memory (ComplexF64)

# 0.75 GiB

Excludes FFT workspaces, grids, masks, ridge backpointers, and other temporary allocations.

Requirements for reducing allocations:

1. Store data in type-stable arrays.
2. Support `Float32` and `Float64`.
3. Reuse FFT plans.
4. Use preallocated buffers and in-place kernels.
5. Fuse instantaneous-frequency calculation with reassignment where possible.
6. Avoid intermediate matrix broadcasts that allocate large arrays.
7. Reuse wavelet responses when sampling rate, signal length, and scales are unchanged.
8. Consider processing scales or output time segments in blocks, with mathematically correct handling of boundary support.
9. Allow users to request only the squeezed result, without retaining the original CWT or its derivative.
10. Benchmark reduced-memory modes separately from full-diagnostics modes.

### 8.4 Multi-threading

Implement CPU threading after the scalar implementation is validated.

Candidate decomposition:

- CWT: parallelize over independent scales, with per-thread FFT scratch buffers.
- Reassignment: parallelize over independent time columns.
- Ridge extraction: parallelize independently tracked signals or candidate processing where algorithmically appropriate.

Avoid nested oversubscription between Julia threads and FFTW threads. Benchmark the combinations `JULIA_NUM_THREADS=1`, `2`, `4`, `8`, and the corresponding FFTW thread settings.

Use a deterministic single-threaded mode as the reference for numerical regression.

### 8.5 GPU support

GPU support should be a later optional extension, not an obstacle to the CPU release.

Use a design compatible with `CUDA.jl`, and consider backend-generic kernels only after profiling justifies the added abstraction.

The Python reference already has optional GPU acceleration, so CPU-versus-GPU comparisons must be clearly distinguished from CPU-versus-CPU comparisons.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub

## 9. Benchmarking Julia against Python

Use [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl?utm_source=chatgpt.com) for Julia timings and an appropriate repeated wall-clock timer for Python. Measure execution after JIT compilation and initial setup. Also report cold-start and planning times separately.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub

### 9.1 Benchmark matrix

| Signal length | Frequency bins | Precision  | Algorithms      |
| ------------- | -------------- | ---------- | --------------- |
| 2,048         | 64, 128        | Float32/64 | CWT, SST        |
| 8,192         | 128, 256       | Float32/64 | CWT, SST, ridge |
| 16,384        | 256, 512       | Float32/64 | Full pipeline   |
| 65,536        | 256, 512       | Float32/64 | Full pipeline   |
| 131,072       | 512, 1024      | Float32/64 | Full pipeline   |

The benchmark configuration must also include the number of scales, not just output frequency bins.

Track the following separately:

- FFT-based CWT time.
- CWT derivative calculation time.
- Reassignment time.
- Dynamic-programming ridge time.
- Reconstruction time.
- End-to-end pipeline time.
- Peak memory usage.
- Allocated bytes.
- Precomputation time.
- Thread count.

Define the reported speedup as

\\[ S=\frac{t\_{\rm Python}}{t\_{\rm Julia}}. \\]

Report median runtime and variability across repeated trials.

### 9.2 Performance goals

Use these as aspirational targets rather than expected measured results.

| Benchmark                                       | Target                                       |
| ----------------------------------------------- | -------------------------------------------- |
| Single-thread Julia versus single-thread Python | At least competitive                         |
| Repeated CWT with reused plans                  | 1.5–3× speedup                               |
| Reassignment kernel                             | 2–5× speedup                                 |
| Ridge extraction                                | 2–5× speedup where algorithm permits         |
| End-to-end CPU pipeline                         | At least 2× speedup                          |
| Repeated workload allocations                   | Much lower than an allocating implementation |

Any failure to reach these targets should trigger profiling and investigation, not numerical shortcuts or changing configurations to make a comparison favorable.

Always publish actual results, even if Python is faster.

Also measure scaling with threads and signal length. Plot runtime versus \\(N\\), memory versus \\(N\\), and speedup versus \\(N\\). For algorithms with different asymptotic complexities, comparisons should span sufficient signal sizes to reveal the effect.

## 10. Scientific application: BNS post-merger signals

The initial package must remain general-purpose. Nevertheless, include a physically motivated benchmark derived from post-merger gravitational-wave signal structure.

Model a synthetic waveform as

\\[ h(t)=\sum\_{k=1}^{3} A_k(t) \cos\left[ 2\pi\int_0^t f_k(\tau)\\,d\tau+\phi\_{k,0} \right] +\epsilon(t). \\]

Here, the three components represent a dominant post-merger oscillation and two secondary spectral components, without assuming that every BNS simulation exhibits exactly three independently identifiable modes.

Implement synthetic tests with:

- A dominant mode near a few kHz.
- Secondary modes of weaker amplitude.
- Frequencies that evolve gradually and nonlinearly.
- Different damping times.
- Temporarily overlapping frequencies.
- Noise and early-time transients.

These should be configurable scenarios, not hard-coded universal BNS frequency relations.

### 10.1 GW-specific diagnostics

Beyond standard frequency RMSE, calculate the ridge tracking error relative to the known synthetic frequency curves.

For an extracted ridge \\(f_k^{\rm rec}(t)\\),

\\[ \Delta f_k(t)= f_k^{\rm rec}(t)-f_k^{\rm true}(t). \\]

Compare reconstruction amplitude and phase independently, since a mode may be reconstructed with the correct frequency but an incorrect amplitude.

For scientific use, preserve phase continuity and avoid smoothing that biases rapidly evolving physical frequencies.

### 10.2 Relevance to eventual symbolic regression

The eventual BNS analysis pipeline may look like:

Numerical BNS waveform

CoRe data or another numerical-relativity source

Synchrosqueezing.jl

Time-frequency localization

Ridge extraction and mode reconstruction

Instantaneous frequencies, amplitudes, phases

Quality control and physical classification

Mode assignment, uncertainty, simulation metadata

SymbolicRegression.jl

Empirical models for mode-frequency evolution

The package's responsibility ends at scientifically validated time-frequency estimates and reconstructed components. Classification and symbolic regression belong in the downstream BNS research project.

## 11. Development phases and acceptance gates

Use sequential milestones. Each phase should produce functioning code, unit tests, and documentation before proceeding.

Phase 0 — Repository and scientific conventions

Create the Julia package, CI, documentation structure, test infrastructure, and Python reference environment. Document normalization, FFT convention, frequency units, scale grids, and array orientations.

Gate: package loads, CI passes, and deterministic Python reference signals are generated.

Phase 1 — CWT and wavelets

Implement Morlet, FFT CWT, spectral derivative, scale generation, padding, inverse CWT, and analytical tests.

Gate: monochromatic and chirp tests pass; CWT matches Python under aligned conventions; inverse CWT works.

Phase 2 — First-order WSST

Implement the phase transform, thresholding, frequency grid, complex coefficient deposition, and inverse WSST.

Gate: exact tone-frequency localization, reassignment conservation, full reconstruction, and Python comparisons pass.

Phase 3 — Ridge extraction and component recovery

Implement maximum-energy, DP, banded DP, multiridge suppression, ridge reconstruction, and residual diagnostics.

Gate: known synthetic modes are tracked and reconstructed; ridge DP agrees with brute-force solutions.

Phase 4 — CPU optimization and Python benchmarks

Introduce reusable workspaces, cached wavelets, fused kernels, threading, advanced ridge optimizations, and detailed benchmark reports.

Gate: performance measurements are reproducible, optimized code agrees with the reference implementation, and no data races occur.

Phase 5 — STFT, FSST and second-order methods

Implement the STFT family and inverse, then the higher-order reassignment operators and rapidly varying chirp benchmarks.

Gate: first-order STFT tests pass, second-order chirp recovery is analytically validated, and reconstruction remains correct.

Phase 6 — Advanced multiridge and BNS validation

Add crossing-aware ridge tracking, uncertainty diagnostics, physically informed synthetic waveform tests, and optional GPU work.

Gate: known BNS-like synthetic modes can be identified with quantified errors, and package behavior at crossings is documented.

## 12. Agent execution rules

These instructions should govern the development process.

# Implementation contract for Synchrosqueezing.jl

You are developing `Synchrosqueezing.jl`, a standalone, high-performance scientific Julia package for synchrosqueezed time-frequency transforms, instantaneous-frequency estimation, multiridge extraction, and oscillatory component reconstruction.

The package must be mathematically correct, numerically tested, modular, and optimized for large-scale scientific workloads. It must support reproducible comparisons against the Python package `ssqueezepy`.

## Non-negotiable requirements

1. **Use Julia for the entire production runtime.** Python may be used only for reference-data generation and benchmark comparisons.
2. **Establish mathematical conventions before implementation.** Document Fourier transforms, wavelet normalization, scale-to-frequency conversion, sampling units, reassignment quadrature, and inversion constants.
3. **Implement independent analytical tests.** Agreement with Python is not sufficient proof of correctness.
4. **Build independently testable modules.** Wavelets, CWT, phase estimation, synchrosqueezing, ridge extraction, reconstruction, and performance utilities must be separated.
5. **Avoid speculative optimization.** Establish a correct reference implementation and benchmark before optimizing.
6. **Use idiomatic, type-stable Julia.** Support Float32 and Float64, parametric types, in-place operations, and reusable workspaces where appropriate.
7. **Never silently change numerical conventions.** All configurable mathematical choices must be explicit and recorded in metadata.
8. **Document limitations.** This includes time-frequency resolution, modal separation, mode crossing, endpoint effects, threshold sensitivity, and reconstruction bandwidth.
9. **Maintain test coverage.** Every new mathematical routine must have tests before it is integrated into higher-level algorithms.
10. **Do not falsify performance claims.** Benchmark identical mathematical workloads and report results even when Python is faster.

## Initial implementation order

Begin with repository creation, package registration metadata, source layout, CI, test infrastructure, and documentation.

Next implement frequency grids, FFT conventions, signal padding, Morlet wavelets, and FFT-based CWT.

Implement the CWT spectral time derivative and verify it against analytical signals.

Then implement first-order synchrosqueezing, complex coefficient reassignment, and inverse reconstruction.

Implement dynamic-programming ridge extraction, first for a single ridge and subsequently for multiple components.

Introduce component reconstruction and residual diagnostics.

Build the Python comparison harness early, but use it only after matching mathematical conventions.

Perform performance optimization after numerical validation.

Do not start second-order SST or GPU acceleration before the first-order CPU implementation is fully functioning and benchmarked.

## Required outputs at each phase

- Implemented Julia source files.
- Updated unit and analytical tests.
- Updated user and mathematical documentation.
- A summary of implemented functionality.
- Test results, including failures and unresolved discrepancies.
- Performance measurements where relevant.
- A concise description of outstanding scientific or technical limitations.

## Definition of completion

The first usable release is complete when a researcher can provide a uniformly sampled real or complex signal, compute a first-order WSST, extract continuous frequency ridges, reconstruct individual components, and quantify reconstruction errors without invoking Python.

All essential results must be independently tested on analytical signals, and comparison against a pinned `ssqueezepy` release must be reproducible.

The broader project is complete when FSST, second-order methods, multiridge tracking, comprehensive validation, and performance benchmarks meet their respective documented acceptance criteria.

**Start with Phase 0 and Phase 1. Do not attempt to implement the entire package in a single monolithic change.**

## 13. Primary references for the agent

These should be included in `docs/src/theory.md` and the package's scientific bibliography.

1. Daubechies, Lu & Wu (2011) — Synchrosqueezed wavelet transforms: An empirical mode decomposition-like tool. The mathematical basis of wavelet SST and its reconstruction properties.&#x20;

   [image](https://www.google.com/s2/favicons?domain=https://www.sciencedirect.com\&sz=32)

   ScienceDirect

2. Thakur et al. (2013) — The Synchrosqueezing algorithm for time-varying spectral analysis: robustness properties and new paleoclimate applications. Robustness and practical SST methodology.&#x20;

   [image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

   GitHub

3. Iatsenko, McClintock & Stefanovska — On the extraction of instantaneous frequencies from ridges in time-frequency representations of signals. Ridge extraction, direct versus path-based methods, and tracking considerations.&#x20;

   [image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

   GitHub

4. Oberlin, Meignen & Perrier — Second-order synchrosqueezing transform or invertible reassignment? Towards ideal time-frequency representations. Use this as a starting reference for deriving higher-order STFT synchrosqueezing.

5. [ssqueezepy source code](https://github.com/OverLordGoldDragon/ssqueezepy?utm_source=chatgpt.com) — Practical algorithms, implementation details, reference output, and profiling comparison.

6. [ContinuousWavelets.jl](https://github.com/UCD4IDS/ContinuousWavelets.jl?utm_source=chatgpt.com) — Julia-side reference for Fourier-domain wavelet implementation and numerical comparison.

## 14. Final architectural recommendations

The most consequential decisions are these:

First, prioritize WSST before FSST. This gives you a direct route to a usable CWT-based pipeline for nonstationary oscillations and an established Python reference to validate against.

Second, treat ridge extraction as a core algorithm rather than post-processing. For BNS waveforms, obtaining a narrow synchrosqueezed spectrum is only the first step; extracting a continuous, physically meaningful, reconstructible frequency trajectory is the actual scientific deliverable.

Third, distinguish ridge detection from mode identification. A maximum-energy ridge is not necessarily the same physical mode throughout a signal. This is particularly important if dominant and secondary modes approach or cross each other.

Fourth, make the implementation reusable at scale. For repeated waveform studies, cached transform plans, precomputed wavelets, in-place kernels, and low-allocation ridge tracking may be more valuable than optimizing an isolated one-time transform.

Finally, use the extracted modes as measurable scientific outputs.

The strongest initial milestone is therefore a scientifically validated, high-performance Julia WSST implementation with exact-convention Python comparisons, dynamic-programming ridge extraction, and quantitatively tested reconstruction. Once that foundation works, SST2 and more sophisticated mode tracking become controlled extensions rather than a source of intertwined numerical uncertainties.

# Revised architecture: Synchrosqueezing.jl using the Julia ecosystem

I would change the original development plan substantially: rather than implement the continuous wavelet transform, wavelet families, signal processing utilities, and FFT infrastructure ourselves, the package should build on Julia's existing scientific ecosystem.

The genuinely new functionality should be the synchrosqueezing algorithms, ridge extraction, and mode reconstruction, together with their integration and performance optimizations.

## 1. Packages identified

The two packages you mentioned are complementary, but they serve different purposes.

[Free High-Quality Julia Language Icon for Creative Design](https://images.openai.com/static-rsc-4/pVnhGG__O40YOf0yVBhe3pW_rAmzv2g20zvmEtLN9flK0BnY3kMFvED_0BpX6TEP9gP8_UdylgkSP64BiAtIxVzWAHq4G-3N1Jp565_JhoWw0HAHdrBAQHWWkDuRaU5fKwM8w69umSIDNuaI2Hvibt33oN98sABH_Lrswt6xHwc?purpose=inline)

[ContinuousWavelets.jl](https://github.com/UCD4IDS/ContinuousWavelets.jl?utm_source=chatgpt.com) — primary CWT backend

Strongest candidate

Provides Fourier-domain continuous wavelet transforms, several continuous wavelet families, scale configuration, boundary treatments, and inverse transforms. It already uses FFT-based convolution and supports complex wavelet coefficients.

This is the natural starting point for a wavelet-based SST implementation.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub

+1

[Free High-Quality Julia Language Icon for Creative Design](https://images.openai.com/static-rsc-4/pVnhGG__O40YOf0yVBhe3pW_rAmzv2g20zvmEtLN9flK0BnY3kMFvED_0BpX6TEP9gP8_UdylgkSP64BiAtIxVzWAHq4G-3N1Jp565_JhoWw0HAHdrBAQHWWkDuRaU5fKwM8w69umSIDNuaI2Hvibt33oN98sABH_Lrswt6xHwc?purpose=inline)

[Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl?utm_source=chatgpt.com) — discrete wavelets and preprocessing

Provides DWT, wavelet packets, filter banks, wavelet-based denoising, and noise-estimation utilities, including GPU support for several discrete transforms.

Useful for preprocessing, but it does not replace the continuous analytic transform required for conventional WSST.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub

[Free High-Quality Julia Language Icon for Creative Design](https://images.openai.com/static-rsc-4/pVnhGG__O40YOf0yVBhe3pW_rAmzv2g20zvmEtLN9flK0BnY3kMFvED_0BpX6TEP9gP8_UdylgkSP64BiAtIxVzWAHq4G-3N1Jp565_JhoWw0HAHdrBAQHWWkDuRaU5fKwM8w69umSIDNuaI2Hvibt33oN98sABH_Lrswt6xHwc?purpose=inline)

[ContinuousWavelet.jl](https://hmatuschek.github.io/ContinuousWavelet.jl/stable/?utm_source=chatgpt.com) — alternative CWT backend

This separate package provides planned wavelet transforms, precomputed analysis wavelets, wavelet synthesis, Morlet and Paul/Cauchy families, and coherence analysis.

Its planning model is interesting for processing many waveforms with the same configuration. It deserves a direct performance and accuracy comparison against `ContinuousWavelets.jl`.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://hmatuschek.github.io\&sz=32)

ContinuousWavelet

+1

### Broader package shortlist

| Package                                                                                           | Reuse in our project                                            | Recommendation               |
| ------------------------------------------------------------------------------------------------- | --------------------------------------------------------------- | ---------------------------- |
| [ContinuousWavelets.jl](https://github.com/UCD4IDS/ContinuousWavelets.jl?utm_source=chatgpt.com)  | CWT, analytic wavelet kernels, inverse transforms               | Core backend candidate       |
| [Wavelets.jl](https://github.com/JuliaDSP/Wavelets.jl?utm_source=chatgpt.com)                     | DWT, wavelet denoising, noise estimation                        | Optional preprocessing       |
| [DSP.jl](https://github.com/JuliaDSP/DSP.jl?utm_source=chatgpt.com)                               | Window functions, signal processing, STFT/spectrogram utilities | FSST infrastructure          |
| [FFTW.jl](https://github.com/JuliaMath/FFTW.jl?utm_source=chatgpt.com)                            | FFTs, inverse FFTs, planning and spectral derivatives           | Core                         |
| [StatsBase.jl](https://github.com/JuliaStats/StatsBase.jl?utm_source=chatgpt.com)                 | Statistical utilities and sampling                              | Optional                     |
| [Distances.jl](https://github.com/JuliaStats/Distances.jl?utm_source=chatgpt.com)                 | Distances and loss functions                                    | Optional; not a ridge solver |
| [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl?utm_source=chatgpt.com)          | Small fixed-size calculations                                   | Only if profiling justifies  |
| [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl?utm_source=chatgpt.com)  | CPU loop optimization                                           | Optional                     |
| [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl?utm_source=chatgpt.com) | Backend-independent parallel kernels                            | Future acceleration          |
| [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl?utm_source=chatgpt.com)                             | NVIDIA GPU computing                                            | Optional extension           |
| [NPZ.jl](https://github.com/JuliaIO/NPZ.jl?utm_source=chatgpt.com)                                | Reading Python-generated NumPy reference files                  | Test dependency              |
| [BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl?utm_source=chatgpt.com)          | Performance benchmarking                                        | Development dependency       |
| [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl?utm_source=chatgpt.com)                | API documentation                                               | Documentation dependency     |

Not every candidate should be a hard dependency. A relatively small package with well-designed extension points is preferable to one that loads a dozen unrelated numerical packages.

## 2. What we should reuse versus implement

| Functionality                       | Existing implementation                              | Our responsibility                               |
| ----------------------------------- | ---------------------------------------------------- | ------------------------------------------------ |
| Morlet wavelets                     | ContinuousWavelets.jl                                | Interface, normalization validation              |
| General continuous wavelet families | ContinuousWavelets.jl                                | Select SST-compatible analytic families          |
| CWT calculation                     | ContinuousWavelets.jl                                | Wrap backend, expose metadata                    |
| Inverse CWT                         | ContinuousWavelets.jl                                | Validate SST compatibility                       |
| Wavelet denoising                   | Wavelets.jl                                          | Optional preprocessing interface                 |
| FFT operations                      | FFTW\.jl                                             | FFT reuse for derivative calculations            |
| STFT windows                        | DSP.jl                                               | Implement only missing FSST-specific derivatives |
| Instantaneous frequency             | Not identified as a reusable complete implementation | Implement                                        |
| SST coefficient reassignment        | Not identified as a reusable complete implementation | Implement                                        |
| Second-order reassignment           | Not identified as reusable                           | Implement                                        |
| Continuous ridge extraction         | No verified suitable Julia package yet               | Implement a small dedicated solver               |
| Multiple-ridge tracking             | No verified suitable package yet                     | Implement                                        |
| Ridge-band reconstruction           | No verified suitable package yet                     | Implement                                        |
| Full SST inversion                  | Backend-dependent                                    | Implement using backend normalization            |
| Plotting                            | Existing plotting ecosystem                          | Provide optional plotting recipes                |

The important architectural decision is to wrap existing wavelet transforms and add SST-specific operations, rather than duplicating the transform implementation.

## 3. Important findings about ContinuousWavelets.jl

Its documentation confirms several features that materially change the implementation plan.

### 3.1 CWT and inversion already exist

The documented API includes:

```
using ContinuousWavelets

config = wavelet(
    Morlet(π);
    averagingType=NoAve(),
    β=2
)

W = cwt(x, config)
```

It also exposes frequency information through `getMeanFreq`, generates wavelet representations through `computeWavelets`, and provides `icwt` with inverse-frame choices including `NaiveDelta()`, `PenroseDelta()`, and `DualFrames()`.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://github.com\&sz=32)

GitHub

This means we should remove the following from our original implementation scope:

- Custom Morlet implementation as the default.
- Custom FFT-based CWT as the default.
- Reimplementation of scale generation where the backend already supports the required discretization.
- Custom inverse CWT, unless a specific missing capability is demonstrated.

However, `icwt` reconstructs from CWT coefficients. It is not automatically an inverse synchrosqueezed transform. We still need to derive and implement the inverse SST corresponding to the backend's normalization.

### 3.2 We need an SST adapter

`ContinuousWavelets.jl` returns coefficients in a documented time × scale layout. Our initial specification used frequency × time.

Rather than modify the backend, introduce an adapter that exposes a consistent view of its output.

```
abstract type AbstractCWTBackend end

struct ContinuousWaveletsBackend <: AbstractCWTBackend end
```

The adapter should expose:

```
compute_cwt(backend, x, config)
get_scales(backend, config, n)
get_frequencies(backend, config, n, fs)
get_normalization(backend, config)
inverse_cwt(backend, W, config)
```

These are proposed functions for our own package, not existing `ContinuousWavelets.jl` functions.

For frequency and scale metadata, reuse documented backend methods rather than independently rebuilding its grid. Validate whether returned frequencies are dimensionless, normalized to sampling, or already in Hz.

The adapter must also handle:

- Mother versus father-wavelet coefficients.
- Nonanalytic or approximately analytic wavelets.
- Boundary treatment.
- Frequency ordering.
- Quadrature weights.
- Complex signal compatibility.
- Transform coefficient layout.

For performance, avoid transposing full matrices unnecessarily. Use `PermutedDimsArray`, direct backend-compatible loops, or a suitable internal storage representation.

### 3.3 The main missing interface is the CWT time derivative

For first-order SST we require both

\\[ W(a,t) \quad\text{and}\quad \partial_t W(a,t). \\]

The existing library provides `cwt`, but I have not verified a public API that directly returns the exact derivative transform needed by SST.

The agent should investigate the implementation of `computeWavelets` and the internal Fourier-domain convolution, then establish whether the derivative can be computed by multiplying the existing wavelet response by \\(i\omega\\).

There are three implementation routes, in order of preference:

1. Public backend extension: reuse `ContinuousWavelets.jl`'s public Fourier-wavelet data or transform interfaces to calculate both \\(W\\) and \\(\partial_tW\\).
2. Upstream contribution: add a general-purpose derivative CWT or in-place transform API to `ContinuousWavelets.jl`, ideally contributing it to the upstream repository.
3. Local specialized derivative backend: implement only the missing derivative calculation using FFTW and the backend's exact wavelet kernels and normalization.

A finite-difference derivative of the CWT can serve as a validation baseline. It should not automatically become the default algorithm because its accuracy depends strongly on temporal resolution.

### 3.4 Check admissibility and analyticity

The `ContinuousWavelets.jl` documentation explicitly warns when a selected wavelet has significant undesirable frequency support. This matters directly for SST because phase-based reassignment assumes sufficiently well-behaved analytic coefficients.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://ucd4ids.github.io\&sz=32)

ContinuousWavelets.jl

The package should therefore validate SST suitability:

```
validate_sst_wavelet(config)
```

This should assess whether the mother wavelet is appropriate for positive-frequency reassignment, whether the frequency coverage is sufficient, and whether the selected reconstruction convention is compatible with the desired inverse SST.

Do not expose every `ContinuousWavelets.jl` wavelet as automatically SST-compatible.

## 4. DSP.jl eliminates much of the STFT work

I also confirmed that [DSP.jl implements a complex STFT](https://docs.juliadsp.org/dev/periodograms/?utm_source=chatgpt.com) directly, not merely spectrograms.

Its documented function is:

```
DSP.Periodograms.stft(
    signal,
    window_length,
    overlap;
    onesided=true,
    nfft=...,
    fs=...,
    window=...
)
```

It returns STFT coefficients. By contrast, `spectrogram` provides a power representation and is not the right input to synchrosqueezing because phase information is lost.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://docs.juliadsp.org\&sz=32)

periodogram estimation · DSP.jl

For FSST, reuse `DSP.jl` for window construction and STFT calculation, then add the required derivative-window transforms and phase reassignment.

One limitation needs testing: standard STFT interfaces often calculate coefficients at a subsampled set of time centers. For accurate reassignment and reconstruction, the desired hop length, centering, phase convention, and synthesis normalization must be compatible with the FSST formulation.

Do not assume the default `DSP.stft` settings will be suitable for inverse FSST.

## 5. Alternative ContinuousWavelet.jl backend

[ContinuousWavelet.jl](https://hmatuschek.github.io/ContinuousWavelet.jl/stable/?utm_source=chatgpt.com) is worth including in an early backend comparison, because its `ContinuousWaveletTransform` object precomputes analysis kernels for repeated transformations.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://hmatuschek.github.io\&sz=32)

ContinuousWavelet

That offers a potentially useful architecture for processing thousands of numerical-relativity waveforms with common sampling and transform parameters.

I would not make it a permanent dependency yet. First benchmark it against `ContinuousWavelets.jl` on the same waveforms, precision, wavelet families, and effective frequency resolutions.

The winner should be chosen based on four properties: mathematical suitability for SST, availability of transform-derivative information, inverse consistency, and repeated-transform throughput.

## 6. Revised repository structure

The updated repository should be notably smaller than the original proposal.

Synchrosqueezing.jl/

├── Project.toml ├── README.md ├── src/ │   ├── Synchrosqueezing.jl │   ├── types.jl │   ├── config.jl │   │ │   ├── backends/ │   │   ├── interface.jl │   │   ├── continuous_wavelets.jl │   │   ├── dsp_stft.jl │   │   └── derivative_support.jl │   │ │   ├── reassignment/ │   │   ├── phase.jl │   │   ├── first_order.jl │   │   ├── second_order.jl │   │   └── deposition.jl │   │ │   ├── ridges/ │   │   ├── costs.jl │   │   ├── dynamic_programming.jl │   │   ├── multiridge.jl │   │   └── refinement.jl │   │ │   ├── reconstruction/ │   │   ├── inverse_sst.jl │   │   ├── components.jl │   │   └── residual.jl │   │ │   └── performance/ │       ├── workspace.jl │       └── kernels.jl │ ├── ext/ │   ├── SynchrosqueezingWaveletsExt.jl │   └── SynchrosqueezingCUDAExt.jl │ ├── test/ │   ├── analytical/ │   ├── backends/ │   ├── reassignment/ │   ├── ridges/ │   ├── reconstruction/ │   └── python_reference/ │ ├── benchmarks/ │   ├── backends/ │   ├── python/ │   └── pipeline/ │ ├── examples/ │   ├── basic_sst.jl │   ├── chirp_reconstruction.jl │   └── bns_postmerger.jl │ └── docs/     ├── theory.md     ├── backends.md     ├── conventions.md     └── validation.md

The production code no longer needs general-purpose wavelet implementations, separate FFT-based CWT machinery, its own window library, or general DSP routines.

Instead, it consists primarily of adapters, reassignment operators, ridge solvers, reconstruction routines, and optimized kernels.

## 7. Revised agent instructions

The following should supersede the relevant implementation sections of the previous specification.

# Synchrosqueezing.jl — Ecosystem-first implementation directive

## Objective

Develop a standalone Julia package implementing first- and second-order synchrosqueezing, instantaneous-frequency estimation, ridge tracking, and oscillatory-component reconstruction.

The package must reuse mature Julia packages wherever appropriate. Do not duplicate continuous wavelet transforms, discrete wavelet transforms, FFTs, STFT routines, or digital signal processing utilities without demonstrating that existing implementations cannot meet our numerical or performance requirements.

The primary Python reference remains `ssqueezepy`.

## Task 1 — Conduct a dependency audit

Before implementing the SST algorithms, investigate:

- `UCD4IDS/ContinuousWavelets.jl`
- `JuliaDSP/Wavelets.jl`
- `hmatuschek/ContinuousWavelet.jl`
- `JuliaDSP/DSP.jl`
- `JuliaMath/FFTW.jl`
- `UCD4IDS/WaveletsExt.jl`

Also search for maintained Julia packages implementing phase reassignment, ridge tracking, second-order SST, or quadratic dynamic-programming distance transforms.

For each relevant package, inspect its source, documentation, project compatibility, testing status, exposed API, and availability in the Julia package registry.

Produce `docs/dependency_audit.md` documenting:

1. What functionality is provided.
2. Which documented public methods can be reused.
3. Whether the implementation handles complex-valued data.
4. What coefficient layout and normalization it uses.
5. Whether its frequency grids and wavelet responses are accessible.
6. Whether derivative transforms can be computed without duplicating the CWT.
7. Whether the library offers in-place transforms, plans, or workspaces.
8. Whether inverse transformations are available and numerically accurate.
9. What interfaces are missing.
10. Whether the package should be a required dependency, optional extension, development dependency, or rejected.

**Do not begin by implementing our own CWT.**

## Task 2 — Establish a baseline with ContinuousWavelets.jl

Use `ContinuousWavelets.jl` as the provisional primary backend.

Implement a minimal experimental script that:

- Constructs a Morlet wavelet.
- Computes the CWT of a constant-frequency sinusoid.
- Computes the CWT of a linear chirp.
- Obtains its frequency information.
- Inspects the transform coefficient layout.
- Reconstructs both inputs with `icwt`.
- Records reconstruction errors and elapsed times.

Test `NoAve()` and appropriate averaging configurations. Determine how the father-wavelet channel affects full reconstruction and whether it must be excluded from SST reassignment.

Verify the wavelet's frequency localization, approximate analyticity, and normalization.

## Task 3 — Investigate the CWT derivative

Inspect the source implementing the Fourier-domain wavelet convolution.

Determine whether the transform kernels can be accessed through public, supported APIs.

Investigate computing

```
∂W/∂t = IFFT(iω × FFT(x) × wavelet_response)
```

without rebuilding the wavelet definitions or scale grid.

Compare this spectral derivative against analytical pure-tone results and centered finite differences.

Prefer a small, well-scoped upstream enhancement to `ContinuousWavelets.jl` if its existing interfaces are insufficient. Do not depend on undocumented internal functions without explicitly isolating and version-testing them.

Deliver a numerical derivative benchmark and a written recommendation for the production implementation.

## Task 4 — Compare CWT backends

Benchmark `ContinuousWavelets.jl` and `ContinuousWavelet.jl` for repeated wavelet transformations.

Use equal signal lengths, equivalent wavelet parameters, matched frequency coverage, and comparable scale resolutions.

Test at least:

- 2,048 samples.
- 8,192 samples.
- 16,384 samples.
- 65,536 samples.

Report:

- Initial construction time.
- Warm transform time.
- Allocations and memory consumption.
- Reconstruction errors.
- Supported precision types.
- Complex-signal behavior.
- Access to derivative-supporting kernels.

Do not choose a backend solely by runtime. Mathematical suitability for SST and interface stability take precedence.

## Task 5 — Create the backend interface

Create `src/backends/interface.jl`.

Define an interface separating the SST algorithms from their input time-frequency transforms.

Required operations:

- Compute complex transform coefficients.
- Obtain physical frequency coordinates and/or wavelet scales.
- Obtain coefficient quadrature weights.
- Calculate the time derivative of the transform.
- Determine reconstruction normalization.
- Recover the complete signal using the backend's inverse where applicable.
- Obtain metadata about padding and edge reliability.

Use multiple dispatch and parametric types.

The public user-facing API should not expose backend-specific array conventions.

Avoid copying entire coefficient matrices merely to transpose them.

## Task 6 — Implement only the missing SST machinery

Implement:

- Stable phase-based instantaneous-frequency estimation.
- Magnitude thresholding and validity masks.
- First-order reassignment.
- Nearest-bin deposition.
- Conservative linear-frequency-bin deposition.
- Second-order reassignment in a later phase.
- Full SST reconstruction.
- Ridge-localized reconstruction.

All weighting and normalization must follow the selected backend's mathematical conventions.

Verify the redistribution of complex coefficients independently of ridge extraction.

## Task 7 — Investigate ridge solver reuse

Before implementing dynamic programming, search for existing Julia optimization or sequence-tracking packages supporting the required Viterbi-style recurrence.

The desired cost function combines local time-frequency energy with a frequency-continuity penalty.

Compare existing algorithms against:

- Straightforward dynamic programming.
- Banded dynamic programming.
- Exact quadratic-distance-transform optimization on regular grids.

If no package provides a suitable efficient solver, implement the ridge-specific recurrence directly.

Avoid introducing a large optimization framework merely for a small specialized recurrence.

Validate every optimized solver against exhaustive solutions on small test matrices.

## Task 8 — Reuse DSP.jl for STFT-based SST

Implement an STFT adapter based on `DSP.Periodograms.stft`.

Investigate how to obtain the derivative and window-moment transforms needed for first- and second-order FSST.

Verify window centering, phase conventions, temporal sampling, FFT normalization, and inverse compatibility.

Do not recreate a general-purpose STFT library.

## Task 9 — Keep Wavelets.jl optional

Use `Wavelets.jl` for optional denoising, thresholding, or preprocessing.

Do not confuse its discrete wavelet transform with the continuous wavelet representation used by standard WSST.

Keep it out of the critical SST runtime if the core algorithms do not require it.

## Task 10 — Preserve Python comparisons

Use `ssqueezepy` as the external validation reference.

First compare Julia's CWT backend with Python under matched wavelet definitions, scale grids, normalizations, and padding.

Then compare the instantaneous-frequency estimator, reassignment results, ridges, and reconstructions.

Where exact internal conventions cannot be matched, compare both implementations against independent analytical test signals.

Do not introduce Python as a runtime dependency.

## Task 11 — Optimize the new algorithms, not mature dependencies prematurely

Prioritize performance profiling of:

- Spectral CWT time derivatives.
- Complex coefficient reassignment.
- Ridge dynamic programming.
- Multiple-ridge extraction.
- Ridge-band reconstruction.

Only investigate replacing existing CWT or STFT kernels if profiling demonstrates that their implementation is the dominant bottleneck and their public interfaces prevent effective optimization.

Prefer improvements or contributions to upstream Julia libraries where possible.

## Immediate deliverables

Complete the following before implementing the production SST:

1. `docs/dependency_audit.md`
2. `experiments/continuouswavelets_baseline.jl`
3. `experiments/continuouswavelet_comparison.jl`
4. `experiments/cwt_derivative.jl`
5. `benchmarks/backend_comparison.jl`
6. `docs/backend_decision.md`

The backend decision document must explain which library was selected, the mathematical conventions it uses, how derivative transforms will be obtained, and the expected integration constraints.

After these deliverables, implement the minimum viable first-order WSST, followed by ridge extraction and reconstruction.

**Core principle: reuse existing numerical transformations; develop only the new algorithms and integration required for synchrosqueezing.**

## 8. My recommended first implementation

I would choose this dependency architecture:

Synchrosqueezing.jl

Public API, configuration, result types

ContinuousWavelets.jl

CWT / wavelets / inverse CWT

DSP.jl

STFT / windows / signal utilities

Our new numerical algorithms

Phase estimation · SST1/SST2 · Coefficient reassignment

Ridge extraction and reconstruction

Dynamic programming · Multiridge tracking · Mode recovery

Optional: Wavelets.jl preprocessing

FFTW\.jl acceleration

There is one particular research direction worth investigating early: adding native CWT derivative support upstream in `ContinuousWavelets.jl`. That would benefit this package and potentially other users developing reassignment, instantaneous-frequency, or higher-order time-frequency methods.

It would also let `Synchrosqueezing.jl` focus on its actual mathematical contribution rather than become another general wavelet-transform library.

The most useful next step for the agent is therefore not to start coding the SST directly, but to produce a tested `ContinuousWavelets.jl` integration prototype and determine precisely what functionality is missing.
