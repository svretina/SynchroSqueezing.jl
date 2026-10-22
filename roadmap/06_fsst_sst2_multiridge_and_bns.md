# Milestone 06 — Second-Order WSST and BNS Waveform Validation

## Objective

Extend the validated first-order `SynchroSqueezing.jl` implementation with second-order wavelet synchrosqueezing (WSST2), targeting oscillatory components with rapidly evolving instantaneous frequencies.

clone https://github.com/phamduonghung/WSST2
and use it as a reference to verify our results, also you can see how things 
are implemented there.

Preserve the existing CWT, reconstruction, ridge-tracking, amplitude-envelope, damping, and activity-detection infrastructure.

Evaluate whether WSST2 provides measurable scientific improvements over first-order WSST on synthetic chirps and numerical BNS post-merger waveforms.

STFT-based synchrosqueezing is explicitly excluded from this milestone.

## Phase A — Mathematical derivation

Study the second-order wavelet synchrosqueezing formulation of Pham and Meignen.

Derive all required operators using the exact L1-normalized CWT and Fourier conventions established in Milestone 03A.

Identify:

- The complex first-order frequency estimator.
- The wavelet-based time-reassignment operator.
- The chirp-rate estimator.
- The second-order instantaneous-frequency estimator.
- Required spectral wavelet derivatives and moment-weighted transforms.
- Conditions for numerical validity.
- Thresholds for ill-conditioned denominators.
- Fallback behavior when the second-order estimator is invalid.

Do not derive WSST2 by applying a second numerical derivative to the first-order frequency ridge.

Document all normalization constants, units, and sign conventions.

## Phase B — Reference implementation

Implement a transparent scalar/reference WSST2 estimator.

Keep second-order frequency estimation separate from complex coefficient deposition.

Reuse:

- Existing wavelet definitions.
- CWT backend and transform metadata.
- Frequency grids.
- Scale quadrature.
- Threshold policies.
- Nearest/linear deposition kernels.
- Inverse WSST routines.

Expose the corrected reassignment frequencies for diagnostic inspection.

Preserve the original complex CWT coefficients during reassignment.

## Phase C — Public API

Allow explicit selection of reassignment order:

    result1 = wsst(signal; order=1, ...)
    result2 = wsst(signal; order=2, ...)

The `order` parameter refers to synchrosqueezing order, not the order of a generalized Morse wavelet.

The result metadata must record the selected order, estimator thresholds, fallback policy, and wavelet configuration.

## Phase D — Analytical verification

Begin with isolated analytic signals:

1. Constant-frequency complex tone.
2. Constant-frequency damped tone.
3. Linear chirp with known chirp rate.
4. Quadratic-frequency chirp.
5. Exponentially evolving frequency.
6. Amplitude-modulated chirp.
7. Damped chirp.
8. Two well-separated chirps.

For linear chirps, verify that the second-order estimator reproduces the analytical instantaneous frequency under the assumptions of the selected formulation and numerical resolution.

Verify chirp-rate estimates independently.

Include positive and negative chirp rates.

Use multiple sampling rates, wavelet parameters, scale densities, Float32/Float64 precision, and signal lengths.

Report interior accuracy separately from boundary accuracy.

## Phase E — Numerical stability

Test:

- Near-zero CWT coefficients.
- Ill-conditioned time-reassignment derivatives.
- Abrupt amplitude changes.
- Weak components.
- Frequencies near the transform's reliable band limits.
- Finite-signal boundary effects.
- Strong damping.
- Close but still separable components.

Use documented validity thresholds and a clear fallback policy.

Do not allow NaN, Inf, or extreme ill-conditioned estimates to generate spurious frequency-bin assignments.

Measure the proportion of coefficients receiving valid second-order reassignment versus falling back to first-order reassignment or being excluded.

## Phase F — Reconstruction and scientific regression

WSST2 must preserve the weighted complex coefficient sum under complete valid deposition.

Reuse the standard inverse from Milestone 03A.

Validate:

- Full oscillatory reconstruction.
- Ridge-band component reconstruction.
- Complex amplitude envelopes.
- Instantaneous phase.
- Damping-time estimates.
- Activity-aware termination.
- Residual and ownership conservation.

Run the existing Milestone 03–05 scientific regression suite with both reassignment orders where applicable.

No regression of first-order WSST is permitted.

## Phase G — Quantitative SST1 versus SST2 comparison

Compare both methods under identical wavelet, scale, frequency-grid, and reconstruction configurations.

Measure:

- Instantaneous-frequency RMSE and bias.
- Ridge concentration and effective width.
- Component reconstruction error.
- Amplitude-envelope error.
- Damping-time error.
- Supported activity duration.
- Boundary error.
- Fraction of invalid reassignment estimates.
- Runtime and allocations.

Do not conclude that WSST2 is superior merely because its spectrogram appears sharper.

Identify a quantitative chirp-rate regime where second-order reassignment provides meaningful improvement.

Also identify regimes where first-order WSST is equal or better.

## Phase H — Independent reference validation

Use the published wavelet-based WSST2 implementation by Pham and Meignen.

Study its operator conventions, wavelet normalization, derivative formulas, and validity conditions.

Compare against its analytical examples and reproduce selected numerical experiments where practical.

Account explicitly for differences in frequency units, scale definitions, padding, and Fourier conventions.

Keep external code as a validation dependency rather than a Julia runtime dependency.

## Phase I — BNS post-merger validation

Apply both first- and second-order WSST to representative numerical BNS post-merger waveforms.

Begin with well-separated frequency features and user-specified frequency-search intervals.

For each selected mode, recover:

- Instantaneous frequency trajectory.
- Complex reconstructed component.
- Amplitude envelope.
- Reliable support interval.
- Damping-time or local damping-rate estimates.
- Reconstruction uncertainty and parameter sensitivity.

Check wavelet bandwidth, scale density, sampling, reconstruction-band width, and boundary sensitivity.

Use available published empirical frequency relations to define plausible initial search intervals, but do not impose these relations as exact mode trajectories.

Treat physical identification of spectral features as a separate scientific inference problem.

## Phase J — Performance

Benchmark WSST1 and WSST2 separately.

Report:

- CWT and derivative calculations.
- Additional second-order operators.
- Frequency reassignment.
- End-to-end execution.
- Allocations and peak working memory.
- Single-thread and supported multithread performance.

Reuse plans, backend wavelet information, and workspaces where appropriate.

Do not optimize second-order kernels until the scalar implementation passes analytical validation.

## Acceptance gate

Milestone 06 is complete when:

1. WSST2 is mathematically derived under the established CWT conventions.
2. An independent scalar reference implementation exists.
3. Linear-chirp and chirp-rate analytical tests pass.
4. Ill-conditioned coefficients are handled explicitly.
5. Complex coefficient conservation and standard inverse reconstruction remain correct.
6. First-order WSST numerical behavior is unchanged.
7. Component reconstruction and damping extraction work with second-order reassignment.
8. Quantitative comparisons identify where WSST2 helps and where it does not.
9. Results are compared against an independent published WSST2 reference.
10. Representative BNS waveform analyses are documented with uncertainties and limitations.
11. Performance costs are measured and reported.

## Explicit non-goals

- STFT-based FSST or FSST2.
- GPU implementation.
- General-purpose joint multiridge optimization.
- Automatic physical identification of every BNS spectral feature.
- Detector-noise parameter estimation.
- New empirical BNS waveform models or symbolic regression.

Those may be considered in later research milestones if supported by results.

# mark the stft backend as a todo in the repo readme.