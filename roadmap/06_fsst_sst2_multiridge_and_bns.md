# Milestone 06 — FSST, SST2, advanced multiridge, and BNS validation

## Objective

Extend a validated WSST pipeline in controlled stages: reuse `DSP.jl` for STFT infrastructure, add higher-order reassignment only from a consistent derivation, and validate difficult multiridge/BNS-like cases with quantified uncertainty.

## Scope

- An STFT adapter based on `DSP.Periodograms.stft`, with explicitly tested window centering, phase convention, time centers, hop size, FFT normalization, derivative-window transforms, and inverse compatibility.
- First-order FSST, initially allowing hop size one where required for faithful reconstruction.
- Second-order SST only after deriving all window-moment/scale-time operators from one convention and documenting them in the theory documentation. Validate against exact quadratic-phase chirps; never approximate SST2 as a naive second numerical derivative.
- Crossing-aware joint multiridge strategies: separation penalties, intentional crossing policy, birth/death states, missing-observation penalties, max-rate constraints, and optional mode priors.
- BNS-like synthetic scenarios with damped multi-kHz modes, nonlinear frequency evolution, temporary overlap, transients, and noise; do not encode them as universal physical relations.
- Optional GPU/`CUDA.jl` or backend-generic acceleration only after a profile justifies it, preferably as package extensions.

## Validation and documentation

Show first-order STFT inversion before SST2. Demonstrate improved concentration or frequency recovery for linear chirps under stated assumptions, using analytical metrics rather than appearance. Evaluate phase continuity, amplitude recovery, ridge error, uncertainty/confidence diagnostics, crossing behavior, and reliable-region masks for BNS-like cases. Include limitations involving resolution, threshold sensitivity, endpoints, overlap, and mode identification.

## Acceptance gate

FSST reconstruction passes its own conventions-based tests; SST2 has analytical chirp validation and documented stability conditions; advanced multiridge/BNS tests report quantified outcomes rather than unsupported identification claims. Optional GPU results, if any, are clearly separated from CPU comparisons.

## Completion definition

The first usable release is already achieved after Milestones 00–05 when users can compute a first-order WSST, extract continuous ridges, reconstruct components, and quantify error entirely in Julia. This milestone advances the broader research package rather than blocking that release.
