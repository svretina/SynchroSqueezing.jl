# Milestone 05 — CPU optimization and reproducible benchmarks

## Objective

Optimize only validated work: derivative evaluation, complex deposition, ridge DP, multiridge extraction, and ridge-band reconstruction. Mature backend transforms are replaced only if profiling shows a material, unresolvable bottleneck.

## Scope

- Reusable immutable plans/configuration and mutable workspaces, with documented sharing/thread-safety rules.
- Cached backend wavelet information where supported, preallocated arrays, type-stable kernels, and reduced-diagnostic modes.
- Controlled fusion of phase estimation and deposition after scalar-reference equivalence tests.
- CPU threading only after deterministic single-thread results are established: independent time columns for deposition, and backend-supported safe decomposition for CWT work.
- `BenchmarkTools.jl` timing after JIT/warm-up; separately report cold start, planning, allocations, peak memory, and each pipeline stage.

## Benchmark protocol

Cover signal lengths 2,048 through 131,072, appropriate scale and output-bin counts, Float32/Float64, and CWT, WSST, ridge, and full-pipeline workloads. Compare equal mathematical configurations with pinned Python reference timings. Report median and spread, thread count, FFTW-thread configuration, memory, allocations, and speedup `t_python / t_julia`.

Test one, two, four, and eight Julia threads where available, while avoiding FFTW/Julia oversubscription. Publish results even when Python is faster. Keep the production path deterministic when run single-threaded and verify that threaded runs contain no shared-state races.

## Acceptance gate

Benchmark reports are reproducible, optimized results remain within established numerical tolerances of the scalar implementation, plans/workspaces demonstrate allocation reductions on repeated workloads, and thread-safety tests pass.

## Explicit non-goals

- Performance targets are aspirational, not evidence of correctness.
- No GPU requirement and no changing of configurations to favor a language comparison.
