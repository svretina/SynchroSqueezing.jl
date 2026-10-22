# Performance and scientific regression

Milestone 05 treats numerical behavior as part of the performance contract.
An implementation is not accepted merely because its WSST coefficient matrix
is close to a previous result: ridge choices, ownership, activity decisions,
termination, and damping estimates are checked as well.

## Frozen reference

`validation/scientific_reference_05a.json` records a deterministic two-mode
case and its environment, transform conventions, configuration, coefficient
summaries, ridge paths, masks, reconstructions, activity decisions, damping
times, local decay rates, and ownership conservation. The executable fixture
in `test/fixtures/scientific_reference_05a.jl` is checked by the full test
suite.

The regression compares the fused WSST path with the readable scalar phase
and deposition implementation. It also compares the constrained ridge solver
with the dense dynamic-programming reference. Continuous summaries use
recorded tolerances; downstream discrete paths and masks must match exactly.
A changed discrete decision is therefore visible even when aggregate matrix
errors remain small.

The reference signal is deterministic and has zero added noise. A random seed
is nevertheless recorded so noise can be enabled without making a new,
underspecified fixture.

## Initial profile and first optimization

An initial 1024-sample, 256-bin profile identified first-order ridge dynamic
programming as the dominant stage:

| Stage | Median time | Allocated memory |
|---|---:|---:|
| CWT and spectral derivative | 3.29 ms | 11.90 MiB |
| Phase estimation | 4.03 ms | 0.86 MiB |
| WSST after CWT | 6.34 ms | 4.10 MiB |
| Single constrained ridge, before | 521.24 ms | 6.92 MiB |
| Single constrained ridge, after | 19.56 ms | 6.18 MiB |

These measurements used Julia 1.13.0-rc1, one Julia thread, one FFTW thread,
and the same runtime-owned input. They are evidence for optimization order,
not portable performance promises.

The optimized solver preserves the same recurrence and tie-breaking. It
limits loops to each column's eligible row span, caches the time-independent
continuity objective under a 64 MiB memory cap, avoids allocating full index
lists for boundary diagnostics, and retains the dense solver as an independent
audit path. Existing exhaustive tests and the frozen end-to-end reference
require identical ridge indices and total cost.

## Reuse and batch analysis

`CWTPlan` is an immutable reusable plan containing the backend configuration
and wavelet bank. Construct it once with `plan_cwt`, then call `cwt(plan,
signal)` for compatible signals. A single `WSSTResult` should be passed to
`extract_ridges`; all requested frequency ranges then reuse the same CWT,
derivative, WSST coefficients, time metadata, and reconstruction convention.
`reconstruct_components` performs batched ownership bookkeeping for those
ridges.

The plan may be shared as read-only data. Current transform and analysis calls
allocate their own result arrays, so independent calls do not share a mutable
package workspace. Arrays stored inside plans and result objects remain Julia
arrays and must not be mutated concurrently by user code. Explicit reusable
mutable workspaces and threaded kernels are later Milestone 05 work; no
thread-safety or scaling claim is made for them yet.

## Reproducible profiling

Run the isolated full-pipeline benchmark environment:

```sh
julia --project=benchmarks -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
julia --project=benchmarks benchmarks/full_pipeline.jl
```

The generated JSON includes planning and first-call costs separately from
warmed BenchmarkTools trials. Every stage reports minimum, quartiles, median,
maximum, sample count, allocation count, and allocated bytes. It also records
Julia, package, Julia-thread, and FFTW-thread versions/configuration.

Python end-to-end speedup is intentionally not reported yet. The existing
`ssqueezepy` reference validates transform behavior, but it does not perform
the same ridge, ownership, activity, and damping workflow. A speed ratio would
therefore describe different scientific workloads.
