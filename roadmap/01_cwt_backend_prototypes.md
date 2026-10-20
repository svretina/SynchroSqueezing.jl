# Milestone 01 — CWT backend prototypes and selection

## Objective

Establish a measured, mathematically suitable CWT backend. `ContinuousWavelets.jl` is the provisional primary candidate; `ContinuousWavelet.jl` is a planned-transform alternative worth comparing for repeated analyses.

## Deliverables

- `experiments/continuouswavelets_baseline.jl`:
  - construct a Morlet wavelet;
  - transform a pure sinusoid and a linear chirp;
  - inspect coefficient shape/order, scales, frequencies, padding, and analytic behavior;
  - test `NoAve()` and appropriate averaging choices;
  - reconstruct with `icwt` and record error and elapsed time.
- `experiments/continuouswavelet_comparison.jl` and `benchmarks/backend_comparison.jl` comparing both backends at 2,048, 8,192, 16,384, and 65,536 samples.
- `docs/backend_decision.md`, choosing a backend based first on SST suitability and interface stability, then inverse consistency and repeated-transform throughput.

## Required measurements

Use comparable wavelet parameters, effective frequency coverage, scale resolution, precision, and signal inputs. Report construction time, warm transform time, allocations, memory, inverse reconstruction error, Float32/Float64 behavior, complex-signal behavior, coefficient layout, frequency units, and access to derivative-supporting kernels.

Determine whether a father-wavelet channel exists, how it participates in full reconstruction, and whether it must be excluded from SST reassignment. Validate that the selected analysis wavelet is sufficiently analytic and admissible for positive-frequency phase reassignment.

## Acceptance gate

`docs/backend_decision.md` names the provisional production CWT backend, documents its conventions and integration constraints, and explains why the alternative is not selected yet. The chosen backend reconstructs the baseline signals within empirically justified interior-region tolerances.

## Explicit non-goals

- Do not reimplement the CWT, custom Morlet kernels, general scale grids, or inverse CWT merely to duplicate backend features.
- Do not yet claim Python coefficient parity; matching definitions comes after the backend contract is known.
