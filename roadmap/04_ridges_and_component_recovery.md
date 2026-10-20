# Milestone 04 — Ridges and component recovery

## Objective

Turn a valid WSST into continuous, reconstructible modes with transparent limits at crossings and overlaps.

## Scope

- Maximum-energy-per-time baseline.
- Single-ridge dynamic programming with local energy cost and a physical- or log-frequency continuity penalty.
- A direct reference recurrence validated against exhaustive small cases, followed by configurable banded DP and, only on suitable regular transformed grids, an exact quadratic-distance-transform optimization.
- Frequency priors, frequency ranges, maximum physical jump rate, and optional prior curves.
- Sequential multi-ridge extraction by suppression/exclusion bands; document its failure modes at crossings.
- Ridge-band reconstruction with fixed bins or Hz bandwidth, full reconstruction, residual reconstruction, and overlap ownership rules that prevent double counting.

## Validation

Test two tones, separated moving chirps, crossing chirps, noisy signals, abrupt changes, and BNS-like synthetic mixtures. Confirm that DP matches brute force on small matrices, respects allowed jumps, and is deterministic. Test ridge recovery separately on precomputed time-frequency matrices to isolate solver behavior from transform differences.

Measure instantaneous-frequency RMSE, phase/amplitude error, ridge coverage, reconstruction error, SRER, residual energy, edge-mask exclusions, and coefficient-exclusion proportions. Treat crossings and insufficient modal separation as diagnostics and documented ambiguity, never as silently correct tracking.

## Acceptance gate

Known synthetic modes are tracked and reconstructed with quantified errors; DP agrees with brute-force reference cases; residual/ownership behavior is defined; the package clearly reports limitations near crossings and boundaries.

## Explicit non-goals

- Do not add a broad optimization framework unless it demonstrably supplies the required Viterbi-style recurrence.
- No second-order reassignment or GPU work.
