# Predictive ridge tracking

Predictive tracking extends the frequency-constrained first-order tracker with
a locally constant chirp-rate model. It is intended for modes whose
instantaneous frequency is smoother than competing time-frequency structures.
The original first-order DP remains the reference when curvature assumptions
are inappropriate.

## Strategies

| Strategy | API | Guarantee | Typical use |
|:--|:--|:--|:--|
| First-order DP | `extract_ridge(...; method=:dp)` | Global optimum of the first-order objective | Default; abrupt or weakly modeled motion |
| Exact second-order DP | `extract_predictive_ridge(...; method=:exact)` | Global optimum of the second-order objective | Offline validation and moderate constrained ranges |
| Prediction-window pruning | `extract_predictive_ridge(...; method=:pruned)` | Global optimum only over retained transitions | Faster smooth-mode tracking with an explicit approximation |

All strategies enforce the inclusive `freqrange`. A strong coefficient outside
the supplied interval can never attract the tracker.

## Constant-velocity predictor

For three consecutive frequencies, the physical slopes are

```math
v_{n-1}=\frac{f_{n-1}-f_{n-2}}{t_{n-1}-t_{n-2}},\qquad
v_n=\frac{f_n-f_{n-1}}{t_n-t_{n-1}}.
```

The predicted next frequency is

```math
f_{\mathrm{pred},n}=f_{n-1}+v_{n-1}(t_n-t_{n-1}).
```

This definition works for nonuniform time coordinates. A linear chirp has
`v_n-v_{n-1}=0`, including on a nonuniform time grid when represented exactly.

## Second-order objective

The exact solver minimizes

```math
\sum_n D(r_n,n)
+\lambda_f\sum_{n=2}^{N_t}
  \left(\frac{f_{r_n}-f_{r_{n-1}}}{s_f}\right)^2
+\lambda_c\sum_{n=3}^{N_t}
  \left(\frac{v_n-v_{n-1}}{s_c}\right)^2.
```

`D` is the same per-column normalized negative-log energy used by the
first-order tracker. Parameter units are explicit:

- `frequency_scale`: Hz per sample transition.
- `curvature_scale`: Hz/s, scaling the change between consecutive physical
  frequency slopes.
- `max_rate`: optional hard frequency rate in Hz/s.
- `max_chirp_rate_change`: optional hard bound on `abs(v_n-v_{n-1})` in Hz/s.
- `prediction_bandwidth`: fixed half-width in Hz around the predicted next
  frequency, used only by the pruned method.

The curvature penalty, hard chirp-rate-change constraint, and prediction
window are independent. A penalty changes the objective; a hard physical
constraint changes the admissible problem; predictive pruning is a heuristic
candidate restriction.

## Exact pair-state DP

The exact state is `(i, k)`, the two preceding frequency indices. A transition
`(i, k) -> (k, j)` adds the observation at `j`, the `k -> j` first-order cost,
and the physical slope-change cost. The first sample retains every eligible
start state, the second uses only the ordinary first-order transition, and
curvature begins at the third sample.

The unbanded runtime is ``O(N_t M^3)`` and predecessor storage is
``O(N_t M^2)``, where `M` is the number of bins inside `freqrange`. A finite
`max_rate` bands the candidate `j` states and can reduce runtime substantially.
This cost is why narrow, physically motivated search intervals remain
important.

Ties are deterministic: the lowest-frequency predecessor and lexicographically
lowest terminal pair win. Tests compare the complete path and cost with
independent exhaustive enumeration.

## Prediction-window pruning

```julia
ridge = extract_predictive_ridge(
    result;
    freqrange = (35.0, 65.0),
    method = :pruned,
    prediction_bandwidth = 2.0,
    curvature_penalty = 1.0,
    curvature_scale = 20.0,
)
```

The pruned solver keeps only predecessors whose constant-velocity prediction
falls within `prediction_bandwidth` of a candidate. It constructs the
predecessor interval directly on the sorted physical-frequency grid, reducing
work as the window narrows. The bandwidth remains fixed unless the user
changes it; it is not an estimated covariance or confidence interval.

Pruning can exclude the optimum of the unpruned second-order objective. A wide
window should be compared against `method=:exact` when establishing a safe
setting. A narrow window can reject genuine abrupt slope changes or make the
problem infeasible.

## Startup and priors

Optional `initial_frequency` and `initial_slope` values add soft quadratic
costs controlled by separate penalties and scales. They never expand
`freqrange` or override hard constraints. Incorrect strong priors can bias the
path and should be included in sensitivity tests.

One- and two-sample inputs use the first-order dense solver because curvature
cannot be identified. The returned `RidgeResult.config` records this fallback
as a `RidgeConfig` rather than claiming predictive inference.

## Choosing a strategy

Use first-order DP when abrupt physical changes are plausible, the signal is
short, or there is no defensible chirp-rate scale. Use exact second-order DP to
validate a smooth-mode model and tune curvature parameters on a manageable
frequency interval. Use prediction pruning only after comparing it with exact
DP on representative data.

A smoother trajectory is not automatically a more accurate trajectory.
Compare RMSE, coverage, boundary contacts, reconstruction error, and runtime.
At crossings, none of these independent single-ridge strategies guarantees
physical mode identity.
