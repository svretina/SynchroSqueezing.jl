# Frequency-constrained ridge tracking

Milestone 04A treats prior frequency knowledge as a hard scientific
constraint. A ridge search must supply an inclusive `freqrange`; bins outside
that interval cannot enter the path:

```julia
ridge = extract_ridge(
    result;
    freqrange = (35.0, 65.0),
    method = :dp,
    continuity = :physical,
    penalty = 2.0,
    continuity_scale = 2.0,
    max_rate = 40.0,
)
```

`freqrange`, `max_rate`, and reconstruction bandwidth are independent:

- `freqrange` is the inclusive region in Hz where the tracker may search.
- `max_rate` is an optional hard transition limit in Hz/s.
- `bandwidth_hz` is the half-width used later to collect coefficients.

The package never expands a range to reach a nearby grid point. It rejects a
range containing no grid bins.

## Objective and solvers

For eligible bin ``j`` and time ``n``, the observation cost is

```math
D(j,n)=-\log\left(\frac{|T_{j,n}|^2}{\max_{k\in\mathcal F_n}|T_{k,n}|^2}
                  +\epsilon\right).
```

A path minimizes the sum of observation costs and

```math
\lambda\left(\frac{q(f_j)-q(f_k)}{s}\right)^2,
```

where `q(f)=f` for `continuity=:physical` and `q(f)=log(f)` for
`continuity=:log`. `continuity_scale` is therefore measured in Hz per sample
transition or log-Hz per sample transition. This fixed-step convention is
separate from `max_rate`, which uses the actual time difference and Hz/s.

The methods are:

- `:maximum`: independently choose the strongest eligible bin in each column.
- `:dense`: evaluate every admissible predecessor; this is the correctness
  reference.
- `:dp`: restrict predecessor iteration with the physical `max_rate` band and
  is the production default. Without a rate constraint it solves the same
  dense transition set.

Both DP implementations return the global optimum of this discrete objective.
Equal costs choose the lower-frequency predecessor and endpoint.

## Activity, boundaries, and priors

Optimization returns a path even in a zero- or low-energy interval. The
`RidgeResult` therefore keeps `valid`, `active`, and
`at_frequency_boundary` separate. `activity_threshold` is a squared-magnitude
fraction of the largest finite coefficient in the complete input matrix; it
does not alter the optimized path. Inactive samples are excluded from
component masks and from frequency-error metrics.

`valid_mask` can exclude unreliable coefficient locations. The optional
`frequency_envelope` accepts one interval per time sample, or a function of
time, and intersects that interval with the mandatory global `freqrange`.
`initial_prior` and `end_prior` add explicit per-bin endpoint costs.

A path touching `at_frequency_boundary` should be interpreted carefully: the
physical mode may be leaving an overly narrow supplied range. This flag is not
proof of failure, but it prevents silent interpretation of a clipped path as
an unconstrained measurement.

## Multiple known modes

Track known modes independently and preserve their supplied order:

```julia
ridges = extract_ridges(
    result;
    freqranges = [(35.0, 65.0), (90.0, 125.0)],
)
```

Pass one `RidgeConfig` or named tuple per range through `configs` when modes
need different penalties or rates. Overlapping trajectories generate a
warning and can be inspected with `ridge_overlap`. Independent tracking does
not establish physical mode identity at crossings.

## Reconstruction and ownership

The search interval does not determine the reconstruction band:

```julia
component = reconstruct_component(result, ridge; bandwidth_hz = 4.0)
waveform = component.signal
```

`ridge_mask` and `reconstruct_component` support fixed-Hz, fixed-bin,
log-relative, or custom masks. Bounds are inclusive. The returned mask and
weights remain inspectable.

For several components, overlapping masks require explicit ownership:

```julia
parts = reconstruct_components(
    result,
    ridges;
    bandwidth_hz = [4.0, 6.0],
    ownership = :exclusive,
)
```

`:exclusive` assigns a selected bin to the nearest active ridge, breaking ties
by ridge order. `:fractional` partitions it by inverse frequency distance.
The nonnegative component weights plus `residual_weights` sum to one at every
SST bin. Consequently, `sum(parts.components) + parts.residual` reproduces the
complete represented oscillatory WSST synthesis without double counting.
Neither policy claims a unique physical separation where modes overlap.

`ridge_metrics` reports frequency RMSE, bias, maximum error, coverage, and
boundary fraction over active samples. `component_metrics` reports relative
waveform error and SRER in dB without amplitude or phase fitting.

## Limitations

Frequency constraints help isolate weak known modes but cannot prove a mode is
present. Activity thresholds are deterministic evidence rules, not calibrated
probabilities. A coarse grid may also make a physically reasonable `max_rate`
infeasible because moving one bin in one sample already exceeds the rate.
Crossings, overlapping ranges, and reconstruction bandwidth remain modeling
choices that require sensitivity studies.

For locally smooth chirps, see [Predictive ridge tracking](@ref) for the exact
second-order and prediction-window-pruned strategies. The first-order solver
on this page remains the reference when a curvature model is not justified.
