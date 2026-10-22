# Julia and Python performance

This page separates two questions that are easy to conflate:

1. Where does the Julia pipeline spend its time?
2. Which operations can be compared strictly with `ssqueezepy==0.6.6`?

All quoted measurements are observations on one machine, not portable
performance promises. The machine was an Intel Core Ultra 7 155U running
Julia 1.13.0-rc1 with one Julia thread and one FFTW thread. Python was also
limited to one OpenMP, BLAS, and Numba thread, with GPU execution disabled.

## Julia stage profile (Milestone 05A)

The frozen 1024-sample profile in `validation/performance_05a.json` used 16
voices per octave and a 384-bin WSST grid. Planning was kept separate from
steady-state execution.

| Julia stage | Median | Allocated per call |
|---|---:|---:|
| CWT planning | 13.44 ms | 7.72 MiB |
| Backend CWT | 2.88 ms | 11.28 MiB |
| Spectral derivative | 2.02 ms | 4.53 MiB |
| Phase plus complex reassignment | 10.55 ms | 6.13 MiB |
| One constrained ridge | 10.52 ms | 9.25 MiB |
| Two constrained ridges | 21.18 ms | 18.49 MiB |
| Component reconstruction/ownership | 6.11 ms | 9.24 MiB |
| Envelope and activity analysis | 0.040 ms | 218 KiB |
| Damping estimation | 0.025 ms | 138 KiB |
| Full workflow after planning | 42.40 ms | 50.06 MiB |

The first ridge optimization reduced the representative single-ridge median
from 521.24 ms to 19.56 ms (26.65×) while preserving ridge indices and total
cost exactly. The table above is a later, narrower-frequency workload, which
is why its single-ridge time is lower than that optimization audit.

## Equivalence audit against Python (Milestone 05B)

The comparison uses the pinned environment under `reference/python/.venv`:
Python 3.13.9, NumPy 2.3.3, and `ssqueezepy==0.6.6`. The benchmark adapter was
derived from the installed source rather than from similarly named APIs.

| Operation | Julia implementation | Python implementation | Strict status |
|---|---|---|---|
| CWT | reused `CWTPlan` daughters and backend `cwt` | `ssqueezepy.cwt`, cached wavelet, no padding | comparable |
| CWT + physical derivative | `cwt(plan, x)` | `ssqueezepy.cwt(..., derivative=true)` | comparable |
| Phase | signed `imag(dW/W)/(2π)` | absolute-valued phase kernel | not strictly comparable |
| Reassignment | trapezoidal finite-bank log weights | constant `log(2)/nv` endpoint weights | not strictly comparable |
| Full first-order WSST | signed phase plus trapezoidal endpoints | absolute phase plus constant endpoints | not strictly comparable |
| Real standard inverse | `inverse_wsst(..., method=:standard)` | `issq_cwt` | comparable |
| Fixed supplied mask | `reconstruct_region` | equivalent NumPy masked sum with Python admissibility | comparable kernel; Python has no arbitrary-mask API |
| Analytic-complex inverse | supported | real synthesis only | unsupported in Python |
| Ridge tracking | constrained physical DP | different objective | excluded |

The generalized Morse filters themselves match. `ContinuousWavelets.jl`
samples a positive-frequency coordinate whose last point is ``N/2+1``; the
adapter therefore uses

```math
a_{py} = \frac{\omega_{peak} N}{2\pi s_{jl}}
         \frac{N/2+1}{N/2}.
```

After scale-order conversion, the Float64 CWT and derivative pass the
coefficient-level gate. `ssqueezepy` halves the analytic Nyquist filter bin;
the deterministic audit signal has zero Nyquist content so this convention
does not contaminate the execution comparison.

The phase and WSST differences are retained as failed equivalence results in
the JSON. Their measured times are useful for profiling each implementation,
but no Julia/Python speed ratio is published for them. In particular, the
requested primary full-WSST ratio is intentionally absent until both packages
execute the same signed phase estimator and finite-scale quadrature.

### Representative measured result

The committed report uses 2048 Float64 samples, 16 voices per octave (160
scales/bins), 12 warmed repetitions, and one CPU thread. The ratio is Python
median divided by Julia median, so values below one mean Python was faster.

| Operation | Julia median | Python median | Python / Julia | Gate |
|---|---:|---:|---:|---|
| CWT | 6.15 ms | 3.30 ms | 0.54× | pass; Python faster |
| CWT plus derivative | 7.03 ms | 6.99 ms | 0.99× | pass; comparable on this run |
| Phase estimation | 16.15 ms | 2.10 ms | — | fail; different sign rule |
| Reassignment | 23.86 ms | 7.82 ms | — | fail; 1.07% relative L2 |
| Complete first-order WSST | 32.94 ms | 18.94 ms | — | fail; 1.07% relative L2 |
| Standard real inverse | 0.188 ms | 1.637 ms | 8.73× | pass; Julia faster |
| Fixed-band supplied mask | 0.464 ms | 0.527 ms | 1.14× | pass; Julia faster |

The CWT coefficient relative error was ``1.45\times10^{-15}``; the derivative
error was ``1.54\times10^{-15}``. Standard and fixed-band inverse errors were
both below ``6.1\times10^{-16}``. The inverse timing uses each package's public
standard synthesis path; `issq_cwt` recomputes its admissibility value per
call, while the Julia result carries that metadata. This is an API-level cost,
not evidence that raw array summation is eight times faster in Julia.

The near-equal complete-WSST medians are shown to make the failed comparison
visible, not to imply parity under matched mathematics. Profiling identifies
Python's compiled phase/reassignment kernels and Julia's scalar fused
reassignment loop as the main implementation difference, while the failed
gate prevents turning that observation into a language speedup claim.

## Reproduce

Instantiate both pinned environments, then run the representative audit:

```sh
julia --project=benchmarks -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
reference/python/.venv/bin/pip install -r reference/python/requirements.txt
JULIA_NUM_THREADS=1 julia --project=benchmarks \
  benchmarks/julia_python_comparison.jl \
  2048 16 12 validation/julia_python_benchmarks.json Float64
```

Arguments are signal length, voices per octave, repetitions, JSON path, and
precision. This makes the requested workload matrix reproducible without
mixing results from different processes. Suggested runs are:

```sh
# 160, 320, and 640-scale representative cases
for q in 16 32 64; do
  JULIA_NUM_THREADS=1 julia --project=benchmarks \
    benchmarks/julia_python_comparison.jl 2048 "$q" 12 \
    "benchmarks/results/python_2048_q${q}.json" Float64
done

# Length scaling; use a fresh process for every precision and size
for n in 2048 8192 32768 131072; do
  for precision in Float32 Float64; do
    JULIA_NUM_THREADS=1 julia --project=benchmarks \
      benchmarks/julia_python_comparison.jl "$n" 16 12 \
      "benchmarks/results/python_${n}_${precision}.json" "$precision"
  done
done
```

The 131072-sample case needs substantial memory because both complex CWT and
derivative matrices are live. Run it separately. Each JSON report includes
minimum, median, interquartile range, cold-call timing, allocations or traced
Python memory, versions, CPU/thread configuration, numerical errors, and a
ratio only for gates that passed. Julia allocated bytes and Python
`tracemalloc` peaks are different metrics and should not be compared as if
they measured the same allocator boundary.
