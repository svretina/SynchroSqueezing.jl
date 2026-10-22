# WSST benchmarks

Use the isolated benchmark environment so BenchmarkTools is not a production
dependency:

```sh
julia --project=benchmarks -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=benchmarks benchmarks/wsst.jl
```

The script reports warmed runtimes and allocations for phase estimation,
nearest and linear WSST after CWT, the complete CWT+WSST pipeline, standard
whole-matrix synthesis, arbitrary-mask standard synthesis, and the separate
FFT-based finite-bank frame inverse.

Frequency-constrained ridge solvers can be benchmarked independently with:

```sh
julia --project=benchmarks benchmarks/ridges.jl
```

The optional first argument sets the number of frequency bins. The script
reports maximum-energy, dense-DP, and physically banded-DP costs separately.

Compare first-order, exact second-order, and prediction-window-pruned costs
with:

```sh
julia --project=benchmarks benchmarks/predictive_ridges.jl
```

Optional arguments set the number of frequency bins and time samples. Exact
second-order storage grows quadratically with the number of bins inside the
search range; use representative constrained ranges when benchmarking.

## Full scientific pipeline

Milestone 05A adds a stage-by-stage benchmark that covers planning, backend
CWT execution, spectral differentiation, phase estimation and reassignment,
single- and multiple-ridge extraction, reconstruction and ownership,
envelope/activity analysis, damping estimation, and the complete workflow:

```sh
julia --project=benchmarks -e 'using Pkg; Pkg.resolve(); Pkg.instantiate()'
julia --project=benchmarks benchmarks/full_pipeline.jl
```

Optional arguments are the signal length, samples per stage, and JSON output
path. For example:

```sh
julia --project=benchmarks benchmarks/full_pipeline.jl \
    4096 30 benchmarks/results/full_pipeline_4096.json
```

The JSON report records the full runtime distribution, allocations, allocated
bytes, Julia and FFTW thread counts, package versions, planning time, and
first-call cost. Steady-state evaluations use `evals=1`; the input arrays are
owned at runtime and reused without embedding constants in generated code.
The multiridge and downstream stages deliberately reuse one CWT and one WSST.

`process_maxrss_bytes` is the process-wide high-water mark, not a per-stage
peak. Run each thread configuration in a fresh process when comparing memory
or Julia thread counts. Do not compare the reported end-to-end time with a
Python implementation unless it performs the same transform, ridge,
reconstruction, activity, and damping operations.

## Julia versus pinned ssqueezepy

Run the equivalence-gated comparison in a single-threaded fresh process:

```sh
JULIA_NUM_THREADS=1 julia --project=benchmarks \
    benchmarks/julia_python_comparison.jl \
    2048 16 12 validation/julia_python_benchmarks.json Float64
```

The arguments select signal length, voices per octave, repetitions, output,
and precision. The driver creates typed temporary fixtures, invokes the pinned
`reference/python/.venv` outside every timed region, and combines numerical
gates and native timings in one JSON document. It emits a cross-language
ratio only for an operation whose gate passed. See the documentation page
"Julia and Python performance" for the audited function mapping and workload
matrix commands.

## Second-order WSST

Run the Milestone 06 first- versus second-order comparison with:

```bash
julia --project=benchmarks benchmarks/wsst2.jl
```

It writes `validation/wsst2_06.json`, including analytical chirp regimes,
interior and boundary accuracy, concentration, coefficient validity,
reconstruction, runtime, allocations, thread counts, and package versions.
