#!/usr/bin/env julia

using BenchmarkTools
using Dates
using FFTW
using JSON
using LinearAlgebra
using Pkg
using Statistics
using SynchroSqueezing

const DEFAULT_OUTPUT = joinpath(@__DIR__, "..", "validation", "julia_python_benchmarks.json")
const PYTHON = joinpath(@__DIR__, "..", "reference", "python", ".venv", "bin", "python")
const PYTHON_SCRIPT = joinpath(@__DIR__, "python_comparison.py")

function benchmark_summary(operation; samples::Int = 12, seconds::Real = 30)
    operation()
    GC.gc()
    benchmark = @benchmarkable $operation() evals = 1
    trial = run(benchmark; samples, seconds)
    times = sort!(Float64.(trial.times) ./ 1e6)
    estimate = median(trial)
    return Dict(
        "minimum_ms" => minimum(times),
        "median_ms" => Statistics.median(times),
        "q25_ms" => Statistics.quantile(times, 0.25),
        "q75_ms" => Statistics.quantile(times, 0.75),
        "samples" => length(times),
        "allocations" => estimate.allocs,
        "allocated_bytes" => estimate.memory,
    )
end

function synthetic_signal(::Type{T}, n::Int, fs::T) where {T<:AbstractFloat}
    time = T.(0:n-1) ./ fs
    duration = T(n) / fs
    # Integer-period carriers and envelopes keep DC and Nyquist leakage at
    # roundoff, isolating the filter-bank comparison from endpoint conventions.
    envelope1 = T(0.75) .+ T(0.20) .* cos.(T(2π) .* time ./ duration)
    envelope2 = T(0.65) .+ T(0.15) .* sin.(T(4π) .* time ./ duration)
    signal = envelope1 .* cos.(T(2π * 31) .* time) .+
             envelope2 .* cos.(T(2π * 67) .* time .+ T(0.35)) .+
             T(0.45) .* cos.(T(2π * 103) .* time .- T(0.2))
    signal .-= sum(signal) / T(n)
    return signal
end

function dtype_string(::Type{Float32})
    return "<f4"
end
function dtype_string(::Type{Float64})
    return "<f8"
end
function dtype_string(::Type{ComplexF32})
    return "<c8"
end
function dtype_string(::Type{ComplexF64})
    return "<c16"
end
dtype_string(::Type{Bool}) = "|u1"

function write_array!(arrays, root, name, array)
    materialized = Array(array)
    file = name * ".bin"
    open(joinpath(root, file), "w") do stream
        write(stream, materialized)
    end
    arrays[name] = Dict(
        "file" => file,
        "dtype" => dtype_string(eltype(materialized)),
        "shape" => collect(size(materialized)),
    )
    return materialized
end

function fixture_manifest(root, plan, transformed, squeezed, signal, threshold, precision)
    arrays = Dict{String,Any}()
    phase = instantaneous_frequency(transformed; threshold)
    inverse = inverse_wsst(squeezed; mode = :oscillatory, method = :standard)
    mask = repeat((squeezed.frequencies .>= 28) .& (squeezed.frequencies .<= 82), 1,
                  length(signal))
    band_inverse = reconstruct_region(squeezed, mask; method = :standard)
    for (name, value) in (
        "signal" => signal,
        "scales" => transformed.scales,
        "frequencies" => transformed.frequencies,
        "cwt" => transformed.coefficients,
        "derivative" => transformed.derivative,
        "phase" => phase.values,
        "wsst" => squeezed.coefficients,
        "inverse" => inverse,
        "mask" => mask,
        "band_inverse" => band_inverse,
    )
        write_array!(arrays, root, name, value)
    end
    relative = precision === Float32 ? 2e-4 : 5e-11
    absolute = precision === Float32 ? 2e-4 : 5e-11
    tolerances = Dict(name => Dict("relative_l2" => relative, "maximum_absolute" => absolute)
                      for name in ("cwt", "derivative", "phase", "standard_inverse",
                                   "fixed_band_inverse"))
    # Reassignment is deliberately held to the same gate. If the packages'
    # finite-scale quadratures differ, it is reported as unsupported.
    tolerances["reassignment"] = Dict("relative_l2" => relative, "maximum_absolute" => absolute)
    tolerances["complete_wsst"] = copy(tolerances["reassignment"])
    precision === Float32 && (tolerances["derivative"]["maximum_absolute"] = 2e-3)
    manifest = Dict(
        "arrays" => arrays,
        "configuration" => Dict(
            "samples" => length(signal),
            "sample_rate_hz" => plan.sample_rate,
            "voices_per_octave" => plan.backend_config.Q,
            "frequency_bins" => length(transformed.frequencies),
            "precision" => string(precision),
            "wavelet" => "generalized Morse gamma=3 beta=20, analytic L1/bandpass",
            "boundary" => "periodic/no padding",
            "absolute_threshold" => threshold.value,
            "admissibility" => real(transformed.admissibility),
            "reassignment" => "nearest center",
        ),
        "tolerances" => tolerances,
    )
    open(joinpath(root, "manifest.json"), "w") do stream
        JSON.print(stream, manifest, 2)
        println(stream)
    end
    return (; phase, inverse, mask, band_inverse)
end

function package_versions()
    selected = Dict{String,String}()
    for dependency in values(Pkg.dependencies())
        dependency.name in ("BenchmarkTools", "ContinuousWavelets", "FFTW", "JSON") || continue
        selected[dependency.name] = string(dependency.version)
    end
    return selected
end

function run_python(root, output, samples)
    isfile(PYTHON) || error("pinned Python environment is missing: $PYTHON")
    command = `$PYTHON $PYTHON_SCRIPT --fixture $root --output $output --samples $samples`
    environment = copy(ENV)
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                "NUMBA_NUM_THREADS")
        environment[key] = "1"
    end
    environment["SSQ_GPU"] = "0"
    run(setenv(command, environment))
    return JSON.parsefile(output)
end

function speedup_table(julia_timings, python_result)
    comparable = python_result["strict_comparability"]
    ratios = Dict{String,Any}()
    for (stage, julia_result) in julia_timings
        python_result_stage = python_result["timings"][stage]
        ratios[stage] = Dict(
            "strictly_comparable" => comparable[stage],
            "julia_median_ms" => julia_result["median_ms"],
            "python_median_ms" => python_result_stage["median_ms"],
            "python_over_julia" => comparable[stage] ?
                python_result_stage["median_ms"] / julia_result["median_ms"] : nothing,
        )
    end
    return ratios
end

function main(output_path = DEFAULT_OUTPUT; n::Int = 2048, q::Int = 16,
              precision::Type{T} = Float64, samples::Int = 12) where {T<:AbstractFloat}
    FFTW.set_num_threads(1)
    fs = T(512)
    signal = synthetic_signal(T, n, fs)
    threshold = AbsoluteThreshold(T(1e-8))

    first_plan_ms = 1e3 * @elapsed plan = plan_cwt(n; fs, voices_per_octave = q,
                                                   precision)
    first_complete_ms = 1e3 * @elapsed begin
        transformed = cwt(plan, signal)
        squeezed = wsst(transformed; frequency_grid = transformed.frequencies,
                        threshold, preserve_residual = false)
    end
    extras = nothing
    python_result = nothing
    mktempdir() do root
        extras = fixture_manifest(root, plan, transformed, squeezed, signal, threshold, precision)
        python_output = joinpath(root, "python.json")
        python_result = run_python(root, python_output, samples)
    end

    centered = signal .- sum(signal) / length(signal)
    frequencies = transformed.frequencies
    julia_timings = Dict{String,Any}(
        "cwt" => benchmark_summary(
            () -> SynchroSqueezing._CW.cwt(centered, plan.backend_config, plan.daughters); samples),
        "cwt_and_derivative" => benchmark_summary(() -> cwt(plan, signal); samples),
        "phase_estimation" => benchmark_summary(
            () -> instantaneous_frequency(transformed; threshold); samples),
        "reassignment" => benchmark_summary(
            () -> wsst(transformed; frequency_grid = frequencies, threshold,
                       preserve_residual = false); samples),
        "complete_first_order_wsst" => benchmark_summary(
            () -> begin
                current = cwt(plan, signal)
                wsst(current; frequency_grid = frequencies, threshold,
                     preserve_residual = false)
            end; samples = max(4, samples ÷ 2)),
        "standard_inverse_wsst" => benchmark_summary(
            () -> inverse_wsst(squeezed; mode = :oscillatory, method = :standard); samples),
        "fixed_band_reconstruction" => benchmark_summary(
            () -> reconstruct_region(squeezed, extras.mask; method = :standard); samples),
    )

    result = Dict(
        "schema_version" => 1,
        "milestone" => "05B",
        "generated_at_utc" => string(Dates.now(Dates.UTC)),
        "primary_metric" => "complete_first_order_wsst when equivalence gate passes",
        "configuration" => python_result["configuration"],
        "environment" => Dict(
            "julia_version" => string(VERSION),
            "julia_threads" => Threads.nthreads(),
            "fftw_threads" => FFTW.get_num_threads(),
            "cpu" => Sys.cpu_info()[1].model,
            "logical_cores" => Sys.CPU_THREADS,
            "package_versions" => package_versions(),
        ),
        "julia" => Dict(
            "cold_current_process" => Dict(
                "plan_ms" => first_plan_ms,
                "complete_after_plan_ms" => first_complete_ms,
            ),
            "timings" => julia_timings,
            "process_maxrss_bytes" => Sys.maxrss(),
            "resident_array_bytes" => Dict(
                "reusable_plan" => Base.summarysize(plan),
                "signal" => Base.summarysize(signal),
                "cwt_coefficients" => Base.summarysize(transformed.coefficients),
                "cwt_derivative" => Base.summarysize(transformed.derivative),
                "wsst_coefficients" => Base.summarysize(squeezed.coefficients),
                "fixed_band_mask" => Base.summarysize(extras.mask),
            ),
        ),
        "python" => python_result,
        "comparison" => speedup_table(julia_timings, python_result),
        "methodology" => [
            "Both implementations run with one application, FFT/BLAS, and Numba thread.",
            "Julia BenchmarkTools and Python perf_counter_ns measure warmed native calls; process startup and fixture I/O are excluded.",
            "The same typed signal, sampled Fourier response, scales, physical derivative units, threshold, and nearest-bin centers are used.",
            "A ratio is emitted only when the numerical equivalence gate for that operation passes.",
            "Julia allocation bytes and Python tracemalloc peaks use different accounting and must not be compared as identical metrics.",
        ],
    )
    mkpath(dirname(output_path))
    open(output_path, "w") do stream
        JSON.print(stream, result, 2)
        println(stream)
    end
    println("saved ", output_path)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    n = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2048
    q = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 16
    samples = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 12
    output = length(ARGS) >= 4 ? abspath(ARGS[4]) : DEFAULT_OUTPUT
    precision = length(ARGS) >= 5 && ARGS[5] == "Float32" ? Float32 : Float64
    main(output; n, q, precision, samples)
end
