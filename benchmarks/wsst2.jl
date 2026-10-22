#!/usr/bin/env julia

using BenchmarkTools
using Dates
using FFTW
using JSON
using LinearAlgebra: norm
using Pkg
using Statistics
using SynchroSqueezing

const DEFAULT_OUTPUT = joinpath(@__DIR__, "..", "validation", "wsst2_06.json")

function benchmark_summary(operation; samples = 10, seconds = 30)
    operation()
    GC.gc()
    benchmark = @benchmarkable $operation() evals = 1
    trial = run(benchmark; samples, seconds)
    times = sort!(Float64.(trial.times) ./ 1e6)
    estimate = median(trial)
    return Dict(
        "minimum_ms" => minimum(times),
        "median_ms" => Statistics.median(times),
        "maximum_ms" => maximum(times),
        "samples" => length(times),
        "allocations" => estimate.allocs,
        "allocated_bytes" => estimate.memory,
    )
end

function dominant_frequency(result)
    rows = [argmax(abs2.(@view result.coefficients[:, column]))
            for column in axes(result.coefficients, 2)]
    return result.frequencies[rows]
end

function frequency_metrics(estimate, truth, indices)
    error = estimate[indices] .- truth[indices]
    return Dict(
        "rmse_hz" => sqrt(mean(abs2, error)),
        "bias_hz" => mean(error),
        "maximum_error_hz" => maximum(abs, error),
    )
end

function median_frequency_width(result, truth, indices)
    widths = Float64[]
    for column in indices
        energy = abs2.(@view result.coefficients[:, column])
        total = sum(energy)
        total > 0 || continue
        push!(widths, sqrt(sum(energy .* (result.frequencies .- truth[column]) .^ 2) / total))
    end
    return median(widths)
end

function chirp_case(plan, times, initial_frequency, chirp_rate, common)
    phase = 2π .* (initial_frequency .* times .+ 0.5chirp_rate .* times .^ 2)
    signal = cos.(phase)
    operators = second_order_operators(plan, signal)
    estimate = second_order_frequency(
        operators;
        threshold = RelativeThreshold(1e-5),
    )
    first = wsst(operators.cwt; order = 1, common...)
    second = wsst(operators; order = 2, common...)
    truth = initial_frequency .+ chirp_rate .* times
    interior = 129:(length(times) - 128)
    boundary = vcat(1:128, (length(times) - 127):length(times))
    first_frequency = dominant_frequency(first)
    second_frequency = dominant_frequency(second)

    local_rates = Float64[]
    for column in interior
        row = argmax(abs.(@view operators.cwt.coefficients[:, column]))
        estimate.second_order_valid[row, column] || continue
        push!(local_rates, estimate.chirp_rate[row, column])
    end
    first_reconstruction = inverse_wsst(first; mode = :full, method = :frame)
    second_reconstruction = inverse_wsst(second; mode = :full, method = :frame)
    return Dict(
        "chirp_rate_hz_per_s" => chirp_rate,
        "first_order" => Dict(
            "interior_frequency" => frequency_metrics(first_frequency, truth, interior),
            "boundary_frequency" => frequency_metrics(first_frequency, truth, boundary),
            "median_frequency_width_hz" => median_frequency_width(first, truth, interior),
            "full_reconstruction_relative_l2" => norm(first_reconstruction - signal) / norm(signal),
        ),
        "second_order" => Dict(
            "interior_frequency" => frequency_metrics(second_frequency, truth, interior),
            "boundary_frequency" => frequency_metrics(second_frequency, truth, boundary),
            "median_frequency_width_hz" => median_frequency_width(second, truth, interior),
            "full_reconstruction_relative_l2" => norm(second_reconstruction - signal) / norm(signal),
            "median_chirp_rate_hz_per_s" => isempty(local_rates) ? NaN : median(local_rates),
            "chirp_rate_absolute_error_hz_per_s" =>
                isempty(local_rates) ? NaN : abs(median(local_rates) - chirp_rate),
            "second_order_valid_fraction" =>
                second.reassignment_metadata.second_order_valid_fraction,
            "fallback_fraction" => second.reassignment_metadata.fallback_fraction,
            "excluded_fraction" => second.reassignment_metadata.excluded_fraction,
        ),
    )
end

function package_versions()
    selected = Dict{String,String}()
    for dependency in values(Pkg.dependencies())
        dependency.name in ("BenchmarkTools", "ContinuousWavelets", "FFTW", "JSON") ||
            continue
        selected[dependency.name] = string(dependency.version)
    end
    return selected
end

function main(output_path = DEFAULT_OUTPUT)
    n = 1024
    fs = 512.0
    times = collect(0:n-1) ./ fs
    plan = plan_cwt(n; fs, voices_per_octave = 32)
    common = (
        freqrange = (20.0, 220.0),
        nfrequencies = 900,
        threshold = RelativeThreshold(1e-5),
        preserve_residual = true,
    )
    comparisons = [chirp_case(plan, times, 35.0, rate, common)
                   for rate in (0.0, 20.0, 80.0)]

    benchmark_rate = 80.0
    signal = cos.(2π .* (35.0 .* times .+ 0.5benchmark_rate .* times .^ 2))
    transformed = cwt(plan, signal)
    operators = second_order_operators(plan, signal)
    stages = Dict{String,Any}(
        "cwt_and_first_derivative" => benchmark_summary(() -> cwt(plan, signal)),
        "additional_second_order_operators_end_to_end" =>
            benchmark_summary(() -> second_order_operators(plan, signal)),
        "second_order_frequency_from_operators" => benchmark_summary(
            () -> second_order_frequency(operators; threshold = RelativeThreshold(1e-5)),
        ),
        "wsst1_reassignment_after_cwt" => benchmark_summary(
            () -> wsst(transformed; order = 1, common...),
        ),
        "wsst2_reassignment_after_operators" => benchmark_summary(
            () -> wsst(operators; order = 2, common...),
        ),
        "wsst1_end_to_end" => benchmark_summary(
            () -> wsst(plan, signal; order = 1, common...); samples = 6,
        ),
        "wsst2_end_to_end" => benchmark_summary(
            () -> wsst(plan, signal; order = 2, common...); samples = 6,
        ),
    )

    result = Dict(
        "schema_version" => 1,
        "milestone" => "06",
        "generated_at_utc" => string(Dates.now(Dates.UTC)),
        "reference" => Dict(
            "implementation" => "Pham-Meignen WSST2",
            "repository" => "https://github.com/phamduonghung/WSST2",
            "commit" => "2721d07e943abcc6faf7bfcaf974b5dfad597882",
        ),
        "environment" => Dict(
            "julia_version" => string(VERSION),
            "package_versions" => package_versions(),
            "julia_threads" => Threads.nthreads(),
            "fftw_threads" => FFTW.get_num_threads(),
            "process_maxrss_bytes" => Sys.maxrss(),
        ),
        "configuration" => Dict(
            "samples" => n,
            "sample_rate_hz" => fs,
            "voices_per_octave" => 32,
            "frequency_range_hz" => [20.0, 220.0],
            "frequency_bins" => 900,
            "threshold_relative" => 1e-5,
            "interior_samples" => [129, n - 128],
            "boundary_samples_per_side" => 128,
            "boundary" => "periodic",
        ),
        "chirp_rate_regime" => comparisons,
        "performance" => stages,
        "notes" => [
            "Both orders use the same plan, wavelet bank, threshold, output grid, and deposition rule.",
            "WSST2 currently has a scalar reference implementation and no package-level threaded kernel.",
            "Peak RSS is process-wide; allocated bytes are per warmed BenchmarkTools evaluation.",
            "Raw MATLAB matrices are not compared because normalization, scale order, units, and boundaries differ; the independently published equations and linear-chirp invariant are used instead.",
        ],
    )
    mkpath(dirname(output_path))
    open(output_path, "w") do stream
        JSON.print(stream, result, 2)
        println(stream)
    end
    println(JSON.json(result, 2))
    println("saved ", output_path)
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    output = isempty(ARGS) ? DEFAULT_OUTPUT : abspath(first(ARGS))
    main(output)
end
