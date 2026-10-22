#!/usr/bin/env julia

using BenchmarkTools
using Dates
using FFTW
using JSON
using Pkg
using Statistics
using SynchroSqueezing

const DEFAULT_SAMPLES = 2048
const DEFAULT_OUTPUT = joinpath(@__DIR__, "results", "full_pipeline.json")

function benchmark_summary(operation; samples = 20, seconds = 30)
    GC.gc()
    benchmark = @benchmarkable $operation() evals = 1
    trial = run(benchmark; samples, seconds)
    times = sort!(Float64.(trial.times) ./ 1e6)
    estimate = median(trial)
    return Dict(
        "minimum_ms" => minimum(times),
        "median_ms" => Statistics.median(times),
        "maximum_ms" => maximum(times),
        "q25_ms" => Statistics.quantile(times, 0.25),
        "q75_ms" => Statistics.quantile(times, 0.75),
        "samples" => length(times),
        "allocations" => estimate.allocs,
        "allocated_bytes" => estimate.memory,
    )
end

function synthetic_case(n, fs)
    times = collect(0:n-1) ./ fs
    first = exp.(-times ./ 0.65) .* cos.(2π .* (34 .* times .+ 4 .* times .^ 2))
    second = 0.65 .* exp.(-times ./ 2.4) .*
             cos.(2π .* (92 .* times .- 2 .* times .^ 2))
    return times, first .+ second
end

function pipeline_after_planning(plan, signal, times)
    transformed = cwt(plan, signal)
    squeezed = wsst(
        transformed;
        freqrange = (20.0, 115.0),
        nfrequencies = 384,
        threshold = RelativeThreshold(1e-8),
        preserve_residual = true,
    )
    ridges = extract_ridges(
        squeezed;
        freqranges = [(27.0, 58.0), (75.0, 105.0)],
        penalty = 1.0,
        continuity_scale = 1.0,
        max_rate = 40.0,
        activity_threshold = 1e-8,
    )
    components = reconstruct_components(
        squeezed,
        ridges;
        bandwidth_hz = 6.0,
        ownership = :exclusive,
        analytic = true,
    )
    features = [component_features(
        components.components[index],
        times;
        valid_mask = ridges[index].valid .& ridges[index].active,
    ) for index in eachindex(ridges)]
    supports = [mode_support(
        feature;
        relative_amplitude_threshold = 0.05,
        persistence_seconds = 0.02,
        boundary_mask = .!ridges[index].at_frequency_boundary,
    ) for (index, feature) in enumerate(features)]
    fit_window = (0.2, min(1.5, times[end] - 0.1))
    fits = [fit_exponential_damping(feature; window = fit_window)
            for feature in features]
    return (; transformed, squeezed, ridges, components, features, supports, fits)
end

function package_versions()
    dependencies = Pkg.dependencies()
    selected = Dict{String,String}()
    for dependency in values(dependencies)
        dependency.name in ("BenchmarkTools", "ContinuousWavelets", "FFTW", "JSON") ||
            continue
        selected[dependency.name] = string(dependency.version)
    end
    return selected
end

function main(output_path = DEFAULT_OUTPUT)
    n = isempty(ARGS) ? DEFAULT_SAMPLES : parse(Int, ARGS[1])
    samples = length(ARGS) < 2 ? 20 : parse(Int, ARGS[2])
    n >= 256 || throw(ArgumentError("signal length must be at least 256"))
    samples > 0 || throw(ArgumentError("benchmark sample count must be positive"))
    fs = 512.0
    times, signal = synthetic_case(n, fs)

    cold_plan_seconds = @elapsed plan = plan_cwt(n; fs, voices_per_octave = 16)
    cold_pipeline_seconds = @elapsed pipeline = pipeline_after_planning(plan, signal, times)
    transformed = pipeline.transformed
    squeezed = pipeline.squeezed
    ridges = pipeline.ridges
    components = pipeline.components
    features = pipeline.features
    centered = signal .- sum(signal) / length(signal)

    stages = Dict{String,Any}()
    stages["cwt_initialization_and_planning"] = benchmark_summary(
        () -> plan_cwt(n; fs, voices_per_octave = 16); samples,
    )
    stages["backend_cwt_execution"] = benchmark_summary(
        () -> SynchroSqueezing._CW.cwt(centered, plan.backend_config, plan.daughters); samples,
    )
    stages["spectral_differentiation"] = benchmark_summary(
        () -> SynchroSqueezing._spectral_derivative(centered, plan); samples,
    )
    stages["phase_estimation"] = benchmark_summary(
        () -> instantaneous_frequency(transformed; threshold = RelativeThreshold(1e-8)); samples,
    )
    stages["phase_and_complex_reassignment"] = benchmark_summary(
        () -> wsst(
            transformed;
            freqrange = (20.0, 115.0),
            nfrequencies = 384,
            threshold = RelativeThreshold(1e-8),
            preserve_residual = true,
        ); samples,
    )
    ridge_options = (
        penalty = 1.0,
        continuity_scale = 1.0,
        max_rate = 40.0,
        activity_threshold = 1e-8,
    )
    stages["single_ridge_dynamic_programming"] = benchmark_summary(
        () -> extract_ridge(squeezed; freqrange = (27.0, 58.0), ridge_options...); samples,
    )
    stages["multiple_frequency_constrained_ridges"] = benchmark_summary(
        () -> extract_ridges(
            squeezed;
            freqranges = [(27.0, 58.0), (75.0, 105.0)],
            ridge_options...,
        ); samples,
    )
    stages["component_reconstruction_and_ownership"] = benchmark_summary(
        () -> reconstruct_components(
            squeezed,
            ridges;
            bandwidth_hz = 6.0,
            ownership = :exclusive,
            analytic = true,
        ); samples,
    )
    stages["envelope_and_activity_analysis"] = benchmark_summary(
        () -> [
            begin
                feature = component_features(
                    components.components[index],
                    times;
                    valid_mask = ridges[index].valid .& ridges[index].active,
                )
                mode_support(
                    feature;
                    relative_amplitude_threshold = 0.05,
                    persistence_seconds = 0.02,
                    boundary_mask = .!ridges[index].at_frequency_boundary,
                )
            end for index in eachindex(ridges)
        ]; samples,
    )
    stages["damping_estimation"] = benchmark_summary(
        () -> [fit_exponential_damping(feature; window = (0.2, 0.8))
               for feature in features]; samples,
    )
    stages["end_to_end_after_planning"] = benchmark_summary(
        () -> pipeline_after_planning(plan, signal, times); samples = max(5, samples ÷ 2),
    )

    result = Dict(
        "schema_version" => 1,
        "milestone" => "05A",
        "generated_at_utc" => string(Dates.now(Dates.UTC)),
        "julia_version" => string(VERSION),
        "package_versions" => package_versions(),
        "configuration" => Dict(
            "samples" => n,
            "sample_rate_hz" => fs,
            "voices_per_octave" => 16,
            "wsst_frequency_bins" => 384,
            "julia_threads" => Threads.nthreads(),
            "fftw_threads" => FFTW.get_num_threads(),
            "benchmark_evals" => 1,
            "requested_samples_per_stage" => samples,
        ),
        "cold_current_process" => Dict(
            "planning_seconds" => cold_plan_seconds,
            "pipeline_after_planning_seconds" => cold_pipeline_seconds,
            "note" => "Includes first-call compilation in this Julia process; package startup is excluded.",
        ),
        "process_maxrss_bytes" => Sys.maxrss(),
        "stages" => stages,
        "methodology" => [
            "Each steady-state stage is warmed before BenchmarkTools records it.",
            "Every evaluation invokes a function over runtime-owned arrays with evals=1.",
            "Planning and first-call compilation are reported separately.",
            "The two ridge modes reuse one CWT and one WSST result.",
            "Allocated bytes are per evaluation; process max RSS is a process-wide high-water mark.",
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
    output = length(ARGS) < 3 ? DEFAULT_OUTPUT : abspath(ARGS[3])
    main(output)
end
