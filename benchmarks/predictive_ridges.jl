#!/usr/bin/env julia

using BenchmarkTools
using SynchroSqueezing

function main()
    nfrequency = isempty(ARGS) ? 64 : parse(Int, first(ARGS))
    ntime = length(ARGS) < 2 ? 256 : parse(Int, ARGS[2])
    frequencies = collect(range(20.0, 100.0; length = nfrequency))
    times = collect(range(0.0, 2.0; length = ntime))
    truth = 35.0 .+ 15.0 .* times
    coefficients = ComplexF64[
        exp(-0.5 * ((frequency - truth[column]) / 1.0)^2)
        for frequency in frequencies, column in eachindex(times)
    ]
    common = (
        times = times,
        freqrange = (28.0, 72.0),
        max_rate = 2 * (frequencies[2] - frequencies[1]) / (times[2] - times[1]),
        frequency_penalty = 0.1,
        frequency_scale = 1.0,
        curvature_penalty = 1.0,
        curvature_scale = 20.0,
    )

    println("first-order banded DP")
    display(@benchmark extract_ridge(
        $coefficients,
        $frequencies;
        times = $times,
        freqrange = (28.0, 72.0),
        max_rate = $(common.max_rate),
        penalty = 0.1,
        continuity_scale = 1.0,
    ))
    println("exact second-order DP")
    display(@benchmark extract_predictive_ridge(
        $coefficients,
        $frequencies;
        $common...,
        method = :exact,
    ))
    println("prediction-window-pruned second-order DP")
    display(@benchmark extract_predictive_ridge(
        $coefficients,
        $frequencies;
        $common...,
        method = :pruned,
        prediction_bandwidth = 2.0,
    ))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
