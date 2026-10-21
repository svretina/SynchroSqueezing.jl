#!/usr/bin/env julia

using BenchmarkTools
using SynchroSqueezing

function main()
    nfrequency = isempty(ARGS) ? 320 : parse(Int, first(ARGS))
    ntime = 2048
    frequencies = collect(range(20.0, 180.0; length = nfrequency))
    times = collect(0:ntime-1) ./ 1024.0
    truth = 60.0 .+ 25.0 .* times
    coefficients = ComplexF64[
        exp(-0.5 * ((frequency - truth[column]) / 2.0)^2)
        for frequency in frequencies, column in eachindex(times)
    ]
    common = (
        times = times,
        freqrange = (40.0, 130.0),
        penalty = 1.0,
        continuity_scale = 1.0,
        max_rate = 600.0,
    )
    ridge = extract_ridge(coefficients, frequencies; method = :dp, common...)
    mask = abs.(reshape(frequencies, :, 1) .- reshape(ridge.frequencies, 1, :)) .<= 4

    println("maximum-energy baseline")
    display(@benchmark extract_ridge($coefficients, $frequencies; method=:maximum, $common...))
    println("dense dynamic programming")
    display(@benchmark extract_ridge($coefficients, $frequencies; method=:dense, $common...))
    println("banded dynamic programming")
    display(@benchmark extract_ridge($coefficients, $frequencies; method=:dp, $common...))
    println("mask materialization")
    display(@benchmark copy($mask))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
