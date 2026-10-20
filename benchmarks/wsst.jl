#!/usr/bin/env julia

using BenchmarkTools
using SynchroSqueezing

function main()
    fs = 2048.0
    n = isempty(ARGS) ? 4096 : parse(Int, first(ARGS))
    time = (0:n-1) ./ fs
    signal = cos.(2π .* (80.0 .* time .+ 60.0 .* time .^ 2))
    plan = plan_cwt(n; fs, voices_per_octave = 16)
    transformed = cwt(plan, signal)
    squeezed = wsst(
        transformed;
        freqrange = (40.0, 400.0),
        nfrequencies = 512,
        preserve_residual = true,
    )
    ridge_frequency = 80.0 .+ 120.0 .* time
    component_mask = abs.(reshape(squeezed.frequencies, :, 1) .-
                          reshape(ridge_frequency, 1, :)) .<= 8.0

    println("phase estimation")
    display(@benchmark instantaneous_frequency($transformed; threshold=RelativeThreshold(1e-6)))
    println("nearest reassignment after CWT")
    display(@benchmark wsst($transformed; freqrange=(40.0, 400.0), nfrequencies=512,
                            reassignment=:nearest, preserve_residual=false))
    println("linear reassignment after CWT")
    display(@benchmark wsst($transformed; freqrange=(40.0, 400.0), nfrequencies=512,
                            reassignment=:linear, preserve_residual=false))
    println("end-to-end CWT and WSST")
    display(@benchmark wsst($plan, $signal; freqrange=(40.0, 400.0),
                            nfrequencies=512, preserve_residual=false))
    println("standard inverse (frequency reduction and scalar synthesis)")
    display(@benchmark inverse_wsst($squeezed))
    println("masked standard component inverse")
    display(@benchmark reconstruct_region($squeezed, $component_mask))
    println("finite-bank frame inverse")
    display(@benchmark inverse_wsst($squeezed; mode=:full, method=:frame))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
