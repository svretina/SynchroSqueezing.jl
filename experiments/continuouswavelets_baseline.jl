#!/usr/bin/env julia

using ContinuousWavelets
using LinearAlgebra
using Printf

include("support.jl")
using .BackendExperimentSupport

function run_baseline(; n::Int = 2048, fs::Float64 = 1024.0)
    signals = analytical_signals(n, fs)
    config = wavelet(Morlet(π); averagingType = NoAve(), β = 2, fsample = fs)
    frequencies = getMeanFreq(n, config)

    tone_coefficients = cwt(signals.tone, config)
    chirp_coefficients = cwt(signals.chirp, config)
    tone_reconstruction = vec(icwt(tone_coefficients, config, PenroseDelta()))
    chirp_reconstruction = vec(icwt(chirp_coefficients, config, PenroseDelta()))

    complex_coefficients = cwt(complex.(signals.tone), config)
    float32_coefficients = cwt(Float32.(signals.tone), config)

    println("ContinuousWavelets.jl baseline")
    println("  package version: ", Base.pkgversion(ContinuousWavelets))
    println("  signal: n=$n, fs=$(fs) Hz")
    println("  configuration: Morlet(π), NoAve(), β=2")
    println("  coefficient layout: time × scale = ", size(tone_coefficients))
    println("  coefficient element type: ", eltype(tone_coefficients))
    @printf("  mean-frequency range: %.6f to %.6f Hz (%d scales; ascending=%s)\n",
        first(frequencies), last(frequencies), length(frequencies), string(issorted(frequencies)))
    @printf("  tone peak frequency: %.6f Hz\n", peak_frequency(tone_coefficients, frequencies))
    @printf("  tone PenroseDelta relative reconstruction error: %.6e\n",
        relative_error(signals.tone, tone_reconstruction))
    @printf("  chirp PenroseDelta relative reconstruction error: %.6e\n",
        relative_error(signals.chirp, chirp_reconstruction))
    println("  complex input element type: ", eltype(complex_coefficients))
    println("  Float32 input element type: ", eltype(float32_coefficients))
    @printf("  warm CWT time (tone): %.6f s\n", warm_seconds(() -> cwt(signals.tone, config)))

    return (; config, frequencies, tone_coefficients, chirp_coefficients)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_baseline()
end
