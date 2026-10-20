#!/usr/bin/env julia

using BenchmarkTools
using ContinuousWavelet
using ContinuousWavelets
using LinearAlgebra
using Printf

include("support.jl")
using .BackendExperimentSupport

function continuouswavelets_configuration(n::Int, fs::Float64)
    config = ContinuousWavelets.wavelet(
        ContinuousWavelets.Morlet(π);
        averagingType = ContinuousWavelets.NoAve(), β = 2, fsample = fs,
    )
    return config, ContinuousWavelets.getMeanFreq(n, config)
end

function trial_summary(f)
    f() # compile and warm up
    trial = @benchmark $f() samples = 10 evals = 1
    estimate = median(trial)
    return (; seconds = estimate.time / 1e9, bytes = estimate.memory, allocations = estimate.allocs)
end

"""
Compare nominally matched Morlet CWTs.

`ContinuousWavelet.jl` documents unit center frequency for its Morlet wavelet,
so `fs ./ frequency` supplies a nominal scale match. This is not an assertion
of equivalent normalization or exact wavelet shape; those are reported rather
than hidden by the comparison.
"""
function compare_backends(; n::Int = 2048, fs::Float64 = 1024.0)
    signals = analytical_signals(n, fs)
    cw_config_seconds = @elapsed cw_config, frequencies = continuouswavelets_configuration(n, fs)
    cw_wavelet_seconds = @elapsed ContinuousWavelets.computeWavelets(n, cw_config)
    cw_coefficients = ContinuousWavelets.cwt(signals.tone, cw_config)
    cw_inverse = vec(ContinuousWavelets.icwt(
        cw_coefficients, cw_config, ContinuousWavelets.PenroseDelta(),
    ))

    # The alternate backend sorts scales internally and accepts Float64 real inputs only.
    nominal_scales = sort(fs ./ frequencies)
    alternate_plan_seconds = @elapsed alternate_plan = ContinuousWaveletTransform(
        MorletWavelet(1.0), nominal_scales,
    )
    alternate_coefficients = transform(alternate_plan, signals.tone)
    alternate_inverse = synthesis(alternate_coefficients)

    cw_timing = trial_summary(() -> ContinuousWavelets.cwt(signals.tone, cw_config))
    alternate_timing = trial_summary(() -> transform(alternate_plan, signals.tone))

    println("CWT backend comparison (nominal scale matching only)")
    println("  signal: n=$n, fs=$(fs) Hz")
    println("  ContinuousWavelets.jl version: ", Base.pkgversion(ContinuousWavelets))
    println("  ContinuousWavelet.jl version: ", Base.pkgversion(ContinuousWavelet))
    println("  ContinuousWavelets layout: ", size(cw_coefficients), ", ", eltype(cw_coefficients))
    println("  ContinuousWavelet layout: ", size(alternate_coefficients), ", ", eltype(alternate_coefficients))
    @printf("  ContinuousWavelets configuration construction: %.6f s\n", cw_config_seconds)
    @printf("  ContinuousWavelets wavelet-kernel construction: %.6f s\n", cw_wavelet_seconds)
    @printf("  ContinuousWavelet plan construction: %.6f s\n", alternate_plan_seconds)
    @printf("  ContinuousWavelets PenroseDelta reconstruction error: %.6e\n",
        relative_error(signals.tone, cw_inverse))
    @printf("  ContinuousWavelet synthesis reconstruction error: %.6e\n",
        relative_error(signals.tone, alternate_inverse))
    @printf("  ContinuousWavelets warm median: %.6f s, %d bytes, %d allocations\n",
        cw_timing.seconds, cw_timing.bytes, cw_timing.allocations)
    @printf("  ContinuousWavelet warm median: %.6f s, %d bytes, %d allocations\n",
        alternate_timing.seconds, alternate_timing.bytes, alternate_timing.allocations)

    return (; cw_timing, alternate_timing, cw_config_seconds, cw_wavelet_seconds,
        alternate_plan_seconds, cw_error = relative_error(signals.tone, cw_inverse),
        alternate_error = relative_error(signals.tone, alternate_inverse))
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_backends()
end
