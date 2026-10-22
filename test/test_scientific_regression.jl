using LinearAlgebra: norm
using Statistics: mean, median

include("fixtures/scientific_reference_05a.jl")

function _fnv1a_bits(values)
    return foldl(values; init = UInt64(1469598103934665603)) do state, value
        (state ⊻ UInt64(value)) * UInt64(1099511628211)
    end
end

@testset "Milestone 05A frozen scientific pipeline" begin
    reference = SCIENTIFIC_REFERENCE_05A
    config = reference.configuration
    n = config.samples
    fs = config.sample_rate_hz
    times = collect(0:n-1) ./ fs
    damping_times = (0.4, 1.8)
    truths = (
        exp.(-times ./ damping_times[1]) .* exp.(2π * 35im .* times),
        0.7 .* exp.(-times ./ damping_times[2]) .* exp.(2π * 85im .* times),
    )
    signal = real.(truths[1] .+ truths[2])

    plan = plan_cwt(n; fs, voices_per_octave = config.voices_per_octave)
    transformed = cwt(plan, signal)
    phase = instantaneous_frequency(
        transformed;
        threshold = RelativeThreshold(1e-8),
    )
    wsst_options = (
        freqrange = config.wsst_freqrange_hz,
        nfrequencies = config.wsst_frequency_bins,
        threshold = RelativeThreshold(1e-8),
        preserve_residual = true,
    )
    optimized = wsst(transformed; wsst_options...)
    scalar = wsst(
        transformed;
        wsst_options...,
        implementation = :reference,
        keep_phase = true,
    )

    @test optimized.coefficients == scalar.coefficients
    @test optimized.residual == scalar.residual
    @test isequal(scalar.phase, phase.values)
    @test scalar.phase_valid == phase.valid
    for field in propertynames(optimized.diagnostics)
        @test getproperty(optimized.diagnostics, field) ==
              getproperty(scalar.diagnostics, field)
    end

    @test sum(abs2, transformed.coefficients) ≈ reference.cwt.sum_abs2 rtol = 2e-12
    @test maximum(abs, transformed.coefficients) ≈ reference.cwt.maximum_abs rtol = 2e-12
    @test transformed.coefficients[10, 100] ≈ reference.cwt.anchor rtol = 2e-12 atol = 1e-15
    @test sum(abs2, transformed.derivative) ≈ reference.derivative.sum_abs2 rtol = 2e-12
    @test maximum(abs, transformed.derivative) ≈ reference.derivative.maximum_abs rtol = 2e-12
    @test transformed.derivative[10, 100] ≈ reference.derivative.anchor rtol = 2e-12 atol = 1e-15
    @test count(phase.valid) == reference.phase.valid_count
    @test _fnv1a_bits(phase.valid) == reference.phase.valid_code
    @test sum(abs2, optimized.coefficients) ≈ reference.wsst.sum_abs2 rtol = 2e-12
    @test maximum(abs, optimized.coefficients) ≈ reference.wsst.maximum_abs rtol = 2e-12
    @test optimized.coefficients[61, 100] ≈ reference.wsst.anchor rtol = 2e-12 atol = 1e-15

    ridge_options = (
        penalty = 0.5,
        continuity_scale = 0.5,
        max_rate = 20.0,
        activity_threshold = 0.0,
    )
    ridges = extract_ridges(
        optimized;
        freqranges = config.ridge_freqranges_hz,
        ridge_options...,
    )
    dense_ridges = [extract_ridge(
        optimized;
        freqrange,
        method = :dense,
        ridge_options...,
    ) for freqrange in config.ridge_freqranges_hz]
    components = reconstruct_components(
        optimized,
        ridges;
        bandwidth_hz = config.reconstruction_bandwidth_hz,
        ownership = :exclusive,
        analytic = true,
    )
    features = [component_features(
        components.components[index],
        times;
        valid_mask = ridges[index].valid .& ridges[index].active,
    ) for index in eachindex(ridges)]
    supports = [mode_support(
        features[index];
        relative_amplitude_threshold = config.amplitude_threshold,
        persistence_seconds = 0.02,
        boundary_mask = .!ridges[index].at_frequency_boundary,
    ) for index in eachindex(ridges)]
    fits = [fit_exponential_damping(
        features[index];
        valid_mask = supports[index].fit_mask,
    ) for index in eachindex(features)]
    rates = [local_damping_rate(feature) for feature in features]

    for index in eachindex(ridges)
        expected = reference.modes[index]
        ridge = ridges[index]
        @test ridge.indices == dense_ridges[index].indices
        @test ridge.total_cost == dense_ridges[index].total_cost
        @test _fnv1a_bits(ridge.indices) == expected.ridge_code
        @test count(ridge.active) == expected.active_count
        @test _fnv1a_bits(ridge.active) == expected.active_code
        @test count(ridge.at_frequency_boundary) == expected.boundary_count
        @test _fnv1a_bits(ridge.at_frequency_boundary) == expected.boundary_code
        target_frequency = index == 1 ? 35.0 : 85.0
        @test sqrt(mean(abs2, ridge.frequencies .- target_frequency)) ≈
              expected.ridge_rmse_hz atol = 1e-14

        relative_error = norm(components.components[index] - truths[index]) / norm(truths[index])
        @test relative_error ≈ expected.reconstruction_relative_error rtol = 2e-11
        @test sum(features[index].amplitude) ≈ expected.envelope_sum rtol = 2e-11
        @test maximum(features[index].amplitude) ≈ expected.envelope_maximum rtol = 2e-11
        @test count(supports[index].supported) == expected.support_count
        @test _fnv1a_bits(supports[index].supported) == expected.support_code
        @test supports[index].stop_time == expected.stop_time
        @test supports[index].stop_reason == expected.stop_reason
        @test fits[index].damping_time ≈ expected.damping_time rtol = 2e-10
        finite_rates = filter(isfinite, rates[index])
        @test median(finite_rates) ≈ expected.local_rate_median rtol = 2e-10
    end

    ownership_sum = copy(components.residual_weights)
    for weights in components.weights
        ownership_sum .+= weights
    end
    ownership_error = maximum(abs, ownership_sum .- 1)
    represented = inverse_wsst(optimized)
    conservation_error = norm(
        real.(sum(components.components) + components.residual) - represented,
    ) / norm(represented)
    @test ownership_error == reference.ownership_error
    # This diagnostic is itself at floating-point roundoff; compare it on an
    # absolute epsilon scale rather than demanding stable relative digits from
    # two values of order 1e-16.
    @test conservation_error ≈ reference.residual_conservation_error atol = 5eps(Float64)
    @test norm(components.residual) ≈ reference.residual_norm rtol = 2e-11
end
