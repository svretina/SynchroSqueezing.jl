@testset "Operational mode support and persistent termination" begin
    fs = 100.0
    time = (0:399) ./ fs
    tau = 0.32
    component = exp.(-time ./ tau) .* exp.(2π * 12im .* time)
    features = component_features(component, time)

    support = mode_support(
        features;
        relative_amplitude_threshold = 0.01,
        persistence_seconds = 0.05,
    )
    predicted = tau * log(100)
    @test abs(support.stop_time - predicted) <= 2 / fs
    @test support.last_supported_time < support.stop_time
    @test support.stop_reason == :relative_amplitude
    @test all(support.supported[1:support.last_supported_index])
    @test !any(support.supported[findfirst(.!support.supported):end])
    @test support.config.reference_amplitude ≈ 1.0

    raw = ComponentReconstruction(component, trues(1, length(time)),
                                  ones(1, length(time)))
    usable = supported_reconstruction(raw, support)
    @test usable.raw === raw
    @test usable.signal[support.supported] == component[support.supported]
    @test all(iszero, usable.signal[.!support.supported])

    dipped = copy(component)
    dipped[101:103] .*= 1e-4
    dipped_features = component_features(dipped, time)
    persistent = mode_support(
        dipped_features;
        relative_amplitude_threshold = 0.01,
        persistence_seconds = 0.05,
        recovery_factor = 1.2,
    )
    @test persistent.stop_time > time[103]
    @test any(persistent.uncertain[101:103])
    @test !any(persistent.fit_mask[101:103])

    prominence = RidgeProminence(
        fill(10.0, length(time)), ones(length(time)), fill(0.1, length(time)),
        trues(length(time)), 1.0, (2.0, 4.0),
    )
    prominence.ratio[121:end] .= 1.0
    combined = mode_support(
        features;
        prominence,
        relative_amplitude_threshold = 1e-3,
        prominence_threshold = 3.0,
        persistence_seconds = 0.04,
    )
    @test combined.stop_time == time[121]
    @test combined.last_supported_index == 120
    @test combined.stop_reason == :insufficient_prominence

    prominence_only = mode_support(
        features;
        prominence,
        relative_amplitude_threshold = nothing,
        prominence_threshold = 3.0,
        persistence_seconds = 0.04,
    )
    @test prominence_only.stop_time == time[121]
    @test prominence_only.config.relative_amplitude_threshold === nothing

    boundary = trues(length(time))
    boundary[1:10] .= false
    boundary[end-9:end] .= false
    bounded = mode_support(
        features;
        relative_amplitude_threshold = 1e-6,
        persistence_seconds = 0,
        boundary_mask = boundary,
    )
    @test findfirst(bounded.fit_mask) == 11
    @test bounded.stop_time == time[end-9]
    @test bounded.stop_reason == :boundary

    @test_throws ArgumentError mode_support(
        features; prominence_threshold = 2.0,
    )
    @test_throws ArgumentError mode_support(
        features; relative_amplitude_threshold = 0,
    )
    @test_throws ArgumentError mode_support(
        features; relative_amplitude_threshold = nothing,
    )

    end_supported = mode_support(
        features; relative_amplitude_threshold = 1e-8,
    )
    @test end_supported.stop_time === nothing
    @test end_supported.stop_reason == :end_of_data

    sensitivity = mode_support_sensitivity(
        features;
        relative_amplitude_thresholds = [0.02, 0.05, 0.10],
        persistence_seconds = 0.03,
    )
    @test length(sensitivity.cases) == 3
    @test sensitivity.damping_time_spread < 1e-10
    @test all(case -> case.damping_fit.damping_time ≈ tau,
              sensitivity.cases)
end

@testset "Nonexponential termination follows the observed envelope" begin
    fs = 200.0
    time = (0:799) ./ fs
    tau1 = 0.7
    tau2 = 0.9
    threshold = 0.05
    amplitude = exp.(-time ./ tau1 .- time .^ 2 ./ (2tau2^2))
    component = amplitude .* exp.(2π * 25im .* time)
    features = component_features(component, time)
    support = mode_support(
        features;
        relative_amplitude_threshold = threshold,
        persistence_seconds = 0.025,
    )
    expected = tau2^2 * (
        -inv(tau1) + sqrt(inv(tau1)^2 + 2log(inv(threshold)) / tau2^2)
    )
    @test abs(support.stop_time - expected) <= 2 / fs
    @test support.stop_reason == :relative_amplitude
    rate = local_damping_rate(features)
    @test rate[500] > rate[100]
end

@testset "Local WSST ridge prominence" begin
    frequencies = collect(0.0:1.0:20.0)
    times = collect(0.0:0.1:0.3)
    coefficients = fill(1.0 + 0im, length(frequencies), length(times))
    coefficients[10:12, :] .= 4.0 + 0im
    convention = ReconstructionConvention(
        ones(length(frequencies)), 1.0 + 0im, ones(3), 0.0, (0.0, 20.0),
        :standard_analytic_l1, :l1, :analytic_complex,
    )
    diagnostics = WSSTDiagnostics(
        zeros(Int, 4), zeros(Int, 4), zeros(Int, 4), zeros(Int, 4), zeros(Int, 4),
        zeros(Int, 4), zeros(Int, 4), zeros(ComplexF64, 4), zeros(ComplexF64, 4),
    )
    metadata = (; sample_rate = 10.0, sample_interval = 0.1)
    result = WSSTResult(
        coefficients, frequencies, collect(-0.5:1.0:20.5), :linear,
        ones(length(frequencies)), :nearest, NoThreshold(),
        ReassignmentMetadata(1, :phase_derivative, :reference, NaN, :none, Inf,
                             0.0, 0.0, 0.0, nothing),
        diagnostics, convention,
        nothing, zeros(ComplexF64, 4), metadata, nothing, nothing, nothing, nothing,
    )
    config = RidgeConfig(activity_threshold = 0.0)
    ridge = RidgeResult(
        fill(11, 4), fill(10.0, 4), times, trues(4), trues(4), falses(4),
        fill(16.0, 4), 0.0, (0.0, 20.0), config,
    )
    prominence = ridge_prominence(
        result, ridge; ridge_bandwidth_hz = 1.0, background_band_hz = (2.0, 4.0),
    )
    @test all(prominence.valid)
    @test prominence.ratio ≈ fill(16.0, 4)
    @test_throws ArgumentError ridge_prominence(
        result, ridge; ridge_bandwidth_hz = 2.0, background_band_hz = (1.0, 4.0),
    )
end
