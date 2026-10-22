using LinearAlgebra: norm
using Statistics: mean, median

function _dominant_wsst_frequency(result)
    rows = [argmax(abs2.(@view result.coefficients[:, column]))
            for column in axes(result.coefficients, 2)]
    return result.frequencies[rows]
end

function _median_wsst_width(result, truth, columns)
    widths = Float64[]
    for column in columns
        energy = abs2.(@view result.coefficients[:, column])
        total = sum(energy)
        total > 0 || continue
        push!(widths, sqrt(sum(energy .* (result.frequencies .- truth[column]) .^ 2) / total))
    end
    return median(widths)
end

@testset "Second-order WSST analytical operators" begin
    for T in (Float32, Float64), rate in (T(60), T(-60))
        n = 512
        fs = T(512)
        time = T.(0:n-1) ./ fs
        initial = rate > 0 ? T(40) : T(150)
        signal = cos.(T(2π) .* (initial .* time .+ T(0.5) .* rate .* time .^ 2))
        plan = plan_cwt(n; fs, voices_per_octave = 32, precision = T)
        operators = second_order_operators(plan, signal)
        estimate = second_order_frequency(
            operators; threshold = RelativeThreshold(T(1e-4)),
        )

        chirp_rates = T[]
        frequency_errors = T[]
        for column in 100:412
            row = argmax(abs.(@view operators.cwt.coefficients[:, column]))
            estimate.second_order_valid[row, column] || continue
            push!(chirp_rates, estimate.chirp_rate[row, column])
            truth = initial + rate * time[column]
            push!(frequency_errors, estimate.values[row, column] - truth)
        end
        @test length(chirp_rates) > 300
        tolerance = T === Float32 ? 2e-3 : 1e-8
        @test abs(median(chirp_rates) - rate) < tolerance
        @test sqrt(mean(abs2, frequency_errors)) < tolerance
    end

    n = 512
    fs = 512.0
    time = (0:n-1) ./ fs
    tone = cos.(2π * 72.0 .* time)
    plan = plan_cwt(n; fs, voices_per_octave = 24)
    operators = second_order_operators(plan, tone)
    estimate = second_order_frequency(operators; threshold = RelativeThreshold(1e-5))
    tone_rates = Float64[]
    for column in 100:412
        row = argmax(abs.(@view operators.cwt.coefficients[:, column]))
        estimate.second_order_valid[row, column] || continue
        push!(tone_rates, estimate.chirp_rate[row, column])
    end
    @test maximum(abs, tone_rates) < 1e-8
end

@testset "WSST2 nonlinear and amplitude-varying signals" begin
    n = 384
    fs = 384.0
    time = (0:n-1) ./ fs
    plan = plan_cwt(n; fs, voices_per_octave = 24)
    cases = NamedTuple[]

    push!(cases, (
        name = :damped_tone,
        signal = exp.(-0.8 .* time) .* cos.(2π * 70.0 .* time),
        frequency = fill(70.0, n),
        omitted = 1:0,
    ))
    quadratic_frequency = 35.0 .+ 30.0 .* time .+ 20.0 .* time .^ 2
    quadratic_phase = 2π .* (35.0 .* time .+ 15.0 .* time .^ 2 .+
                             (20.0 / 3) .* time .^ 3)
    push!(cases, (
        name = :quadratic_chirp,
        signal = cos.(quadratic_phase),
        frequency = quadratic_frequency,
        omitted = 1:0,
    ))
    exponential_frequency = 30.0 .* exp.(0.8 .* time)
    exponential_phase = 2π * 30.0 / 0.8 .* (exp.(0.8 .* time) .- 1)
    push!(cases, (
        name = :exponential_chirp,
        signal = cos.(exponential_phase),
        frequency = exponential_frequency,
        omitted = 1:0,
    ))
    linear_frequency = 40.0 .+ 50.0 .* time
    linear_phase = 2π .* (40.0 .* time .+ 25.0 .* time .^ 2)
    push!(cases, (
        name = :amplitude_modulated_chirp,
        signal = (1 .+ 0.3 .* cos.(2π * 3.0 .* time)) .* cos.(linear_phase),
        frequency = linear_frequency,
        omitted = 1:0,
    ))
    push!(cases, (
        name = :damped_chirp,
        signal = exp.(-1.2 .* time) .* cos.(linear_phase),
        frequency = linear_frequency,
        omitted = 1:0,
    ))
    abrupt_amplitude = ifelse.(time .< 0.55, 1.0, 0.2)
    push!(cases, (
        name = :abrupt_amplitude,
        signal = abrupt_amplitude .* cos.(linear_phase),
        frequency = linear_frequency,
        omitted = 195:230,
    ))

    for case in cases
        operators = second_order_operators(plan, case.signal)
        estimate = second_order_frequency(
            operators;
            threshold = RelativeThreshold(1e-4),
        )
        errors = Float64[]
        for column in 80:304
            column in case.omitted && continue
            row = argmax(abs.(@view operators.cwt.coefficients[:, column]))
            estimate.second_order_valid[row, column] || continue
            push!(errors, estimate.values[row, column] - case.frequency[column])
        end
        @test length(errors) > 180
        @test sqrt(mean(abs2, errors)) < 0.15
    end
end

@testset "WSST2 component reconstruction and damping" begin
    n = 1024
    fs = 512.0
    time = (0:n-1) ./ fs
    damping_time = 1.2
    frequency = 40.0 .+ 50.0 .* time
    phase = 2π .* (40.0 .* time .+ 25.0 .* time .^ 2)
    signal = exp.(-time ./ damping_time) .* cos.(phase)
    plan = plan_cwt(n; fs, voices_per_octave = 32)
    result = wsst(
        plan,
        signal;
        order = 2,
        freqrange = (25.0, 160.0),
        nfrequencies = 700,
        threshold = RelativeThreshold(1e-7),
    )
    ridge = extract_ridge(result; freqrange = (35.0, 145.0), method = :maximum)
    component = reconstruct_component(result, ridge; bandwidth_hz = 6.0, analytic = true)
    features = component_features(component.signal, time; valid_mask = ridge.active)
    fit = fit_exponential_damping(features; window = (0.3, 1.5))
    interior = falses(n)
    interior[154:768] .= true

    @test ridge_metrics(ridge, frequency; expected = interior).rmse < 0.1
    @test component_metrics(real.(component.signal[interior]), signal[interior]).relative_error < 1e-3
    @test abs(fit.damping_time / damping_time - 1) < 0.01
end

@testset "WSST2 concentration, metadata, and reconstruction" begin
    n = 1024
    fs = 512.0
    time = (0:n-1) ./ fs
    rate = 80.0
    truth = 35.0 .+ rate .* time
    signal = cos.(2π .* (35.0 .* time .+ 0.5rate .* time .^ 2))
    plan = plan_cwt(n; fs, voices_per_octave = 32)
    common = (
        freqrange = (20.0, 210.0),
        nfrequencies = 900,
        threshold = RelativeThreshold(1e-5),
        preserve_residual = true,
    )
    first = wsst(plan, signal; order = 1, common...)
    second = wsst(plan, signal; order = 2, keep_phase = true, common...)
    interior = 200:824

    @test first.reassignment_metadata.order == 1
    @test second.reassignment_metadata.order == 2
    @test second.reassignment_metadata.estimator === :time_derivative_wsst2
    @test second.reassignment_metadata.implementation === :reference
    @test 0.7 < second.reassignment_metadata.second_order_valid_fraction <= 1
    @test second.reassignment_metadata.second_order_valid_fraction +
          second.reassignment_metadata.fallback_fraction +
          second.reassignment_metadata.excluded_fraction ≈ 1
    @test second.phase !== nothing
    @test second.diagnostics.accepted_sum ≈ second.diagnostics.deposited_sum atol = 1e-12
    @test _median_wsst_width(second, truth, interior) <
          0.85 * _median_wsst_width(first, truth, interior)
    @test sqrt(mean(abs2, _dominant_wsst_frequency(second)[interior] - truth[interior])) < 0.15

    reconstructed = inverse_wsst(second; mode = :full, method = :frame)
    @test norm(reconstructed - signal) / norm(signal) < 1e-10

    ridge = extract_ridge(
        second;
        freqrange = (50.0, 175.0),
        penalty = 1.0,
        continuity_scale = 1.0,
        max_rate = 120.0,
    )
    component = reconstruct_component(second, ridge; bandwidth_hz = 4.0, analytic = true)
    features = component_features(component.signal, time; valid_mask = ridge.active)
    expected = falses(n)
    expected[interior] .= true
    @test ridge_metrics(ridge, truth; expected).rmse < 0.2
    @test count(features.valid) > n ÷ 2
end

@testset "WSST2 stability and fallback" begin
    n = 256
    fs = 256.0
    time = (0:n-1) ./ fs
    signal = cos.(2π .* (30.0 .* time .+ 20.0 .* time .^ 2))
    plan = plan_cwt(n; fs, voices_per_octave = 16)
    operators = second_order_operators(plan, signal)

    fallback = second_order_frequency(
        operators;
        threshold = RelativeThreshold(1e-5),
        denominator_threshold = 1e6,
        fallback = :first_order,
    )
    excluded = second_order_frequency(
        operators;
        threshold = RelativeThreshold(1e-5),
        denominator_threshold = 1e6,
        fallback = :exclude,
    )
    @test !any(fallback.second_order_valid)
    @test count(fallback.used_fallback) == count(fallback.valid)
    @test !any(excluded.valid)

    excluded_result = wsst(
        operators;
        denominator_threshold = 1e6,
        fallback = :exclude,
        freqrange = (10.0, 100.0),
        nfrequencies = 181,
    )
    @test excluded_result.reassignment_metadata.excluded_fraction ≈ 1

    zero_result = wsst(
        plan,
        zeros(n);
        order = 2,
        freqrange = (10.0, 100.0),
        nfrequencies = 181,
    )
    @test all(iszero, zero_result.coefficients)
    @test all(==(0), zero_result.diagnostics.accepted)

    limited = second_order_frequency(
        operators;
        threshold = RelativeThreshold(1e-5),
        max_abs_chirp_rate = 10.0,
        fallback = :first_order,
    )
    @test count(limited.used_fallback) > 0
    @test all(isfinite, limited.values[limited.valid])

    cwt_result = cwt(plan, signal)
    @test_throws ArgumentError wsst(cwt_result; order = 2)
    @test_throws ArgumentError wsst(operators; order = 1)
    @test_throws ArgumentError second_order_frequency(operators; fallback = :invalid)
    @test_throws ArgumentError second_order_frequency(operators; denominator_threshold = -1)
end
