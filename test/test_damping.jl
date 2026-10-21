using LinearAlgebra: norm

@testset "Analytic component features and damping fits" begin
    fs = 512.0
    n = 1024
    time = (0:n-1) ./ fs

    for (tau, phase, true_frequency) in (
        (0.7, 2π .* 55.0 .* time, fill(55.0, n)),
        (1.2, 2π .* (35.0 .* time .+ 8.0 .* time .^ 2), 35.0 .+ 16.0 .* time),
        (2.5, 2π .* (30.0 .* time .+ 5.0 .* time .^ 2 .+ time .^ 3),
         30.0 .+ 10.0 .* time .+ 3.0 .* time .^ 2),
    )
        exact = 1.4 .* exp.(-time ./ tau) .* exp.(im .* phase)
        features = component_features(exact, time)
        fit = fit_exponential_damping(features; window = (0.1, 1.8))
        nonlinear = fit_exponential_damping(
            features; method = :nonlinear, window = (0.1, 1.8),
        )
        @test fit.status == :ok
        @test fit.damping_time ≈ tau rtol = 1e-10
        @test nonlinear.damping_time ≈ tau rtol = 1e-8
        @test maximum(abs, features.instantaneous_frequency[3:end-2] .-
                           true_frequency[3:end-2]) < 0.1
        local_rate = local_damping_rate(features)
        @test maximum(abs, local_rate[3:end-2] .- inv(tau)) < 1e-10
    end

    real_tone = exp.(-time ./ 1.5) .* cos.(2π .* 48.0 .* time)
    real_features = component_features(real_tone, time)
    real_fit = fit_exponential_damping(real_features; window = (0.2, 1.7))
    @test abs(real_fit.damping_time / 1.5 - 1) < 0.02

    constant_fit = fit_exponential_damping(ones(n), time)
    growing_fit = fit_exponential_damping(exp.(time), time)
    @test constant_fit.status == :nondecaying
    @test growing_fit.status == :nondecaying
    @test isnan(constant_fit.damping_time)
    @test isnan(growing_fit.damping_time)

    amplitude = exp.(-time ./ 0.8)
    valid = trues(n)
    valid[1:100] .= false
    fit = fit_exponential_damping(
        amplitude, time;
        valid_mask = valid,
        window = 200:800,
        weights = range(0.5, 1.5; length = n),
        amplitude_floor = 0.01,
    )
    @test fit.damping_time ≈ 0.8 rtol = 1e-10
    @test !any(fit.selected[1:199])
    @test fit.fit_window == (time[200], time[800])
    @test fit.damping_time_standard_error < 1e-9
end

@testset "Analytic ridge reconstruction and damping sensitivity" begin
    fs = 256.0
    n = 512
    time = (0:n-1) ./ fs
    tau = 1.4
    truth = exp.(-time ./ tau) .* cos.(2π .* (42.0 .* time .+ 2.0 .* time .^ 2))
    true_frequency = 42.0 .+ 4.0 .* time
    result = wsst(
        truth;
        fs,
        freqrange = (30.0, 58.0),
        nfrequencies = 225,
        threshold = RelativeThreshold(1e-10),
        preserve_residual = false,
    )
    ridge = extract_ridge(
        result;
        freqrange = (38.0, 54.0),
        penalty = 1.0,
        continuity_scale = 1.0,
    )
    real_reconstruction = reconstruct_component(result, ridge; bandwidth_hz = 5.0)
    analytic_reconstruction = reconstruct_component(
        result, ridge; bandwidth_hz = 5.0, analytic = true,
    )
    @test eltype(analytic_reconstruction.signal) <: Complex
    @test real.(analytic_reconstruction.signal) ≈ real_reconstruction.signal
    @test ridge_metrics(ridge, true_frequency).rmse < 0.75

    features = component_features(
        analytic_reconstruction.signal,
        time;
        valid_mask = ridge.valid .& ridge.active,
    )
    fit = fit_exponential_damping(features; window = (0.3, 1.6))
    @test abs(fit.damping_time / tau - 1) < 0.08

    sensitivity = damping_sensitivity(
        result,
        ridge;
        bandwidths_hz = [4.0, 5.0, 6.0],
        windows = [(0.3, 1.5), (0.4, 1.6)],
    )
    @test length(sensitivity.cases) == 6
    @test sensitivity.minimum_damping_time <= sensitivity.maximum_damping_time
    @test sensitivity.damping_time_spread >= 0
end

@testset "Noisy and nonexponential envelope diagnostics" begin
    Random.seed!(404)
    fs = 200.0
    time = (0:399) ./ fs
    tau = 0.9
    clean = exp.(-time ./ tau)
    noisy = max.(clean .+ 0.025 .* randn(length(time)), 1e-6)
    log_fit = fit_exponential_damping(noisy, time; window = (0.1, 1.4))
    nonlinear_fit = fit_exponential_damping(
        noisy, time; method = :nonlinear, window = (0.1, 1.4),
    )
    @test abs(log_fit.damping_time / tau - 1) < 0.15
    @test abs(nonlinear_fit.damping_time / tau - 1) < 0.15
    @test isfinite(log_fit.damping_time_standard_error)
    @test isfinite(nonlinear_fit.damping_time_standard_error)

    modulated = clean .* (1 .+ 0.35 .* sin.(2π .* 3 .* time))
    modulated_fit = fit_exponential_damping(modulated, time; window = (0.1, 1.4))
    @test modulated_fit.log_rmse > log_fit.log_rmse

    colored_noise = zeros(length(time))
    innovations = randn(length(time))
    for index in 2:length(time)
        colored_noise[index] = 0.92colored_noise[index - 1] + innovations[index]
    end
    colored = max.(clean .+ 0.006 .* colored_noise, 1e-6)
    colored_fit = fit_exponential_damping(colored, time; window = (0.1, 1.4))
    @test abs(colored_fit.damping_time / tau - 1) < 0.15

    interference = abs.(
        clean .* exp.(2π * 30im .* time) .+
        0.3clean .* exp.(2π * 31.5im .* time)
    )
    interference_fit = fit_exponential_damping(interference, time; window = (0.1, 1.4))
    @test interference_fit.log_rmse > 0.15

    fading = component_features(
        clean .* exp.(2π * 35im .* time), time; amplitude_floor = 0.2,
    )
    fading_fit = fit_exponential_damping(fading; window = (0.0, 1.8))
    @test !any(fading_fit.selected[fading.amplitude .<= 0.2])
    @test fading_fit.damping_time ≈ tau rtol = 1e-10
end
