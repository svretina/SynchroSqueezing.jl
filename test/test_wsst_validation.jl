function dominant_frequency(result)
    rows = [argmax(abs2.(@view result.coefficients[:, column])) for column in axes(result.coefficients, 2)]
    return result.frequencies[rows]
end

@testset "Output-grid convergence" begin
    fs = 256.0
    n = 256
    frequency = 37.3
    _, signal = constant_tone(n, fs, frequency)
    coarse = wsst(signal; fs, freqrange = (30.0, 50.0), nfrequencies = 11)
    fine = wsst(signal; fs, freqrange = (30.0, 50.0), nfrequencies = 101)
    coarse_error = abs(dominant_frequency(coarse)[n ÷ 2] - frequency)
    fine_error = abs(dominant_frequency(fine)[n ÷ 2] - frequency)
    @test fine_error < coarse_error
    @test fine_error <= 0.1 + 10eps(Float64)
end

@testset "Analytical WSST signals" begin
    fs = 512.0
    n = 512
    interior = 65:448
    plan = plan_cwt(n; fs, voices_per_octave = 16)

    t, linear_signal = linear_chirp(n, fs, 35.0, 20.0)
    linear_result = wsst(
        plan,
        linear_signal;
        freqrange = (25.0, 70.0),
        nfrequencies = 181,
        threshold = RelativeThreshold(1e-5),
    )
    linear_exact = 35.0 .+ 20.0 .* t
    linear_estimate = dominant_frequency(linear_result)
    @test sqrt(sum(abs2, linear_estimate[interior] - linear_exact[interior]) / length(interior)) < 1.5

    _, quadratic_signal, quadratic_exact = quadratic_chirp(n, fs, 30.0, 15.0, 5.0)
    quadratic_result = wsst(
        plan,
        quadratic_signal;
        freqrange = (20.0, 70.0),
        nfrequencies = 201,
    )
    quadratic_estimate = dominant_frequency(quadratic_result)
    @test sqrt(sum(abs2, quadratic_estimate[interior] - quadratic_exact[interior]) / length(interior)) < 2.0

    _, exponential_signal, exponential_exact = exponential_chirp(n, fs, 30.0, 0.5)
    exponential_result = wsst(
        plan,
        exponential_signal;
        freqrange = (20.0, 60.0),
        nfrequencies = 161,
    )
    exponential_estimate = dominant_frequency(exponential_result)
    @test sqrt(sum(abs2, exponential_estimate[interior] - exponential_exact[interior]) / length(interior)) < 2.0

    _, am_signal, _ = amplitude_modulated_tone(n, fs, 48.0)
    am_result = wsst(plan, am_signal; freqrange = (35.0, 60.0), nfrequencies = 101)
    am_estimate = dominant_frequency(am_result)
    @test maximum(abs, am_estimate[interior] .- 48.0) <= 0.5
end

@testset "Robustness and deterministic conditions" begin
    fs = 512.0
    n = 512
    t = sample_times(n, fs)
    two_tones = cos.(2π * 40.0 .* t) .+ 0.7 .* cos.(2π * 96.0 .* t)
    result = wsst(two_tones; fs, freqrange = (20.0, 120.0), nfrequencies = 201)
    energy = vec(sum(abs2, result.coefficients; dims = 2))
    strongest = sortperm(energy; rev = true)[1:10]
    @test any(abs.(result.frequencies[strongest] .- 40.0) .<= 0.5)
    @test any(abs.(result.frequencies[strongest] .- 96.0) .<= 0.5)

    rng = MersenneTwister(20261020)
    noisy = two_tones .+ 0.05 .* randn(rng, n)
    noisy_result = wsst(noisy; fs, freqrange = (20.0, 120.0), nfrequencies = 201)
    @test all(isfinite, noisy_result.coefficients)

    for (length_, frequency) in ((511, 5.0), (512, 220.0))
        time = sample_times(length_, fs)
        signal = cos.(2π * frequency .* time)
        low = frequency < 10 ? 1.0 : 180.0
        high = frequency < 10 ? 15.0 : 250.0
        edge_result = wsst(signal; fs, freqrange = (low, high), nfrequencies = 141)
        @test size(edge_result.coefficients, 2) == length_
        @test all(isfinite, edge_result.coefficients)
    end

    scaled = wsst(1e-9 .* two_tones; fs, freqrange = (20.0, 120.0), nfrequencies = 201)
    @test scaled.coefficients ≈ 1e-9 .* result.coefficients rtol = 1e-5 atol = 1e-14
end
