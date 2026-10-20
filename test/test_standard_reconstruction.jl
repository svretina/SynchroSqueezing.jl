import ContinuousWavelets
using LinearAlgebra: norm

frequency_mask(result, curve, half_width) =
    abs.(reshape(result.frequencies, :, 1) .- reshape(curve, 1, :)) .<= half_width

@testset "Morse admissibility and finite synthesis response" begin
    wavelet = ContinuousWavelets.Morse(3.0, 20.0, 1.0)
    analytical = admissibility_constant(wavelet)
    numerical = admissibility_constant(wavelet; method = :quadrature)
    @test analytical > 0
    @test abs(analytical - numerical) / analytical < 1e-12
    @test_throws ArgumentError admissibility_constant(
        ContinuousWavelets.Morse(3.0, 0.0, 1.0),
    )

    errors = Float64[]
    for voices in (4, 8, 16)
        plan = plan_cwt(1024; fs = 1024.0, voices_per_octave = voices)
        @test plan.standard_compatible
        push!(errors, abs(plan.synthesis_response[65] / plan.admissibility - 1))
        first_reliable, last_reliable = plan.reliable_frequency_range
        @test first_reliable <= 64.0 <= last_reliable
    end
    @test errors[3] < errors[2] < errors[1]
    @test errors[end] < 1e-12

    coverage_lower_bounds = [
        plan_cwt(samples; fs = 1024.0).reliable_frequency_range[1]
        for samples in (128, 256, 512, 1024)
    ]
    @test coverage_lower_bounds == [16.0, 8.0, 4.0, 2.0]
    @test issorted(reverse(coverage_lower_bounds))
end

@testset "Direct CWT scale sum and standard SST inverse" begin
    fs = 1024.0
    n = 1024
    time = sample_times(n, fs)
    real_tone = cos.(2π * 64.0 .* time)
    plan = plan_cwt(n; fs)
    transformed = cwt(plan, real_tone)
    integrated = vec(sum(
        transformed.coefficients .* reshape(transformed.squeezing_weights, :, 1);
        dims = 1,
    ))
    direct = 2 .* real.(integrated ./ transformed.admissibility)
    @test norm(direct - real_tone) / norm(real_tone) < 1e-3

    squeezed = wsst(
        transformed;
        freqrange = (40.0, 80.0),
        nfrequencies = 161,
        threshold = NoThreshold(),
        preserve_residual = false,
    )
    @test squeezed.cwt === nothing
    @test squeezed.derivative === nothing
    @test maximum(abs, vec(sum(squeezed.coefficients; dims = 1)) - integrated) < 1e-10
    @test norm(inverse_wsst(squeezed) - real_tone) / norm(real_tone) < 1e-3

    analytic_tone = cis.(2π * 64.0 .* time)
    analytic_result = wsst(
        plan,
        analytic_tone;
        freqrange = (40.0, 80.0),
        nfrequencies = 161,
        threshold = NoThreshold(),
        preserve_residual = false,
    )
    reconstructed = inverse_wsst(analytic_result)
    @test eltype(reconstructed) == ComplexF64
    @test norm(reconstructed - analytic_tone) / norm(analytic_tone) < 1e-3
    @test norm(reconstructed - 2 .* analytic_tone) / norm(analytic_tone) > 0.9
end

@testset "Known two-component reconstruction" begin
    fs = 1024.0
    n = 1024
    time = sample_times(n, fs)
    component1 = cos.(2π * 40.0 .* time)
    component2 = 0.7 .* cos.(2π * 96.0 .* time)
    signal = component1 + component2
    result = wsst(
        signal;
        fs,
        freqrange = (20.0, 120.0),
        nfrequencies = 401,
        threshold = RelativeThreshold(1e-10),
        preserve_residual = false,
    )
    mask1 = frequency_mask(result, fill(40.0, n), 4.0)
    mask2 = frequency_mask(result, fill(96.0, n), 4.0)
    recovered1 = reconstruct_region(result, mask1)
    recovered2 = reconstruct_region(result, mask2)
    error1 = norm(recovered1 - component1) / norm(component1)
    error2 = norm(recovered2 - component2) / norm(component2)
    combined_error = norm(signal - recovered1 - recovered2) / norm(signal)
    @test error1 < 0.01
    @test error2 < 0.01
    @test combined_error < 0.01
    phasor1 = 2sum(recovered1 .* cis.(-2π * 40.0 .* time)) / n
    phasor2 = 2sum(recovered2 .* cis.(-2π * 96.0 .* time)) / n
    @test abs(abs(phasor1) - 1) < 0.01
    @test abs(abs(phasor2) - 0.7) / 0.7 < 0.01
    @test abs(angle(phasor1)) < 0.01
    @test abs(angle(phasor2)) < 0.01
    contamination = abs(2sum(recovered1 .* cis.(-2π * 96.0 .* time)) / n)
    @test contamination < 0.01

    bandwidth_errors = [
        norm(reconstruct_region(result, frequency_mask(result, fill(96.0, n), width)) - component2) /
        norm(component2) for width in (2.0, 4.0, 8.0)
    ]
    @test maximum(bandwidth_errors) < 0.01

    unselected = .!(mask1 .| mask2)
    sst_residual = reconstruct_region(result, unselected)
    @test norm(recovered1 + recovered2 + sst_residual - inverse_wsst(result)) /
          norm(signal) < 1e-12
    @test_throws DimensionMismatch reconstruct_region(result, falses(2, 2))
    @test_throws ArgumentError reconstruct_region(result, mask1; method = :frame)

    analytic1 = cis.(2π * 40.0 .* time)
    analytic2 = 0.7 .* cis.(2π * 96.0 .* time)
    analytic_result = wsst(
        analytic1 + analytic2;
        fs,
        freqrange = (20.0, 120.0),
        nfrequencies = 401,
        threshold = RelativeThreshold(1e-10),
        preserve_residual = false,
    )
    analytic_recovered1 = reconstruct_region(analytic_result, mask1)
    analytic_recovered2 = reconstruct_region(analytic_result, mask2)
    @test norm(analytic_recovered1 - analytic1) / norm(analytic1) < 0.01
    @test norm(analytic_recovered2 - analytic2) / norm(analytic2) < 0.01
end

@testset "Known three-component reconstruction" begin
    fs = 1024.0
    n = 1024
    time = sample_times(n, fs)
    frequencies = (20.0, 80.0, 200.0)
    amplitudes = (1.0, 0.7, 0.5)
    components = [amplitude .* cos.(2π * frequency .* time)
                  for (amplitude, frequency) in zip(amplitudes, frequencies)]
    result = wsst(
        sum(components);
        fs,
        freqrange = (10.0, 250.0),
        nfrequencies = 481,
        threshold = RelativeThreshold(1e-10),
        preserve_residual = false,
    )
    recovered = [
        reconstruct_region(result, frequency_mask(result, fill(frequency, n), 6.0))
        for frequency in frequencies
    ]
    errors = [norm(estimate - exact) / norm(exact)
              for (estimate, exact) in zip(recovered, components)]
    @test maximum(errors) < 0.01
    @test norm(sum(recovered) - sum(components)) / norm(sum(components)) < 0.01
end

@testset "Legacy frame inverse is explicitly separate" begin
    fs = 256.0
    n = 256
    _, signal = constant_tone(n, fs, 32.0)
    legacy = cwt(plan_cwt(n; fs, standard = false), signal)
    result = wsst(legacy; threshold = NoThreshold())
    @test !legacy.standard_compatible
    @test_throws ArgumentError inverse_wsst(result; method = :standard)
    @test norm(inverse_wsst(result; method = :frame, mode = :full) - signal) /
          norm(signal) < 1e-8
end

@testset "Known parallel chirp reconstruction" begin
    fs = 512.0
    n = 1024
    time = sample_times(n, fs)
    frequency1 = 35.0 .+ 8.0 .* time
    frequency2 = 95.0 .+ 5.0 .* time
    component1 = cos.(2π .* (35.0 .* time .+ 4.0 .* time .^ 2))
    component2 = 0.6 .* cos.(2π .* (95.0 .* time .+ 2.5 .* time .^ 2))
    result = wsst(
        component1 + component2;
        fs,
        freqrange = (20.0, 130.0),
        nfrequencies = 441,
        threshold = RelativeThreshold(1e-8),
        preserve_residual = false,
    )
    recovered1 = reconstruct_region(result, frequency_mask(result, frequency1, 5.0))
    recovered2 = reconstruct_region(result, frequency_mask(result, frequency2, 5.0))
    interior = 129:896
    @test norm(recovered1[interior] - component1[interior]) /
          norm(component1[interior]) < 0.03
    @test norm(recovered2[interior] - component2[interior]) /
          norm(component2[interior]) < 0.03
end

@testset "Quadratic, exponential, AM, and approaching components" begin
    fs = 512.0
    n = 1024
    time, component1, frequency1 = quadratic_chirp(n, fs, 30.0, 5.0, 2.0)
    component1 .*= 1 .+ 0.2 .* cos.(2π * 2.0 .* time)
    _, component2, frequency2 = exponential_chirp(n, fs, 100.0, 0.1)
    component2 .*= 0.6
    result = wsst(
        component1 + component2;
        fs,
        freqrange = (20.0, 140.0),
        nfrequencies = 481,
        threshold = RelativeThreshold(1e-8),
        preserve_residual = false,
    )
    recovered1 = reconstruct_region(result, frequency_mask(result, frequency1, 6.0))
    recovered2 = reconstruct_region(result, frequency_mask(result, frequency2, 6.0))
    interior = 129:896
    @test norm(recovered1[interior] - component1[interior]) /
          norm(component1[interior]) < 0.01
    @test norm(recovered2[interior] - component2[interior]) /
          norm(component2[interior]) < 0.01

    approaching_frequency1 = 50.0 .+ 15.0 .* time
    approaching_frequency2 = 100.0 .- 10.0 .* time
    approaching1 = cos.(2π .* (50.0 .* time .+ 7.5 .* time .^ 2))
    approaching2 = 0.7 .* cos.(2π .* (100.0 .* time .- 5.0 .* time .^ 2))
    approaching = wsst(
        approaching1 + approaching2;
        fs,
        freqrange = (30.0, 120.0),
        nfrequencies = 361,
        preserve_residual = false,
    )
    estimate1 = reconstruct_region(
        approaching,
        frequency_mask(approaching, approaching_frequency1, 5.0),
    )
    estimate2 = reconstruct_region(
        approaching,
        frequency_mask(approaching, approaching_frequency2, 5.0),
    )
    error1 = norm(estimate1[interior] - approaching1[interior]) / norm(approaching1[interior])
    error2 = norm(estimate2[interior] - approaching2[interior]) / norm(approaching2[interior])
    @test all(isfinite, estimate1)
    @test all(isfinite, estimate2)
    @test max(error1, error2) > 0.1 # Deliberate separation-limit case.
end
