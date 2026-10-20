@testset "Inverse WSST and residual accounting" begin
    fs = 1024.0
    n = 1024
    t = sample_times(n, fs)
    signal = 0.35 .+ cos.(2π * 64.0 .* t)
    cwt_result = cwt(plan_cwt(n; fs, dc_tolerance = 0.99), signal)

    nearest = wsst(
        cwt_result;
        freqrange = (40.0, 80.0),
        nfrequencies = 81,
        threshold = NoThreshold(),
        reassignment = :nearest,
    )
    linear = wsst(
        cwt_result;
        freqrange = (40.0, 80.0),
        nfrequencies = 81,
        threshold = NoThreshold(),
        reassignment = :linear,
    )

    nearest_full = inverse_wsst(nearest; mode = :full)
    linear_full = inverse_wsst(linear; mode = :full)
    backend_full = inverse_cwt(cwt_result)
    @test norm(nearest_full - signal) / norm(signal) < 1e-3
    @test norm(linear_full - signal) / norm(signal) < 1e-3
    @test maximum(abs, nearest_full - linear_full) < 1e-10
    @test norm(backend_full - signal) / norm(signal) < 1e-3
    @test maximum(abs,
        inverse_wsst(nearest; mode = :oscillatory) .+
        inverse_wsst(nearest; mode = :residual) .- nearest_full,
    ) < 1e-10

    high_threshold = wsst(
        cwt_result;
        freqrange = (40.0, 80.0),
        nfrequencies = 81,
        threshold = AbsoluteThreshold(Inf),
    )
    @test all(iszero, high_threshold.coefficients)
    @test norm(inverse_wsst(high_threshold; mode = :full) - signal) / norm(signal) < 1e-3

    no_residual = wsst(cwt_result; preserve_residual = false)
    @test_throws ArgumentError inverse_wsst(no_residual; mode = :full)
    @test_throws ArgumentError inverse_wsst(no_residual; mode = :residual)
    @test_throws ArgumentError inverse_wsst(nearest; mode = :unknown)
end

@testset "Float32 and Nyquist synthesis conventions" begin
    fs32 = Float32(512)
    n = 512
    time32 = Float32.(0:n-1) ./ fs32
    signal32 = Float32.(0.1 .+ cos.(Float32(2π * 32) .* time32))
    result32 = wsst(
        signal32;
        fs = fs32,
        freqrange = (Float32(20), Float32(50)),
        nfrequencies = 61,
    )
    reconstructed32 = inverse_wsst(result32; mode = :full)
    @test eltype(reconstructed32) == Float32
    @test norm(reconstructed32 - signal32) / norm(signal32) < 2e-3

    nyquist = Float64.((-1) .^ (0:n-1))
    nyquist_result = wsst(nyquist; fs = 512.0, threshold = AbsoluteThreshold(Inf))
    standard_nyquist = inverse_wsst(nyquist_result; mode = :full)
    @test norm(standard_nyquist - nyquist) / norm(nyquist) > 0.1
    @test inverse_wsst(nyquist_result; mode = :full, method = :frame) ≈ nyquist atol = 1e-12
end

@testset "Unrecoverable coefficient reporting" begin
    result = cwt(plan_cwt(128; fs = 128.0), zeros(128))
    coefficients = Matrix(result.coefficients)
    coefficients[1, 1] = ComplexF64(Inf, 0)
    damaged = CWTResult(
        coefficients,
        Matrix(result.derivative),
        result.frequencies,
        result.scales,
        result.frequency_weights,
        result.reconstruction_weights,
        result.squeezing_weights,
        result.synthesis_response,
        result.admissibility,
        result.standard_compatible,
        result.reliable_frequency_range,
        result.valid_rows,
        result.reliability,
        result.metadata,
        result.dc_component,
    )
    squeezed = wsst(damaged)
    @test squeezed.diagnostics.unrecoverable[1] == 1
    @test_throws ArgumentError inverse_wsst(squeezed; mode = :full)
end

@testset "Real DC and analytic-complex reconstruction" begin
    fs = 512.0
    n = 512
    constant_signal = fill(2.5, n)
    constant_result = wsst(constant_signal; fs)
    @test all(iszero, constant_result.coefficients)
    @test inverse_wsst(constant_result; mode = :full) ≈ constant_signal atol = 1e-12

    t = sample_times(n, fs)
    analytic_signal = 0.2 + 0.1im .+ cis.(2π * 32.0 .* t)
    analytic_cwt = cwt(plan_cwt(n; fs), analytic_signal)
    @test analytic_cwt.metadata.signal_kind === :analytic_complex
    analytic_result = wsst(
        analytic_cwt;
        freqrange = (20.0, 50.0),
        nfrequencies = 61,
        threshold = NoThreshold(),
    )
    reconstructed = inverse_wsst(analytic_result; mode = :full)
    @test eltype(reconstructed) == ComplexF64
    @test norm(reconstructed - analytic_signal) / norm(analytic_signal) < 1e-3
    @test norm(inverse_cwt(analytic_cwt) - analytic_signal) / norm(analytic_signal) < 1e-3

    negative_signal = cis.(-2π * 32.0 .* t)
    negative_cwt = cwt(plan_cwt(n; fs), negative_signal)
    @test negative_cwt.metadata.signal_kind === :two_sided_complex
    @test_throws ArgumentError wsst(negative_cwt)
    @test_throws ArgumentError inverse_cwt(negative_cwt)
end
