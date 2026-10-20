@testset "First-order WSST reassignment" begin
    fs = 1024.0
    n = 1024
    frequency = 64.0
    _, tone = constant_tone(n, fs, frequency)
    cwt_result = cwt(plan_cwt(n; fs), tone)

    nearest = wsst(
        cwt_result;
        freqrange = (40.0, 80.0),
        nfrequencies = 81,
        threshold = AbsoluteThreshold(1e-8),
        keep_cwt = true,
        keep_derivative = true,
        keep_phase = true,
    )
    @test size(nearest.coefficients) == (81, n)
    @test issorted(nearest.frequencies)
    @test length(nearest.frequency_edges) == 82
    @test nearest.frequency_spacing === :linear
    @test nearest.cwt === cwt_result.coefficients
    @test nearest.derivative === cwt_result.derivative
    @test nearest.phase !== nothing
    @test all(nearest.diagnostics.accepted .> 0)
    @test nearest.diagnostics.accepted_sum ≈ nearest.diagnostics.deposited_sum atol = 1e-11

    row = argmax(vec(sum(abs2, nearest.coefficients; dims = 2)))
    @test abs(nearest.frequencies[row] - frequency) <= 0.5

    linear = wsst(
        cwt_result;
        freqrange = (40.0, 80.0),
        nfrequencies = 81,
        reassignment = :linear,
        threshold = RelativeThreshold(1e-8),
    )
    @test linear.cwt === nothing
    @test linear.diagnostics.accepted_sum ≈ linear.diagnostics.deposited_sum atol = 1e-11
    @test sum(abs2, linear.coefficients) <= sum(abs2, nearest.coefficients) + 1e-8

    reference = wsst(
        cwt_result;
        freqrange = (40.0, 80.0),
        nfrequencies = 81,
        reassignment = :linear,
        threshold = RelativeThreshold(1e-8),
        implementation = :reference,
    )
    @test reference.coefficients ≈ linear.coefficients atol = 1e-12
    @test reference.residual ≈ linear.residual atol = 1e-12
    @test reference.diagnostics.accepted == linear.diagnostics.accepted

    logarithmic = wsst(
        cwt_result;
        freqrange = (32.0, 128.0),
        nfrequencies = 17,
        frequency_spacing = :log,
    )
    @test logarithmic.frequency_spacing === :log
    @test all(>(0), logarithmic.frequency_edges)

    explicit = wsst(
        cwt_result;
        frequency_grid = [40.0, 52.0, 64.0, 80.0],
        frequency_edges = [34.0, 46.0, 58.0, 72.0, 88.0],
    )
    @test explicit.frequency_spacing === :arbitrary
    @test explicit.frequency_edges == [34.0, 46.0, 58.0, 72.0, 88.0]
end

@testset "Randomized reference and fused kernel agreement" begin
    rng = MersenneTwister(0x5a17)
    fs = 256.0
    n = 129
    base = cwt(plan_cwt(n; fs), randn(rng, n))
    coefficients = randn(rng, ComplexF64, size(base.coefficients))
    target_frequencies = 5 .+ 100 .* rand(rng, size(coefficients))
    derivative = im .* (2π .* target_frequencies) .* coefficients
    synthetic = CWTResult(
        coefficients,
        derivative,
        base.frequencies,
        base.scales,
        base.frequency_weights,
        base.reconstruction_weights,
        base.squeezing_weights,
        base.synthesis_response,
        base.admissibility,
        base.standard_compatible,
        base.reliable_frequency_range,
        trues(size(coefficients, 1)),
        trues(size(coefficients)),
        base.metadata,
        base.dc_component,
    )
    optimized = wsst(
        synthetic;
        freqrange = (5.0, 105.0),
        nfrequencies = 101,
        reassignment = :linear,
        threshold = AbsoluteThreshold(0.2),
    )
    reference = wsst(
        synthetic;
        freqrange = (5.0, 105.0),
        nfrequencies = 101,
        reassignment = :linear,
        threshold = AbsoluteThreshold(0.2),
        implementation = :reference,
    )
    @test optimized.coefficients == reference.coefficients
    @test optimized.residual == reference.residual
    @test optimized.diagnostics.accepted == reference.diagnostics.accepted
    @test optimized.diagnostics.rejected_out_of_range ==
          reference.diagnostics.rejected_out_of_range
end

@testset "WSST precision, zero input, and validation" begin
    fs = Float32(512)
    n = 513
    _, tone = constant_tone(n, fs, Float32(32))
    cwt_result = cwt(plan_cwt(n; fs, precision = Float32), Float32.(tone))
    squeezed = wsst(cwt_result; freqrange = (Float32(20), Float32(50)), nfrequencies = 31)
    @test eltype(squeezed.coefficients) == ComplexF32
    @test eltype(squeezed.frequencies) == Float32

    zero_result = cwt(plan_cwt(n; fs, precision = Float32), zeros(Float32, n))
    zero_sst = wsst(zero_result; freqrange = (Float32(20), Float32(50)), nfrequencies = 31)
    @test all(iszero, zero_sst.coefficients)
    @test all(==(0), zero_sst.diagnostics.accepted)
    @test all(>(0), zero_sst.diagnostics.rejected_magnitude)

    @test_throws ArgumentError wsst(cwt_result; frequency_grid = [20.0, 30.0], freqrange = (20.0, 30.0))
    @test_throws ArgumentError wsst(cwt_result; reassignment = :unsupported)
    @test_throws ArgumentError wsst(cwt_result; frequency_grid = [20.0, 20.0])
    @test_throws ArgumentError wsst(cwt_result; freqrange = (20.0, 256.0))
end
