using LinearAlgebra

@testset "ContinuousWavelets backend adapter" begin
    fs = 1024.0
    n = 1024
    frequency = 64.0
    _, tone = constant_tone(n, fs, frequency)
    plan = plan_cwt(n; fs, voices_per_octave = 8, boundary = :periodic)
    result = cwt(plan, tone)

    @test size(result.coefficients) == (length(plan.frequencies), n)
    @test size(result.derivative) == size(result.coefficients)
    @test issorted(result.frequencies)
    @test all(>(0), result.frequency_weights)
    @test result.metadata.boundary === :periodic
    @test result.metadata.sample_rate == fs
    @test result.metadata.padded_length == n
    @test result.metadata.crop == 1:n
    @test parent(result.coefficients) isa Matrix
    @test parent(result.derivative) isa Matrix
    @test size(result.reliability) == size(result.coefficients)

    centered_difference = (result.coefficients[:, 3:end] .-
                          result.coefficients[:, 1:end-2]) .* (fs / 2)
    interior_derivative = result.derivative[:, 2:end-1]
    @test norm(interior_derivative - centered_difference) / norm(interior_derivative) < 0.04

    recovered = inverse_cwt(result)
    @test length(recovered) == n
    @test norm(recovered - tone) / norm(tone) < 0.05

    estimate = instantaneous_frequency(result; threshold = 1e-8)
    energy = abs2.(result.coefficients)
    row, time = Tuple(CartesianIndices(energy)[argmax(energy)])
    @test estimate.valid[row, time]
    @test abs(estimate.values[row, time] - frequency) < 1e-8
    @test all(isnan, instantaneous_frequency(result; threshold = Inf).values)
end

@testset "CWT plan precision, complex signals, and constraints" begin
    fs32 = Float32(512)
    n = 512
    _, tone = constant_tone(n, fs32, Float32(32))
    plan32 = plan_cwt(n; fs = fs32, precision = Float32)
    result32 = cwt(plan32, Float32.(tone))
    @test eltype(result32.coefficients) == ComplexF32

    analytic_tone = complex.(Float64.(tone))
    plan64 = plan_cwt(n; fs = 512.0)
    complex_result = cwt(plan64, analytic_tone)
    @test eltype(complex_result.coefficients) == ComplexF64

    @test_throws ArgumentError plan_cwt(n; fs = 512.0, boundary = :reflect)
    @test_throws ArgumentError cwt(plan64, Float32.(tone))
    @test_throws DimensionMismatch cwt(plan64, tone[1:end-1])
end
