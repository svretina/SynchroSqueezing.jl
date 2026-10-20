@testset "Standalone instantaneous-frequency estimator" begin
    frequency = 37.0
    coefficient = ComplexF64[1 + 2im 0 + 0im 1e-300 + 1e-300im 1e150 + 1e150im]
    derivative = im * 2π * frequency .* coefficient
    estimate = instantaneous_frequency(coefficient, derivative; threshold = NoThreshold())

    @test estimate.valid[1]
    @test estimate.values[1] ≈ frequency atol = 1e-12
    @test !estimate.valid[2]
    @test isnan(estimate.values[2])
    @test estimate.valid[3]
    @test estimate.values[3] ≈ frequency atol = 1e-12
    @test estimate.valid[4]
    @test estimate.values[4] ≈ frequency atol = 1e-12

    zero_frequency = instantaneous_frequency(
        reshape(ComplexF64[1], 1, 1),
        reshape(ComplexF64[0], 1, 1),
    )
    @test zero_frequency.valid[1]
    @test zero_frequency.values[1] == 0

    negative = instantaneous_frequency(
        reshape(ComplexF64[1], 1, 1),
        reshape(ComplexF64[-im * 2π * frequency], 1, 1),
    )
    @test negative.valid[1]
    @test negative.values[1] ≈ -frequency

    bad_coefficients = reshape(ComplexF64[Inf, NaN + 0im, 1], 1, 3)
    bad_derivatives = reshape(ComplexF64[1, 1, Inf], 1, 3)
    bad = instantaneous_frequency(bad_coefficients, bad_derivatives)
    @test !any(bad.valid)
    @test all(isnan, bad.values)

    relative = instantaneous_frequency(
        reshape(ComplexF64[1, 1e-6], 1, 2),
        reshape(ComplexF64[im, im * 1e-6], 1, 2);
        threshold = RelativeThreshold(1e-3),
    )
    @test relative.valid == BitMatrix([true false])
    @test_throws ArgumentError instantaneous_frequency(
        reshape(ComplexF64[1], 1, 1),
        reshape(ComplexF64[im], 1, 1);
        threshold = RelativeThreshold(1e-3; reference = :scale),
    )
    @test_throws ArgumentError AbsoluteThreshold(-1.0)
    @test_throws ArgumentError AbsoluteThreshold(NaN)
    @test_throws ArgumentError RelativeThreshold(-1.0)
    @test_throws DimensionMismatch instantaneous_frequency(
        zeros(ComplexF64, 2, 2),
        zeros(ComplexF64, 2, 3),
    )
end

@testset "Physical-time phase units" begin
    for fs in (512.0, 1024.0)
        n = Int(fs)
        t = (0:n-1) ./ fs
        signal = cis.(2π * 32.0 .* t)
        result = cwt(plan_cwt(n; fs), signal)
        estimate = instantaneous_frequency(result; threshold = AbsoluteThreshold(1e-8))
        peak = argmax(abs2.(result.coefficients))
        @test estimate.values[peak] ≈ 32.0 atol = 1e-9
    end
end
