@testset "WSST frequency grids" begin
    linear = linear_frequency_grid(10.0, 30.0, 5)
    @test linear == [10.0, 15.0, 20.0, 25.0, 30.0]
    logarithmic = log_frequency_grid(10.0, 160.0, 5)
    @test logarithmic ≈ [10.0, 20.0, 40.0, 80.0, 160.0]
    @test issorted(logarithmic)

    @test_throws ArgumentError linear_frequency_grid(0.0, 30.0, 5)
    @test_throws ArgumentError linear_frequency_grid(10.0, 10.0, 5)
    @test_throws ArgumentError log_frequency_grid(-1.0, 30.0, 5)

    centers = [10.0, 20.0, 40.0]
    edges = [5.0, 15.0, 30.0, 50.0]
    output = zeros(ComplexF64, 3, 1)
    @test SynchroSqueezing._deposit_linear!(output, 30.0, 2 + 3im, centers, edges, 1)
    @test sum(output) == 2 + 3im
    @test sum(centers .* real.(output[:, 1])) / real(sum(output)) == 30.0
    fill!(output, 0)
    @test !SynchroSqueezing._deposit_linear!(output, 4.0, 1 + 0im, centers, edges, 1)
    @test all(iszero, output)
end
