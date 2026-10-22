using Test
using SynchroSqueezing
using Random

include("signals.jl")
include("test_backends.jl")
include("test_phase.jl")
include("test_grids.jl")
include("test_wsst.jl")
include("test_inverse_wsst.jl")
include("test_standard_reconstruction.jl")
include("test_wsst_validation.jl")
include("test_ridges.jl")
include("test_predictive_ridges.jl")
include("test_ridge_reconstruction.jl")
include("test_damping.jl")
include("test_mode_support.jl")
include("test_scientific_regression.jl")

@testset "SynchroSqueezing.jl package foundation" begin
    @test SynchroSqueezing.package_version() == v"0.1.0-DEV"
end

@testset "deterministic analytical signal helpers" begin
    fs = 128.0
    t, x = constant_tone(128, fs, 8.0)
    @test length(t) == 128
    @test length(x) == 128
    @test t[1] == 0.0
    @test t[end] == (length(t) - 1) / fs

    _, chirp = linear_chirp(128, fs, 8.0, 2.0)
    @test eltype(chirp) == Float64
    @test chirp == linear_chirp(128, fs, 8.0, 2.0)[2]
end
