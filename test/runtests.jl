using Test
using SynchroSqueezing

include("signals.jl")

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
