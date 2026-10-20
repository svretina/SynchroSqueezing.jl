using LinearAlgebra
using SynchroSqueezing

fs = 1024.0
n = 1024
time = (0:n-1) ./ fs
signal = 0.25 .+ cos.(2π * 64.0 .* time)

result = wsst(
    signal;
    fs,
    freqrange = (50.0, 75.0),
    nfrequencies = 101,
    preserve_residual = true,
)
oscillatory = inverse_wsst(result; mode = :oscillatory)
residual = inverse_wsst(result; mode = :residual)
reconstructed = inverse_wsst(result; mode = :full)

println("oscillatory norm: ", norm(oscillatory))
println("residual/DC mean: ", sum(residual) / length(residual))
println("full relative error: ", norm(reconstructed - signal) / norm(signal))
