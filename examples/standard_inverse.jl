using LinearAlgebra: norm
using SynchroSqueezing

fs = 1024.0
n = 1024
time = (0:n-1) ./ fs
oscillation = cos.(2π * 64.0 .* time)
signal = 0.2 .+ oscillation

result = wsst(
    signal;
    fs,
    freqrange = (40.0, 80.0),
    nfrequencies = 161,
    threshold = NoThreshold(),
    preserve_residual = true,
)

oscillatory = inverse_wsst(result)
full = inverse_wsst(result; mode = :full)

println("standard oscillatory relative error: ", norm(oscillatory - oscillation) / norm(oscillation))
println("full-input relative error: ", norm(full - signal) / norm(signal))
println("admissibility constant: ", real(result.reconstruction.admissibility))
println("1% reliable frequency range: ", result.reconstruction.reliable_frequency_range, " Hz")
