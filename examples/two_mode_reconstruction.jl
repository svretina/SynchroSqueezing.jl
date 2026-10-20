using LinearAlgebra: norm
using SynchroSqueezing

fs = 1024.0
n = 1024
time = (0:n-1) ./ fs
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

function band_mask(result, frequency, half_width)
    return repeat(abs.(result.frequencies .- frequency) .<= half_width, 1, n)
end

recovered1 = reconstruct_region(result, band_mask(result, 40.0, 4.0))
recovered2 = reconstruct_region(result, band_mask(result, 96.0, 4.0))

println("40 Hz component relative error: ", norm(recovered1 - component1) / norm(component1))
println("96 Hz component relative error: ", norm(recovered2 - component2) / norm(component2))
println("combined relative error: ", norm(recovered1 + recovered2 - signal) / norm(signal))
