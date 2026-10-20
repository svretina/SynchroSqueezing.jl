using SynchroSqueezing

fs = 1024.0
n = 1024
time = (0:n-1) ./ fs
signal = cos.(2π * 64.0 .* time)

result = wsst(
    signal;
    fs,
    freqrange = (40.0, 80.0),
    nfrequencies = 81,
    threshold = RelativeThreshold(1e-6),
)
energy = vec(sum(abs2, result.coefficients; dims = 2))
peak_frequency = result.frequencies[argmax(energy)]

println("dominant WSST frequency: ", peak_frequency, " Hz")
println("coefficient shape: ", size(result.coefficients))
