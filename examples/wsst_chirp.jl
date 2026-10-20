using SynchroSqueezing

fs = 512.0
n = 512
time = (0:n-1) ./ fs
exact_frequency = 30.0 .+ 40.0 .* time
phase = 2π .* (30.0 .* time .+ 20.0 .* time .^ 2)
signal = cos.(phase)

result = wsst(
    signal;
    fs,
    freqrange = (20.0, 80.0),
    nfrequencies = 241,
)
ridge = [result.frequencies[argmax(abs2.(@view result.coefficients[:, column]))]
         for column in axes(result.coefficients, 2)]
interior = 65:448
rmse = sqrt(sum(abs2, ridge[interior] - exact_frequency[interior]) / length(interior))

println("interior instantaneous-frequency RMSE: ", rmse, " Hz")
