#!/usr/bin/env julia

using SynchroSqueezing

"""Run the Milestone 02 spectral-derivative validation on a periodic pure tone."""
function run_derivative_experiment(; n::Int = 1024, fs::Float64 = 1024.0, frequency::Float64 = 64.0)
    time = (0:(n - 1)) ./ fs
    signal = cos.(2π * frequency .* time)
    plan = plan_cwt(n; fs, boundary = :periodic)
    result = cwt(plan, signal)

    centered_difference = (result.coefficients[:, 3:end] .-
                          result.coefficients[:, 1:end-2]) .* (fs / 2)
    spectral_derivative = result.derivative[:, 2:end-1]
    difference_error = sqrt(sum(abs2, spectral_derivative - centered_difference) /
                            sum(abs2, spectral_derivative))

    estimate = instantaneous_frequency(result; threshold = 1e-8)
    peak = argmax(abs2.(result.coefficients))
    row, sample = Tuple(CartesianIndices(result.coefficients)[peak])

    println("CWT derivative experiment")
    println("  layout (frequency × time): ", size(result.coefficients))
    println("  peak frequency estimate (Hz): ", estimate.values[row, sample])
    println("  centered-difference relative error: ", difference_error)
    println("  reliable wavelet rows: ", count(result.valid_rows), "/", length(result.valid_rows))
    println("  inverse relative error: ", sqrt(sum(abs2, inverse_cwt(result) - signal) / sum(abs2, signal)))

    return (; result, difference_error, peak_frequency = estimate.values[row, sample])
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_derivative_experiment()
end
