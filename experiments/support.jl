module BackendExperimentSupport

using LinearAlgebra
using Statistics

export analytical_signals, relative_error, peak_frequency, warm_seconds

"""Return deterministic tone and chirp fixtures on a no-duplicate-endpoint grid."""
function analytical_signals(n::Integer, fs::Real)
    t = collect(0:(n - 1)) ./ float(fs)
    tone = cos.(2π .* 64.0 .* t)
    chirp = cos.(2π .* (32.0 .* t .+ 20.0 .* t .^ 2))
    return (; t, tone, chirp)
end

relative_error(reference, recovered) = norm(vec(recovered) - vec(reference)) / norm(vec(reference))

"""Return the frequency with the largest scale-averaged coefficient energy."""
function peak_frequency(coefficients, frequencies)
    scores = vec(mean(abs2, coefficients; dims = 1))
    return frequencies[argmax(scores)]
end

"""Time one already-compiled operation without including compilation time."""
function warm_seconds(f)
    f()
    return @elapsed f()
end

end # module BackendExperimentSupport
