"""Uniform sample times without a duplicated endpoint."""
sample_times(n::Integer, fs::Real) = collect(0:(n - 1)) ./ float(fs)

"""A deterministic real sinusoid sampled at `fs` Hz."""
function constant_tone(n::Integer, fs::Real, frequency::Real; phase::Real = 0.0)
    t = sample_times(n, fs)
    return t, cos.(2π .* frequency .* t .+ phase)
end

"""A deterministic linear-frequency chirp with `f(t) = f0 + rate * t`."""
function linear_chirp(n::Integer, fs::Real, f0::Real, rate::Real; phase::Real = 0.0)
    t = sample_times(n, fs)
    return t, cos.(2π .* (f0 .* t .+ 0.5 .* rate .* t .^ 2) .+ phase)
end
