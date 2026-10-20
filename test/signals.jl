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

"""A chirp with instantaneous frequency `f0 + rate*t + curvature*t^2`."""
function quadratic_chirp(n::Integer, fs::Real, f0::Real, rate::Real, curvature::Real)
    t = sample_times(n, fs)
    phase = 2π .* (f0 .* t .+ 0.5 .* rate .* t .^ 2 .+ (curvature / 3) .* t .^ 3)
    return t, cos.(phase), f0 .+ rate .* t .+ curvature .* t .^ 2
end

"""An exponential-frequency chirp with `f(t) = f0 * exp(rate*t)`."""
function exponential_chirp(n::Integer, fs::Real, f0::Real, rate::Real)
    t = sample_times(n, fs)
    frequency = f0 .* exp.(rate .* t)
    phase = 2π .* f0 .* expm1.(rate .* t) ./ rate
    return t, cos.(phase), frequency
end

"""A sinusoid with a sinusoidal amplitude envelope."""
function amplitude_modulated_tone(n::Integer, fs::Real, frequency::Real; depth = 0.25, modulation = 3.0)
    t = sample_times(n, fs)
    envelope = 1 .+ depth .* cos.(2π .* modulation .* t)
    return t, envelope .* cos.(2π .* frequency .* t), envelope
end
