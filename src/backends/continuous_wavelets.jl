const _CW = ContinuousWavelets

"""Return the real scalar type associated with a real or complex element type."""
_realtype(::Type{T}) where {T<:Real} = T
_realtype(::Type{Complex{T}}) where {T<:Real} = T

"""Map the currently supported boundary choice to the backend configuration."""
function _continuous_wavelets_boundary(boundary::Symbol)
    boundary === :periodic && return _CW.Periodic
    throw(ArgumentError("only boundary=:periodic is supported by the spectral derivative adapter; got $boundary"))
end

"""Return physical angular frequencies for an `rfft` spectrum."""
function _rfft_angular_frequencies(n::Int, fs::T) where {T<:AbstractFloat}
    bins = T.(0:(n ÷ 2))
    return (T(2π) * fs / T(n)) .* bins
end

"""Validate each wavelet row using the relative DC response of its Fourier daughter."""
function _validate_wavelets(daughters::AbstractMatrix, tolerance::T) where {T<:AbstractFloat}
    peak = vec(maximum(abs, daughters; dims = 1))
    leakage = T.(abs.(daughters[1, :]) ./ peak)
    return WaveletValidation(leakage, BitVector(leakage .< tolerance), tolerance)
end

"""Return midpoint-bin widths for a strictly ascending physical-frequency grid."""
function _frequency_weights(frequencies::Vector{T}) where {T<:AbstractFloat}
    length(frequencies) > 1 || throw(ArgumentError("at least two CWT frequencies are required"))
    all(isfinite, frequencies) && all(>(zero(T)), frequencies) && issorted(frequencies) ||
        throw(ArgumentError("backend returned invalid physical frequencies"))

    steps = diff(frequencies)
    all(>(zero(T)), steps) || throw(ArgumentError("backend frequencies must be strictly ascending"))
    weights = similar(frequencies)
    weights[1] = steps[1]
    weights[end] = steps[end]
    @views weights[2:end-1] .= (steps[1:end-1] .+ steps[2:end]) ./ 2
    return weights
end

"""
    plan_cwt(n; fs, wavelet=ContinuousWavelets.Morlet(π), voices_per_octave=8,
             beta=2, precision=Float64, dc_tolerance=0.01)

Build a periodic `ContinuousWavelets.jl` plan. A periodic boundary is required
in this milestone because it makes the FFT-based time derivative exact under
the recorded transform convention.
"""
function plan_cwt(
    n::Integer;
    fs::Real,
    wavelet = _CW.Morlet(π),
    voices_per_octave::Integer = 8,
    beta::Real = 2,
    precision::Type{T} = Float64,
    dc_tolerance::Real = 0.01,
    boundary::Symbol = :periodic,
    inverse_style = _CW.PenroseDelta(),
) where {T<:AbstractFloat}
    n > 1 || throw(ArgumentError("signal length must exceed one"))
    fs > 0 || throw(ArgumentError("fs must be positive"))
    voices_per_octave > 0 || throw(ArgumentError("voices_per_octave must be positive"))
    0 < dc_tolerance < 1 || throw(ArgumentError("dc_tolerance must lie in (0, 1)"))

    backend_config = _CW.wavelet(
        wavelet;
        Q = voices_per_octave,
        β = beta,
        fsample = T(fs),
        averagingType = _CW.NoAve(),
        boundary = _continuous_wavelets_boundary(boundary),
    )
    daughters, _ = _CW.computeWavelets(n, backend_config; T = T)
    _, _, scales, _ = _CW.getNWavelets(n, backend_config)
    frequencies = T.(_CW.getMeanFreq(daughters, T(fs)))
    frequency_weights = _frequency_weights(frequencies)
    validation = _validate_wavelets(daughters, T(dc_tolerance))

    return CWTPlan(
        ContinuousWaveletsBackend(),
        backend_config,
        daughters,
        T.(scales),
        frequencies,
        frequency_weights,
        inverse_style,
        validation,
        T(fs),
        Int(n),
    )
end

"""Extract the positive-frequency signal spectrum in the backend's convention."""
function _positive_spectrum(signal::AbstractVector{T}) where {T<:Real}
    return rfft(signal)
end

function _positive_spectrum(signal::AbstractVector{T}) where {T<:Complex}
    return @view fft(signal)[1:(length(signal) ÷ 2 + 1)]
end

"""
Compute `∂ₜW` from the same Fourier daughters used by `ContinuousWavelets.jl`.

For the periodic CWT, the backend applies a one-sided daughter bank to the
positive signal spectrum and takes a complex inverse FFT. Multiplication by
`im * ω` before that inverse FFT is therefore the exact spectral time
derivative under the package's recorded Fourier convention.
"""
function _spectral_derivative(
    signal::AbstractVector,
    plan::CWTPlan{ContinuousWaveletsBackend},
)
    n = length(signal)
    x̂ = _positive_spectrum(signal)
    nfrequency = length(x̂)
    nscales = size(plan.daughters, 2)
    output_type = Complex{_realtype(eltype(signal))}
    spectrum = zeros(output_type, n, nscales)
    angular_frequency = _rfft_angular_frequencies(n, plan.sample_rate)

    @views spectrum[1:nfrequency, :] .=
        reshape(x̂, :, 1) .* plan.daughters .* reshape(im .* angular_frequency, :, 1)
    return ifft(spectrum, 1)
end

"""Create result metadata whose conventions are sufficient for inverse CWT."""
function _metadata(plan::CWTPlan{ContinuousWaveletsBackend})
    return CWTMetadata(
        plan.sample_rate,
        inv(plan.sample_rate),
        plan.signal_length,
        plan.signal_length,
        1:plan.signal_length,
        :periodic,
        :continuouswavelets_linf,
        plan.inverse_style,
        plan.backend_config,
    )
end

"""
    cwt(plan::CWTPlan{ContinuousWaveletsBackend}, signal)

Compute the selected backend CWT and an FFT-exact physical-time derivative.
The returned matrices are scale/frequency × time non-copying views.
"""
function cwt(plan::CWTPlan{ContinuousWaveletsBackend}, signal::AbstractVector)
    length(signal) == plan.signal_length || throw(DimensionMismatch(
        "plan length $(plan.signal_length) does not match signal length $(length(signal))",
    ))
    _realtype(eltype(signal)) === typeof(plan.sample_rate) || throw(ArgumentError(
        "signal precision $(_realtype(eltype(signal))) does not match plan precision $(typeof(plan.sample_rate))",
    ))

    backend_coefficients = _CW.cwt(signal, plan.backend_config, plan.daughters)
    backend_derivative = _spectral_derivative(signal, plan)
    coefficients = PermutedDimsArray(backend_coefficients, (2, 1))
    derivative = PermutedDimsArray(backend_derivative, (2, 1))
    reliability = repeat(plan.validation.valid_rows, 1, plan.signal_length)

    return CWTResult(
        coefficients,
        derivative,
        copy(plan.frequencies),
        copy(plan.scales),
        copy(plan.frequency_weights),
        copy(plan.validation.valid_rows),
        reliability,
        _metadata(plan),
    )
end

"""Convenience method that plans and computes a CWT for a single signal."""
function cwt(signal::AbstractVector; fs::Real, kwargs...)
    precision = _realtype(eltype(signal))
    plan = plan_cwt(length(signal); fs, precision, kwargs...)
    return cwt(plan, signal)
end

"""Recover the input signal with the inverse convention stored in `result`."""
function inverse_cwt(result::CWTResult)
    backend_coefficients = parent(result.coefficients)
    recovered = _CW.icwt(
        backend_coefficients,
        result.metadata.backend_config,
        result.metadata.inverse_style,
    )
    return vec(recovered)
end

"""
    instantaneous_frequency(result; threshold=0)

Return the phase-based physical frequency estimate in Hz. Values below the
coefficient threshold or on analytically unsuitable rows are marked invalid and
stored as `NaN`; they are never replaced by zero.
"""
function instantaneous_frequency(result::CWTResult; threshold::Real = 0)
    threshold >= 0 || throw(ArgumentError("threshold must be nonnegative"))
    coefficients = result.coefficients
    derivative = result.derivative
    values = fill(NaN, size(coefficients))
    valid = falses(size(coefficients))
    threshold2 = threshold^2

    for time_index in axes(coefficients, 2), row in axes(coefficients, 1)
        coefficient = coefficients[row, time_index]
        magnitude2 = abs2(coefficient)
        if result.reliability[row, time_index] && isfinite(magnitude2) && magnitude2 > threshold2
            values[row, time_index] = imag(derivative[row, time_index] * conj(coefficient)) /
                                      (2π * magnitude2)
            valid[row, time_index] = isfinite(values[row, time_index])
        end
    end
    return FrequencyEstimate(values, valid)
end
