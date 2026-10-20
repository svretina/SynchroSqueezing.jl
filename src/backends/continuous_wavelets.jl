const _CW = ContinuousWavelets
const _DEFAULT_STANDARD_WAVELET = _CW.Morse(3.0, 20.0, 1.0)

"""Return the real scalar type associated with a real or complex element type."""
_realtype(::Type{T}) where {T<:Real} = T
_realtype(::Type{Complex{T}}) where {T<:Real} = T

"""Classify whether the one-sided analytic transform can represent the input."""
_signal_kind(signal::AbstractVector{<:Real}) = :real

function _signal_kind(signal::AbstractVector{Complex{T}}) where {T<:AbstractFloat}
    spectrum = fft(signal)
    first_negative = length(signal) ÷ 2 + 2
    scale = maximum(abs, spectrum)
    isfinite(scale) || throw(ArgumentError("signal samples must be finite"))
    iszero(scale) && return :analytic_complex
    total_energy = sum(abs2, spectrum ./ scale)
    negative_energy = sum(abs2, @view(spectrum[first_negative:end]) ./ scale)
    tolerance = sqrt(eps(T)) * max(total_energy, floatmin(T))
    return negative_energy <= tolerance ? :analytic_complex : :two_sided_complex
end

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

"""Obtain the documented Penrose-delta synthesis weights for this wavelet bank."""
function _reconstruction_weights(daughters, backend_config, inverse_style, ::Type{T}) where {T<:AbstractFloat}
    inverse_style isa _CW.PenroseDelta || return nothing
    weights = vec(_CW.computeDualWeights(daughters, backend_config))
    maximum(abs, imag.(weights)) <= sqrt(eps(T)) * max(maximum(abs, weights), one(T)) ||
        throw(ArgumentError("the backend returned non-real Penrose reconstruction weights"))
    return T.(real.(weights))
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

"""Trapezoidal quadrature weights in the backend's log-frequency coordinate."""
function _log_frequency_weights(frequencies::Vector{T}) where {T<:AbstractFloat}
    log_frequencies = log.(frequencies)
    steps = diff(log_frequencies)
    weights = similar(frequencies)
    weights[1] = steps[1] / 2
    weights[end] = steps[end] / 2
    @views weights[2:end-1] .= (steps[1:end-1] .+ steps[2:end]) ./ 2
    return weights
end

"""Unnormalized analytic Morse response used by ContinuousWavelets.jl."""
function _morse_response(wavelet::_CW.Morse, frequency::Real)
    frequency >= 0 || return zero(float(frequency))
    wavelet.ga > 0 || throw(ArgumentError(
        "standard reconstruction requires a Morse wavelet with gamma > 0",
    ))
    wavelet.be > 0 || throw(ArgumentError(
        "standard reconstruction requires a Morse wavelet with beta > 0",
    ))
    wavelet.cf > 0 || throw(ArgumentError(
        "standard reconstruction requires a positive Morse center frequency",
    ))
    iszero(frequency) && return zero(float(frequency))
    gamma_parameter = wavelet.ga
    beta_parameter = wavelet.be
    peak_frequency = (beta_parameter / gamma_parameter)^(1 / gamma_parameter)
    frequency_factor = wavelet.cf / peak_frequency
    normalized_frequency = frequency / frequency_factor
    amplitude = 2 * (exp(1) * gamma_parameter / beta_parameter)^(
        beta_parameter / gamma_parameter
    )
    return amplitude * normalized_frequency^beta_parameter *
           exp(-normalized_frequency^gamma_parameter)
end

"""
    admissibility_constant(wavelet; method=:analytic)

Return the positive-frequency direct-synthesis constant
`integral(conj(wavelet_hat(u)) / u, u=0..Inf)` for the unnormalized analytic
Morse response used by the standard backend adapter. `method=:quadrature`
evaluates the same integral independently on a refined logarithmic grid.
"""
function admissibility_constant(
    wavelet::_CW.Morse;
    method::Symbol = :analytic,
    points::Integer = 65_537,
)
    wavelet.ga > 0 || throw(ArgumentError("Morse gamma must be positive for admissibility"))
    wavelet.be > 0 || throw(ArgumentError("Morse beta must be positive for admissibility"))
    wavelet.cf > 0 || throw(ArgumentError(
        "Morse center frequency must be positive for admissibility",
    ))
    if method === :analytic
        ratio = wavelet.be / wavelet.ga
        amplitude = 2 * (exp(1) * wavelet.ga / wavelet.be)^ratio
        return amplitude * gamma(ratio) / wavelet.ga
    elseif method === :quadrature
        points >= 257 && isodd(points) || throw(ArgumentError(
            "quadrature points must be an odd integer of at least 257",
        ))
        center = log(wavelet.cf)
        log_frequencies = range(center - 40, center + 8; length = points)
        step = Base.step(log_frequencies)
        integral = zero(Float64)
        for (index, log_frequency) in enumerate(log_frequencies)
            value = _morse_response(wavelet, exp(log_frequency))
            weight = index == 1 || index == points ? 0.5 : 1.0
            integral += weight * value
        end
        return step * integral
    end
    throw(ArgumentError("admissibility method must be :analytic or :quadrature"))
end

"""Undo the backend's per-filter peak normalization for an analytic Morse bank."""
function _restore_morse_amplitudes!(daughters, scales, omega, wavelet::_CW.Morse)
    for index in eachindex(scales)
        raw_peak = maximum(frequency -> _morse_response(wavelet, frequency / scales[index]), omega)
        @views daughters[:, index] .*= raw_peak
    end
    return daughters
end

function _reliable_frequency_range(response::Vector{T}, admissibility::Complex{T}, fs::T, n::Int;
    tolerance::T = T(0.01),
) where {T<:AbstractFloat}
    relative_error = abs.(response ./ admissibility .- one(T))
    reliable = relative_error .<= tolerance
    reliable[1] = false
    best_start = 0
    best_stop = 0
    current_start = 0
    for index in eachindex(reliable)
        if reliable[index]
            iszero(current_start) && (current_start = index)
        elseif !iszero(current_start)
            current_stop = index - 1
            if current_stop - current_start > best_stop - best_start
                best_start = current_start
                best_stop = current_stop
            end
            current_start = 0
        end
    end
    if !iszero(current_start) && length(reliable) - current_start > best_stop - best_start
        best_start = current_start
        best_stop = length(reliable)
    end
    iszero(best_start) && throw(ArgumentError(
        "the finite wavelet bank has no contiguous frequency region within $(tolerance) synthesis error",
    ))
    bin_width = fs / T(n)
    return (T(best_start - 1) * bin_width, T(best_stop - 1) * bin_width)
end

"""Fourier response of summing the quadrature-weighted backend coefficients."""
function _synthesis_response(daughters::AbstractMatrix, weights::Vector{T}) where {T<:AbstractFloat}
    response = vec(daughters * weights)
    maximum(abs, imag.(response)) <= sqrt(eps(T)) * max(maximum(abs, response), one(T)) ||
        throw(ArgumentError("the analytic wavelet-bank response must be real"))
    return T.(real.(response))
end

"""
    plan_cwt(n; fs, wavelet=nothing, voices_per_octave=16, standard=true,
             beta=2, precision=Float64, dc_tolerance=0.01)

Build a periodic `ContinuousWavelets.jl` plan. A periodic boundary is required
in this milestone because it makes the FFT-based time derivative exact under
the recorded transform convention.
"""
function plan_cwt(
    n::Integer;
    fs::Real,
    wavelet = nothing,
    voices_per_octave::Integer = 16,
    standard::Bool = true,
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

    selected_wavelet = isnothing(wavelet) ?
        (standard ? _DEFAULT_STANDARD_WAVELET : _CW.Morlet(π)) : wavelet
    standard && !(selected_wavelet isa _CW.Morse) && throw(ArgumentError(
        "standard WSST reconstruction currently supports ContinuousWavelets.Morse wavelets",
    ))
    maximum_scale_steps = max(1, floor(Int, voices_per_octave * log2(max(2, n ÷ 2))) - 1)
    backend_config = if standard
        _CW.wavelet(
            selected_wavelet;
            Q = voices_per_octave,
            J = maximum_scale_steps,
            p = Inf,
            frameBound = 0,
            fsample = T(fs),
            averagingType = _CW.NoAve(),
            boundary = _continuous_wavelets_boundary(boundary),
        )
    else
        _CW.wavelet(
            selected_wavelet;
            Q = voices_per_octave,
            β = beta,
            fsample = T(fs),
            averagingType = _CW.NoAve(),
            boundary = _continuous_wavelets_boundary(boundary),
        )
    end
    daughters, omega = _CW.computeWavelets(n, backend_config; T = T)
    _, _, scales, _ = _CW.getNWavelets(n, backend_config)
    standard && _restore_morse_amplitudes!(daughters, scales, omega, selected_wavelet)
    frequencies = if standard
        T.(selected_wavelet.cf .* scales .* (fs / n))
    else
        T.(_CW.getMeanFreq(daughters, T(fs)))
    end
    frequency_weights = _frequency_weights(frequencies)
    validation = _validate_wavelets(daughters, T(dc_tolerance))
    reconstruction_weights = _reconstruction_weights(daughters, backend_config, inverse_style, T)
    squeezing_weights = _log_frequency_weights(frequencies)
    synthesis_response = _synthesis_response(daughters, squeezing_weights)
    admissibility = standard ? Complex{T}(admissibility_constant(selected_wavelet)) : Complex{T}(NaN)
    reliable_frequency_range = standard ?
        _reliable_frequency_range(synthesis_response, admissibility, T(fs), Int(n)) :
        (first(frequencies), last(frequencies))

    return CWTPlan(
        ContinuousWaveletsBackend(),
        backend_config,
        daughters,
        T.(scales),
        frequencies,
        frequency_weights,
        inverse_style,
        reconstruction_weights,
        squeezing_weights,
        synthesis_response,
        admissibility,
        standard,
        reliable_frequency_range,
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
function _metadata(plan::CWTPlan{ContinuousWaveletsBackend}, signal_kind::Symbol)
    return CWTMetadata(
        plan.sample_rate,
        inv(plan.sample_rate),
        plan.signal_length,
        plan.signal_length,
        1:plan.signal_length,
        :periodic,
        plan.standard_compatible ? :analytic_l1 : :continuouswavelets_linf_frame,
        signal_kind,
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
    all(isfinite, signal) || throw(ArgumentError("signal samples must be finite"))

    signal_kind = _signal_kind(signal)
    dc_component = sum(signal) / length(signal)
    centered_signal = signal .- dc_component
    backend_coefficients = _CW.cwt(centered_signal, plan.backend_config, plan.daughters)
    backend_derivative = _spectral_derivative(centered_signal, plan)
    coefficients = PermutedDimsArray(backend_coefficients, (2, 1))
    derivative = PermutedDimsArray(backend_derivative, (2, 1))
    reliability = repeat(plan.validation.valid_rows, 1, plan.signal_length)

    return CWTResult(
        coefficients,
        derivative,
        copy(plan.frequencies),
        copy(plan.scales),
        copy(plan.frequency_weights),
        isnothing(plan.reconstruction_weights) ? nothing : copy(plan.reconstruction_weights),
        copy(plan.squeezing_weights),
        copy(plan.synthesis_response),
        plan.admissibility,
        plan.standard_compatible,
        plan.reliable_frequency_range,
        copy(plan.validation.valid_rows),
        reliability,
        _metadata(plan, signal_kind),
        dc_component,
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
    result.metadata.signal_kind === :two_sided_complex && throw(ArgumentError(
        "the one-sided analytic CWT cannot reconstruct a two-sided complex signal",
    ))
    isnothing(result.reconstruction_weights) && throw(ArgumentError(
        "CWT inversion requires Penrose reconstruction weights",
    ))
    weighted = vec(sum(
        result.coefficients .* reshape(result.reconstruction_weights, :, 1);
        dims = 1,
    ))
    result.metadata.signal_kind === :analytic_complex && return weighted .+ result.dc_component
    quadrature = irfft(im .* rfft(imag.(weighted)), length(weighted))
    return real.(weighted) .+ quadrature .+ result.dc_component
end

"""
    instantaneous_frequency(result; threshold=0)

Return the phase-based physical frequency estimate in Hz. Values below the
coefficient threshold or on analytically unsuitable rows are marked invalid and
stored as `NaN`; they are never replaced by zero.
"""
function instantaneous_frequency(result::CWTResult; threshold = NoThreshold())
    return instantaneous_frequency(
        result.coefficients,
        result.derivative;
        threshold,
        reliability = result.reliability,
    )
end
