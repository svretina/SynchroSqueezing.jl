"""Apply a positive-frequency filter bank and a per-bin spectral multiplier."""
function _positive_filter_transform(
    spectrum::AbstractVector,
    filters::AbstractMatrix,
    multiplier::AbstractVector,
    n::Int,
)
    nfrequency, nscales = size(filters)
    length(spectrum) == nfrequency || throw(DimensionMismatch(
        "positive spectrum and filter bank have incompatible lengths",
    ))
    length(multiplier) == nfrequency || throw(DimensionMismatch(
        "spectral multiplier and filter bank have incompatible lengths",
    ))
    C = promote_type(eltype(spectrum), eltype(filters), eltype(multiplier))
    transformed_spectrum = zeros(C, n, nscales)
    @views transformed_spectrum[1:nfrequency, :] .=
        reshape(spectrum, :, 1) .* filters .* reshape(multiplier, :, 1)
    return PermutedDimsArray(ifft(transformed_spectrum, 1), (2, 1))
end

"""Derivative of the sampled Morse daughter bank with respect to physical Hz."""
function _morse_frequency_derivative_bank(
    plan::CWTPlan{ContinuousWaveletsBackend},
    ::Type{T},
) where {T<:AbstractFloat}
    plan.standard_compatible || throw(ArgumentError(
        "second-order WSST currently requires the standard analytic Morse backend",
    ))
    wavelet = plan.backend_config.waveType
    wavelet isa _CW.Morse || throw(ArgumentError(
        "second-order WSST currently supports ContinuousWavelets.Morse wavelets",
    ))

    nfrequency = size(plan.daughters, 1)
    backend_frequency = T.(_CW.computeOmega(plan.signal_length, nfrequency))
    physical_frequency = _rfft_angular_frequencies(plan.signal_length, plan.sample_rate) ./ T(2π)
    bank = Matrix{T}(undef, size(plan.daughters))
    gamma_parameter = T(wavelet.ga)
    beta_parameter = T(wavelet.be)
    peak_frequency = (beta_parameter / gamma_parameter)^(inv(gamma_parameter))
    frequency_factor = T(wavelet.cf) / peak_frequency

    for scale_index in axes(bank, 2), frequency_index in axes(bank, 1)
        frequency = physical_frequency[frequency_index]
        if iszero(frequency)
            bank[frequency_index, scale_index] = zero(T)
            continue
        end
        normalized = (backend_frequency[frequency_index] /
                      T(plan.scales[scale_index])) / frequency_factor
        response = T(plan.daughters[frequency_index, scale_index])
        bank[frequency_index, scale_index] = response *
            (beta_parameter - gamma_parameter * normalized^gamma_parameter) / frequency
    end
    return bank
end

"""
    second_order_operators(plan, signal)

Compute the auxiliary transforms used by the time-derivative WSST2 formula of
Pham and Meignen. The implementation uses the exact sampled Morse daughter
bank and spectral differentiation; it does not differentiate an extracted
ridge or use finite differences in time.
"""
function second_order_operators(
    plan::CWTPlan{ContinuousWaveletsBackend},
    signal::AbstractVector,
)
    length(signal) == plan.signal_length || throw(DimensionMismatch(
        "plan length $(plan.signal_length) does not match signal length $(length(signal))",
    ))
    T = _realtype(eltype(signal))
    T === typeof(plan.sample_rate) || throw(ArgumentError(
        "signal precision $T does not match plan precision $(typeof(plan.sample_rate))",
    ))
    all(isfinite, signal) || throw(ArgumentError("signal samples must be finite"))

    transformed = cwt(plan, signal)
    centered = signal .- transformed.dc_component
    spectrum = _positive_spectrum(centered)
    angular_frequency = _rfft_angular_frequencies(plan.signal_length, plan.sample_rate)
    daughters = T.(plan.daughters)
    frequency_derivative = _morse_frequency_derivative_bank(plan, T)
    ones_multiplier = ones(T, length(angular_frequency))

    second_derivative = _positive_filter_transform(
        spectrum,
        daughters,
        .-(angular_frequency .^ 2),
        plan.signal_length,
    )
    moment_transform = _positive_filter_transform(
        spectrum,
        frequency_derivative,
        ones_multiplier,
        plan.signal_length,
    ) ./ Complex{T}(0, T(2π))
    moment_derivative = _positive_filter_transform(
        spectrum,
        frequency_derivative,
        Complex{T}.(zero(T), angular_frequency),
        plan.signal_length,
    ) ./ Complex{T}(0, T(2π))

    return SecondOrderOperators(
        transformed,
        second_derivative,
        moment_transform,
        moment_derivative,
    )
end

"""
    second_order_frequency(operators; threshold=RelativeThreshold(1e-6),
                           denominator_threshold=sqrt(eps(T)),
                           fallback=:first_order,
                           max_abs_chirp_rate=Inf)

Estimate instantaneous frequency and chirp rate from the complex WSST2
operators. `fallback=:first_order` uses the first-order estimate wherever the
second-order correction is invalid, including an ill-conditioned
time-reassignment derivative;
`fallback=:exclude` rejects those coefficients.
"""
function second_order_frequency(
    operators::SecondOrderOperators;
    threshold = nothing,
    denominator_threshold = nothing,
    fallback::Symbol = :first_order,
    max_abs_chirp_rate::Real = Inf,
)
    T = eltype(operators.cwt.frequencies)
    threshold === nothing && (threshold = RelativeThreshold(T(1e-6)))
    threshold isa AbstractMagnitudeThreshold || throw(ArgumentError(
        "threshold must be an AbstractMagnitudeThreshold",
    ))
    denominator_threshold === nothing && (denominator_threshold = sqrt(eps(T)))
    fallback in (:first_order, :exclude) || throw(ArgumentError(
        "second-order fallback must be :first_order or :exclude",
    ))
    isfinite(denominator_threshold) && denominator_threshold >= 0 || throw(ArgumentError(
        "denominator_threshold must be finite and nonnegative",
    ))
    (isfinite(max_abs_chirp_rate) || max_abs_chirp_rate == Inf) &&
        max_abs_chirp_rate >= 0 || throw(ArgumentError(
            "max_abs_chirp_rate must be nonnegative and not NaN",
        ))

    cwt_result = operators.cwt
    coefficients = cwt_result.coefficients
    derivative = cwt_result.derivative
    size(operators.second_derivative) == size(coefficients) || throw(DimensionMismatch(
        "second time derivative must match CWT dimensions",
    ))
    size(operators.time_moment) == size(coefficients) || throw(DimensionMismatch(
        "time-moment transform must match CWT dimensions",
    ))
    size(operators.time_moment_derivative) == size(coefficients) || throw(DimensionMismatch(
        "time-moment derivative must match CWT dimensions",
    ))

    magnitude_threshold = T(_threshold_magnitude(threshold, coefficients))
    denominator_floor = T(denominator_threshold)
    chirp_limit = T(max_abs_chirp_rate)
    values = fill(T(NaN), size(coefficients))
    first_order = fill(T(NaN), size(coefficients))
    chirp_rate = fill(T(NaN), size(coefficients))
    nan_complex = Complex{T}(T(NaN), T(NaN))
    time_reassignment = fill(nan_complex, size(coefficients))
    time_derivative = fill(nan_complex, size(coefficients))
    valid = falses(size(coefficients))
    second_order_valid = falses(size(coefficients))
    used_fallback = falses(size(coefficients))
    inverse_i2pi = inv(Complex{T}(0, T(2π)))
    dt = cwt_result.metadata.sample_interval

    for column in axes(coefficients, 2), row in axes(coefficients, 1)
        cwt_result.reliability[row, column] || continue
        w = coefficients[row, column]
        wt = derivative[row, column]
        wtt = operators.second_derivative[row, column]
        moment = operators.time_moment[row, column]
        moment_t = operators.time_moment_derivative[row, column]
        all(isfinite, (w, wt, wtt, moment, moment_t)) || continue
        abs(w) > magnitude_threshold || continue

        first_ratio = wt / w
        first_frequency = real(inverse_i2pi * first_ratio)
        isfinite(first_frequency) || continue
        first_order[row, column] = first_frequency

        moment_ratio = moment / w
        tau = T(column - 1) * dt + moment_ratio
        tau_derivative = one(T) + moment_t / w - moment_ratio * first_ratio
        time_reassignment[row, column] = tau
        time_derivative[row, column] = tau_derivative

        second_valid = abs(tau_derivative) > denominator_floor
        if second_valid
            omega_derivative = inverse_i2pi * (wtt / w - first_ratio^2)
            chirp = omega_derivative / tau_derivative
            corrected = inverse_i2pi * first_ratio - chirp * moment_ratio
            second_frequency = real(corrected)
            real_chirp = real(chirp)
            second_valid = isfinite(chirp) && isfinite(second_frequency) &&
                           second_frequency > zero(T) &&
                           abs(real_chirp) <= chirp_limit
            if second_valid
                values[row, column] = second_frequency
                chirp_rate[row, column] = real_chirp
                valid[row, column] = true
                second_order_valid[row, column] = true
            end
        end

        if !second_valid && fallback === :first_order && first_frequency > zero(T)
            values[row, column] = first_frequency
            valid[row, column] = true
            used_fallback[row, column] = true
        end
    end

    return SecondOrderEstimate(
        values,
        valid,
        second_order_valid,
        used_fallback,
        first_order,
        chirp_rate,
        time_reassignment,
        time_derivative,
    )
end
