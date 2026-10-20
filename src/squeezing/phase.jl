"""Resolve a threshold policy to one absolute coefficient magnitude."""
_threshold_magnitude(::NoThreshold, coefficients) = zero(_realtype(eltype(coefficients)))

function _threshold_magnitude(threshold::AbsoluteThreshold, coefficients)
    threshold.value >= zero(threshold.value) || throw(ArgumentError("absolute threshold must be nonnegative"))
    return threshold.value
end

function _threshold_magnitude(threshold::RelativeThreshold, coefficients)
    threshold.fraction >= zero(threshold.fraction) || throw(ArgumentError("relative threshold must be nonnegative"))
    threshold.reference === :global || throw(ArgumentError(
        "only reference=:global is currently supported for relative thresholds",
    ))
    T = _realtype(eltype(coefficients))
    reference = zero(T)
    for coefficient in coefficients
        magnitude = abs(coefficient)
        isfinite(magnitude) && (reference = max(reference, magnitude))
    end
    return threshold.fraction * reference
end

_threshold_policy(value::AbstractMagnitudeThreshold, ::Type) = value
_threshold_policy(value::Real, ::Type{T}) where {T<:AbstractFloat} = AbsoluteThreshold(T(value))

"""Stable phase ratio for one finite coefficient pair, returned in physical Hz."""
function _phase_frequency(coefficient::Complex{T}, derivative::Complex{T}) where {T<:AbstractFloat}
    scale = max(abs(real(coefficient)), abs(imag(coefficient)))
    isfinite(scale) && scale > zero(T) || return T(NaN)
    normalized_coefficient = coefficient / scale
    normalized_derivative = derivative / scale
    denominator = abs2(normalized_coefficient)
    numerator = imag(normalized_derivative * conj(normalized_coefficient))
    frequency = numerator / (T(2π) * denominator)
    return isfinite(frequency) ? frequency : T(NaN)
end

"""
    instantaneous_frequency(coefficients, derivative; threshold=NoThreshold(), reliability=nothing)

Estimate physical frequency in Hz with a separate validity mask. Relative
thresholds use the global finite maximum coefficient magnitude. Invalid values
are stored as `NaN`, never as zero frequency.
"""
function instantaneous_frequency(
    coefficients::AbstractMatrix{Complex{T}},
    derivative::AbstractMatrix{Complex{T}};
    threshold = NoThreshold(),
    reliability = nothing,
) where {T<:AbstractFloat}
    size(coefficients) == size(derivative) || throw(DimensionMismatch(
        "coefficient and derivative matrices must have identical dimensions",
    ))
    reliability === nothing || size(reliability) == size(coefficients) || throw(DimensionMismatch(
        "reliability mask must match the coefficient dimensions",
    ))
    policy = _threshold_policy(threshold, T)
    magnitude_threshold = T(_threshold_magnitude(policy, coefficients))
    values = fill(T(NaN), size(coefficients))
    valid = falses(size(coefficients))

    for index in eachindex(coefficients, derivative)
        reliable = reliability === nothing || reliability[index]
        coefficient = coefficients[index]
        differentiated = derivative[index]
        if reliable && isfinite(coefficient) && isfinite(differentiated) &&
           abs(coefficient) > magnitude_threshold
            frequency = _phase_frequency(coefficient, differentiated)
            if isfinite(frequency)
                values[index] = frequency
                valid[index] = true
            end
        end
    end
    return FrequencyEstimate(values, valid)
end
