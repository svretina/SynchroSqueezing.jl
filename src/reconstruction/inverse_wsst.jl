"""Apply standard scalar analytic-wavelet synthesis to one scale integral."""
function _standard_synthesis(
    integrated::AbstractVector{Complex{T}},
    convention::ReconstructionConvention{T},
) where {T<:AbstractFloat}
    convention.method === :standard_analytic_l1 || throw(ArgumentError(
        "standard inversion requires an L1-normalized analytic CWT result",
    ))
    isfinite(convention.admissibility) && !iszero(convention.admissibility) ||
        throw(ArgumentError("standard inversion requires a finite nonzero admissibility constant"))
    normalized = integrated ./ convention.admissibility
    if convention.signal_kind === :real
        return T(2) .* real.(normalized)
    elseif convention.signal_kind === :analytic_complex
        return normalized
    end
    throw(ArgumentError("unsupported reconstruction signal kind $(convention.signal_kind)"))
end

"""Apply standard synthesis while retaining the positive-frequency analytic signal."""
function _analytic_standard_synthesis(
    integrated::AbstractVector{Complex{T}},
    convention::ReconstructionConvention{T},
) where {T<:AbstractFloat}
    convention.method === :standard_analytic_l1 || throw(ArgumentError(
        "standard inversion requires an L1-normalized analytic CWT result",
    ))
    isfinite(convention.admissibility) && !iszero(convention.admissibility) ||
        throw(ArgumentError("standard inversion requires a finite nonzero admissibility constant"))
    normalized = integrated ./ convention.admissibility
    if convention.signal_kind === :real
        return T(2) .* normalized
    elseif convention.signal_kind === :analytic_complex
        return normalized
    end
    throw(ArgumentError("unsupported reconstruction signal kind $(convention.signal_kind)"))
end

"""Apply the retained backend-derived diagonal frame synthesis operator."""
function _frame_synthesis(
    weighted::AbstractVector{Complex{T}},
    convention::ReconstructionConvention{T},
) where {T<:AbstractFloat}
    n = length(weighted)
    spectrum = fft(weighted)
    npositive = n ÷ 2 + 1
    length(convention.synthesis_response) == npositive || throw(DimensionMismatch(
        "the stored synthesis response does not match the signal length",
    ))
    @inbounds for index in 1:npositive
        response = convention.synthesis_response[index]
        spectrum[index] = abs(response) > convention.response_floor ?
            spectrum[index] / response : zero(eltype(spectrum))
    end
    @views fill!(spectrum[npositive+1:end], zero(eltype(spectrum)))
    spectrum[1] = zero(eltype(spectrum))
    convention.signal_kind === :real && iseven(n) && (spectrum[npositive] /= T(2))
    analytic = ifft(spectrum)
    if convention.signal_kind === :analytic_complex
        return analytic
    elseif convention.signal_kind === :real
        return T(2) .* real.(analytic)
    end
    throw(ArgumentError("unsupported reconstruction signal kind $(convention.signal_kind)"))
end

function _synthesize(integrated, convention, method::Symbol)
    method === :standard && return _standard_synthesis(integrated, convention)
    method === :frame && return _frame_synthesis(integrated, convention)
    throw(ArgumentError("reconstruction method must be :standard or :frame"))
end

function _sum_masked(coefficients::AbstractMatrix{C}, mask::AbstractMatrix{Bool}) where {C<:Complex}
    size(mask) == size(coefficients) || throw(DimensionMismatch(
        "the reconstruction mask must match the WSST coefficient matrix",
    ))
    integrated = zeros(C, size(coefficients, 2))
    for column in axes(coefficients, 2), row in axes(coefficients, 1)
        mask[row, column] && (integrated[column] += coefficients[row, column])
    end
    return integrated
end

"""
    reconstruct_region(result, mask; method=:standard)

Reconstruct an arbitrary time-frequency region using only the selected,
bin-integrated WSST coefficients and stored reconstruction metadata. The mask
must have the same frequency × time shape as `result.coefficients`. DC and
excluded-coefficient residuals are never included.
"""
function reconstruct_region(
    result::WSSTResult,
    mask::AbstractMatrix{Bool};
    method::Symbol = :standard,
)
    method === :standard || throw(ArgumentError(
        "masked component reconstruction is defined only for method=:standard",
    ))
    return _standard_synthesis(_sum_masked(result.coefficients, mask), result.reconstruction)
end

"""
    inverse_wsst(result; mode=:oscillatory, method=:standard)

Reconstruct the deposited oscillatory coefficients with the standard scalar
admissibility formula by default. `mode=:full` explicitly adds retained,
excluded weighted-CWT contributions and DC; `mode=:residual` returns that
complement alone. `method=:frame` retains the finite-bank diagonal correction
from Milestone 03 for whole-result diagnostics, but is not used for masked
component reconstruction.
"""
function inverse_wsst(
    result::WSSTResult;
    mode::Symbol = :oscillatory,
    method::Symbol = :standard,
)
    mode in (:oscillatory, :residual, :full) || throw(ArgumentError(
        "mode must be :oscillatory, :residual, or :full",
    ))
    deposited = vec(sum(result.coefficients; dims = 1))
    mode === :oscillatory && return _synthesize(deposited, result.reconstruction, method)

    result.residual === nothing && throw(ArgumentError(
        "$(mode) inversion requires a result computed with preserve_residual=true",
    ))
    any(>(0), result.diagnostics.unrecoverable) && throw(ArgumentError(
        "$(mode) inversion is unavailable because nonfinite CWT coefficients were omitted",
    ))
    if method === :frame
        all(response -> abs(response) > result.reconstruction.response_floor,
            @view(result.reconstruction.synthesis_response[2:end])) || throw(ArgumentError(
            "frame inversion requires wavelet-bank coverage at every non-DC Fourier bin",
        ))
    end
    if mode === :residual
        return _synthesize(result.residual, result.reconstruction, method) .+
               result.dc_component
    end
    return _synthesize(deposited .+ result.residual, result.reconstruction, method) .+
           result.dc_component
end
