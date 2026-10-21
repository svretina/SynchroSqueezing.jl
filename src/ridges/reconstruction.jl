function _validate_ridge_for_result(result::WSSTResult, ridge::RidgeResult)
    ntime = size(result.coefficients, 2)
    length(ridge.indices) == ntime || throw(DimensionMismatch(
        "ridge and WSST result must have the same number of time samples",
    ))
    for column in 1:ntime
        ridge.valid[column] || continue
        row = ridge.indices[column]
        1 <= row <= length(result.frequencies) || throw(ArgumentError(
            "ridge index at time $column is outside the WSST grid",
        ))
        result.frequencies[row] == ridge.frequencies[column] || throw(ArgumentError(
            "ridge frequency coordinates do not match the WSST result",
        ))
    end
    return nothing
end

function _select_bandwidth(bandwidth_hz, bandwidth_bins, bandwidth_relative, custom_mask)
    supplied = count(!isnothing, (
        bandwidth_hz,
        bandwidth_bins,
        bandwidth_relative,
        custom_mask,
    ))
    supplied == 1 || throw(ArgumentError(
        "provide exactly one of bandwidth_hz, bandwidth_bins, bandwidth_relative, or custom_mask",
    ))
    return nothing
end

"""
    ridge_mask(result, ridge; bandwidth_hz, bandwidth_bins,
               bandwidth_relative, custom_mask)

Construct an inclusive frequency-by-time reconstruction mask. Only active,
valid ridge samples are selected; inactive trajectory estimates are not
silently synthesized.
"""
function ridge_mask(
    result::WSSTResult,
    ridge::RidgeResult;
    bandwidth_hz = nothing,
    bandwidth_bins = nothing,
    bandwidth_relative = nothing,
    custom_mask = nothing,
)
    _validate_ridge_for_result(result, ridge)
    _select_bandwidth(bandwidth_hz, bandwidth_bins, bandwidth_relative, custom_mask)
    if custom_mask !== nothing
        size(custom_mask) == size(result.coefficients) || throw(DimensionMismatch(
            "custom_mask must match the WSST coefficient matrix",
        ))
        return BitMatrix(custom_mask)
    end
    mask = falses(size(result.coefficients))
    if bandwidth_hz !== nothing
        isfinite(bandwidth_hz) && bandwidth_hz >= 0 || throw(ArgumentError(
            "bandwidth_hz must be finite and nonnegative",
        ))
        for column in axes(mask, 2)
            ridge.valid[column] && ridge.active[column] || continue
            center = ridge.frequencies[column]
            for row in axes(mask, 1)
                mask[row, column] = abs(result.frequencies[row] - center) <= bandwidth_hz
            end
        end
    elseif bandwidth_bins !== nothing
        bandwidth_bins isa Integer && bandwidth_bins >= 0 || throw(ArgumentError(
            "bandwidth_bins must be a nonnegative integer",
        ))
        for column in axes(mask, 2)
            ridge.valid[column] && ridge.active[column] || continue
            center = ridge.indices[column]
            lower = max(first(axes(mask, 1)), center - bandwidth_bins)
            upper = min(last(axes(mask, 1)), center + bandwidth_bins)
            mask[lower:upper, column] .= true
        end
    else
        isfinite(bandwidth_relative) && bandwidth_relative >= 0 || throw(ArgumentError(
            "bandwidth_relative must be finite and nonnegative",
        ))
        for column in axes(mask, 2)
            ridge.valid[column] && ridge.active[column] || continue
            center = ridge.frequencies[column]
            lower = center * exp(-bandwidth_relative)
            upper = center * exp(bandwidth_relative)
            for row in axes(mask, 1)
                mask[row, column] = lower <= result.frequencies[row] <= upper
            end
        end
    end
    return mask
end

function _sum_weighted(coefficients::AbstractMatrix{C}, weights::AbstractMatrix) where {C<:Complex}
    size(weights) == size(coefficients) || throw(DimensionMismatch(
        "ownership weights must match the WSST coefficient matrix",
    ))
    integrated = zeros(C, size(coefficients, 2))
    for column in axes(coefficients, 2), row in axes(coefficients, 1)
        integrated[column] += weights[row, column] * coefficients[row, column]
    end
    return integrated
end

"""Reconstruct one component from an independently configured ridge band.

Set `analytic=true` to retain the complex positive-frequency reconstruction,
including when the original signal was real. Its magnitude is the component
envelope and its real part matches the ordinary real reconstruction.
"""
function reconstruct_component(
    result::WSSTResult,
    ridge::RidgeResult;
    analytic::Bool = false,
    kwargs...,
)
    mask = ridge_mask(result, ridge; kwargs...)
    T = eltype(result.frequencies)
    weights = T.(mask)
    integrated = _sum_weighted(result.coefficients, weights)
    signal = analytic ?
        _analytic_standard_synthesis(integrated, result.reconstruction) :
        _standard_synthesis(integrated, result.reconstruction)
    return ComponentReconstruction(signal, mask, weights)
end

function _component_masks(
    result,
    ridges;
    bandwidth_hz,
    bandwidth_bins,
    bandwidth_relative,
    custom_masks,
)
    count(!isnothing, (bandwidth_hz, bandwidth_bins, bandwidth_relative, custom_masks)) == 1 ||
        throw(ArgumentError(
            "provide exactly one reconstruction-band specification",
        ))
    ncomponents = length(ridges)
    expand(parameter, name) = if parameter isa AbstractVector
        length(parameter) == ncomponents || throw(DimensionMismatch(
            "$name must contain one value per ridge",
        ))
        collect(parameter)
    else
        fill(parameter, ncomponents)
    end
    if custom_masks !== nothing
        length(custom_masks) == ncomponents || throw(DimensionMismatch(
            "custom_masks must contain one mask per ridge",
        ))
        return [ridge_mask(result, ridges[index]; custom_mask = custom_masks[index])
                for index in 1:ncomponents]
    elseif bandwidth_hz !== nothing
        widths = expand(bandwidth_hz, "bandwidth_hz")
        return [ridge_mask(result, ridges[index]; bandwidth_hz = widths[index])
                for index in 1:ncomponents]
    elseif bandwidth_bins !== nothing
        widths = expand(bandwidth_bins, "bandwidth_bins")
        return [ridge_mask(result, ridges[index]; bandwidth_bins = widths[index])
                for index in 1:ncomponents]
    else
        widths = expand(bandwidth_relative, "bandwidth_relative")
        return [ridge_mask(result, ridges[index]; bandwidth_relative = widths[index])
                for index in 1:ncomponents]
    end
end

function _exclusive_weights(result, ridges, masks, ::Type{T}) where {T}
    weights = [zeros(T, size(result.coefficients)) for _ in ridges]
    for column in axes(result.coefficients, 2), row in axes(result.coefficients, 1)
        owner = 0
        best_distance = T(Inf)
        for component in eachindex(ridges)
            masks[component][row, column] || continue
            distance = abs(result.frequencies[row] - ridges[component].frequencies[column])
            if distance < best_distance
                owner = component
                best_distance = distance
            end
        end
        owner == 0 || (weights[owner][row, column] = one(T))
    end
    return weights
end

function _fractional_weights(result, ridges, masks, ::Type{T}) where {T}
    weights = [zeros(T, size(result.coefficients)) for _ in ridges]
    distances = zeros(T, length(ridges))
    for column in axes(result.coefficients, 2), row in axes(result.coefficients, 1)
        eligible_components = Int[]
        for component in eachindex(ridges)
            masks[component][row, column] || continue
            push!(eligible_components, component)
            distances[component] = abs(
                result.frequencies[row] - ridges[component].frequencies[column],
            )
        end
        isempty(eligible_components) && continue
        exact = filter(component -> iszero(distances[component]), eligible_components)
        if !isempty(exact)
            share = inv(T(length(exact)))
            for component in exact
                weights[component][row, column] = share
            end
        else
            inverse_sum = sum(component -> inv(distances[component]), eligible_components)
            for component in eligible_components
                weights[component][row, column] =
                    inv(distances[component]) / inverse_sum
            end
        end
    end
    return weights
end

"""
    reconstruct_components(result, ridges; ownership=:exclusive, ...)

Reconstruct multiple ridge bands without double counting. Ownership is either
deterministic nearest-ridge `:exclusive` assignment or inverse-distance
`:fractional` partitioning. The returned residual contains every unowned SST
coefficient, so component plus residual synthesis equals the full represented
oscillatory synthesis.
"""
function reconstruct_components(
    result::WSSTResult,
    ridges::AbstractVector{<:RidgeResult};
    bandwidth_hz = nothing,
    bandwidth_bins = nothing,
    bandwidth_relative = nothing,
    custom_masks = nothing,
    ownership::Symbol = :exclusive,
    analytic::Bool = false,
)
    isempty(ridges) && throw(ArgumentError("at least one ridge is required"))
    ownership in (:exclusive, :fractional) || throw(ArgumentError(
        "ownership must be :exclusive or :fractional",
    ))
    masks = _component_masks(
        result,
        ridges;
        bandwidth_hz,
        bandwidth_bins,
        bandwidth_relative,
        custom_masks,
    )
    T = eltype(result.frequencies)
    weights = ownership === :exclusive ?
        _exclusive_weights(result, ridges, masks, T) :
        _fractional_weights(result, ridges, masks, T)
    residual_weights = ones(T, size(result.coefficients))
    for component_weights in weights
        residual_weights .-= component_weights
    end
    tolerance = T(32) * eps(T)
    all(weight -> weight >= -tolerance, residual_weights) || error(
        "component ownership exceeded unity",
    )
    residual_weights .= max.(residual_weights, zero(T))
    synthesis = analytic ? _analytic_standard_synthesis : _standard_synthesis
    components = [synthesis(
        _sum_weighted(result.coefficients, component_weights),
        result.reconstruction,
    ) for component_weights in weights]
    residual = synthesis(
        _sum_weighted(result.coefficients, residual_weights), result.reconstruction,
    )
    return ComponentSet(components, masks, weights, residual, residual_weights, ownership)
end
