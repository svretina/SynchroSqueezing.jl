"""
    linear_frequency_grid(fmin, fmax, nbins; precision=Float64)

Return `nbins` ascending positive-frequency bin centers in Hz, including both
endpoints. The matching bin edges are available internally to the reassignment
kernel and are stored in every `WSSTResult`.
"""
function linear_frequency_grid(
    fmin::Real,
    fmax::Real,
    nbins::Integer;
    precision::Type{T} = Float64,
) where {T<:AbstractFloat}
    nbins > 1 || throw(ArgumentError("nbins must exceed one"))
    fmin > 0 || throw(ArgumentError("fmin must be positive"))
    fmax > fmin || throw(ArgumentError("fmax must exceed fmin"))
    return collect(range(T(fmin), T(fmax); length = nbins))
end

"""
    log_frequency_grid(fmin, fmax, nbins; precision=Float64)

Return `nbins` geometrically spaced positive-frequency centers in Hz,
including both endpoints.
"""
function log_frequency_grid(
    fmin::Real,
    fmax::Real,
    nbins::Integer;
    precision::Type{T} = Float64,
) where {T<:AbstractFloat}
    nbins > 1 || throw(ArgumentError("nbins must exceed one"))
    fmin > 0 || throw(ArgumentError("fmin must be positive"))
    fmax > fmin || throw(ArgumentError("fmax must exceed fmin"))
    return exp.(collect(range(log(T(fmin)), log(T(fmax)); length = nbins)))
end

"""Return outer half-step edges for a strictly ascending frequency grid."""
function _frequency_bin_edges(
    frequencies::Vector{T};
    spacing::Symbol = :arbitrary,
) where {T<:AbstractFloat}
    length(frequencies) > 1 || throw(ArgumentError("at least two output frequencies are required"))
    all(isfinite, frequencies) && all(>(zero(T)), frequencies) ||
        throw(ArgumentError("output frequencies must be finite and positive"))
    steps = diff(frequencies)
    all(>(zero(T)), steps) || throw(ArgumentError("output frequencies must be strictly ascending"))

    spacing in (:linear, :log, :arbitrary) || throw(ArgumentError(
        "frequency_spacing must be :linear, :log, or :arbitrary",
    ))
    edges = Vector{T}(undef, length(frequencies) + 1)
    if spacing === :log
        ratios = frequencies[2:end] ./ frequencies[1:end-1]
        maximum(abs.(ratios .- first(ratios))) <= sqrt(eps(T)) * first(ratios) ||
            throw(ArgumentError("frequency_spacing=:log requires a geometric grid"))
        @views edges[2:end-1] .= sqrt.(frequencies[1:end-1] .* frequencies[2:end])
        edges[1] = frequencies[1] / sqrt(first(ratios))
        edges[end] = frequencies[end] * sqrt(last(ratios))
    else
        if spacing === :linear
            maximum(abs.(steps .- first(steps))) <= sqrt(eps(T)) * first(steps) ||
                throw(ArgumentError("frequency_spacing=:linear requires a uniform grid"))
        end
        edges[1] = frequencies[1] - steps[1] / 2
        @views edges[2:end-1] .= (frequencies[1:end-1] .+ frequencies[2:end]) ./ 2
        edges[end] = frequencies[end] + steps[end] / 2
    end
    edges[1] > zero(T) || throw(ArgumentError("the first bin edge must be positive"))
    return edges
end

function _validate_frequency_edges(edges::Vector{T}, centers::Vector{T}) where {T<:AbstractFloat}
    length(edges) == length(centers) + 1 || throw(DimensionMismatch(
        "frequency_edges must contain one more element than frequency centers",
    ))
    all(isfinite, edges) && all(>(zero(T)), edges) && all(>(zero(T)), diff(edges)) ||
        throw(ArgumentError("frequency edges must be finite, positive, and strictly ascending"))
    for index in eachindex(centers)
        edges[index] <= centers[index] <= edges[index + 1] || throw(ArgumentError(
            "each frequency center must lie inside its corresponding bin",
        ))
    end
    return edges
end
