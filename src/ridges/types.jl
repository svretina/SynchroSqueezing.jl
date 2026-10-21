"""Configuration for deterministic single-ridge tracking.

`continuity_scale` is expressed in Hz per time step for `:physical`
continuity and in log-Hz per time step for `:log`. `max_rate`, when present,
is an independent hard constraint expressed in Hz/s.
"""
struct RidgeConfig{T<:AbstractFloat}
    method::Symbol
    continuity::Symbol
    penalty::T
    continuity_scale::T
    max_rate::Union{Nothing,T}
    epsilon::T
    activity_threshold::T
end

function RidgeConfig(
    ::Type{T} = Float64;
    method::Symbol = :dp,
    continuity::Symbol = :physical,
    penalty::Real = 1,
    continuity_scale::Real = 1,
    max_rate::Union{Nothing,Real} = nothing,
    epsilon::Real = eps(T),
    activity_threshold::Real = 0,
) where {T<:AbstractFloat}
    method in (:maximum, :dense, :dp) || throw(ArgumentError(
        "ridge method must be :maximum, :dense, or :dp",
    ))
    continuity in (:physical, :log) || throw(ArgumentError(
        "ridge continuity must be :physical or :log",
    ))
    converted_penalty = T(penalty)
    converted_scale = T(continuity_scale)
    converted_epsilon = T(epsilon)
    converted_activity = T(activity_threshold)
    isfinite(converted_penalty) && converted_penalty >= zero(T) || throw(ArgumentError(
        "ridge penalty must be finite and nonnegative",
    ))
    isfinite(converted_scale) && converted_scale > zero(T) || throw(ArgumentError(
        "continuity_scale must be finite and positive",
    ))
    isfinite(converted_epsilon) && converted_epsilon > zero(T) || throw(ArgumentError(
        "epsilon must be finite and positive",
    ))
    isfinite(converted_activity) && converted_activity >= zero(T) || throw(ArgumentError(
        "activity_threshold must be finite and nonnegative",
    ))
    converted_rate = if max_rate === nothing
        nothing
    else
        rate = T(max_rate)
        isfinite(rate) && rate >= zero(T) || throw(ArgumentError(
            "max_rate must be finite and nonnegative",
        ))
        rate
    end
    return RidgeConfig{T}(
        method,
        continuity,
        converted_penalty,
        converted_scale,
        converted_rate,
        converted_epsilon,
        converted_activity,
    )
end

"""A constrained ridge path and its separate data-support diagnostics."""
struct RidgeResult{T<:AbstractFloat,C}
    indices::Vector{Int}
    frequencies::Vector{T}
    times::Vector{T}
    valid::BitVector
    active::BitVector
    at_frequency_boundary::BitVector
    local_energy::Vector{T}
    total_cost::T
    freqrange::Tuple{T,T}
    config::C
end

"""Configuration for second-order predictive ridge tracking.

The exact method uses `method=:exact`. `method=:pruned` additionally imposes
the heuristic fixed-Hz `prediction_bandwidth` around the constant-velocity
prediction. `frequency_scale` is in Hz per step, `curvature_scale` and
`max_chirp_rate_change` are in Hz/s, and `max_rate` is in Hz/s.
"""
struct PredictiveRidgeConfig{T<:AbstractFloat}
    method::Symbol
    frequency_penalty::T
    frequency_scale::T
    curvature_penalty::T
    curvature_scale::T
    max_rate::Union{Nothing,T}
    max_chirp_rate_change::Union{Nothing,T}
    prediction_bandwidth::Union{Nothing,T}
    initial_frequency::Union{Nothing,T}
    initial_frequency_penalty::T
    initial_frequency_scale::T
    initial_slope::Union{Nothing,T}
    initial_slope_penalty::T
    initial_slope_scale::T
    epsilon::T
    activity_threshold::T
end

function PredictiveRidgeConfig(
    ::Type{T} = Float64;
    method::Symbol = :exact,
    frequency_penalty::Real = 1,
    frequency_scale::Real = 1,
    curvature_penalty::Real = 1,
    curvature_scale::Real = 1,
    max_rate::Union{Nothing,Real} = nothing,
    max_chirp_rate_change::Union{Nothing,Real} = nothing,
    prediction_bandwidth::Union{Nothing,Real} = nothing,
    initial_frequency::Union{Nothing,Real} = nothing,
    initial_frequency_penalty::Real = 0,
    initial_frequency_scale::Real = 1,
    initial_slope::Union{Nothing,Real} = nothing,
    initial_slope_penalty::Real = 0,
    initial_slope_scale::Real = 1,
    epsilon::Real = eps(T),
    activity_threshold::Real = 0,
) where {T<:AbstractFloat}
    method in (:exact, :pruned) || throw(ArgumentError(
        "predictive ridge method must be :exact or :pruned",
    ))
    finite_nonnegative(value, name) = begin
        converted = T(value)
        isfinite(converted) && converted >= zero(T) || throw(ArgumentError(
            "$name must be finite and nonnegative",
        ))
        converted
    end
    finite_positive(value, name) = begin
        converted = T(value)
        isfinite(converted) && converted > zero(T) || throw(ArgumentError(
            "$name must be finite and positive",
        ))
        converted
    end
    optional_nonnegative(value, name) = value === nothing ? nothing :
        finite_nonnegative(value, name)
    optional_finite(value, name) = if value === nothing
        nothing
    else
        converted = T(value)
        isfinite(converted) || throw(ArgumentError("$name must be finite"))
        converted
    end
    converted_bandwidth = optional_nonnegative(
        prediction_bandwidth,
        "prediction_bandwidth",
    )
    method === :pruned && converted_bandwidth === nothing && throw(ArgumentError(
        "method=:pruned requires prediction_bandwidth in Hz",
    ))
    method === :exact && converted_bandwidth !== nothing && throw(ArgumentError(
        "prediction_bandwidth is a heuristic restriction and requires method=:pruned",
    ))
    return PredictiveRidgeConfig{T}(
        method,
        finite_nonnegative(frequency_penalty, "frequency_penalty"),
        finite_positive(frequency_scale, "frequency_scale"),
        finite_nonnegative(curvature_penalty, "curvature_penalty"),
        finite_positive(curvature_scale, "curvature_scale"),
        optional_nonnegative(max_rate, "max_rate"),
        optional_nonnegative(max_chirp_rate_change, "max_chirp_rate_change"),
        converted_bandwidth,
        optional_finite(initial_frequency, "initial_frequency"),
        finite_nonnegative(initial_frequency_penalty, "initial_frequency_penalty"),
        finite_positive(initial_frequency_scale, "initial_frequency_scale"),
        optional_finite(initial_slope, "initial_slope"),
        finite_nonnegative(initial_slope_penalty, "initial_slope_penalty"),
        finite_positive(initial_slope_scale, "initial_slope_scale"),
        finite_positive(epsilon, "epsilon"),
        finite_nonnegative(activity_threshold, "activity_threshold"),
    )
end

"""One ridge-band reconstruction with inspectable selection and ownership."""
struct ComponentReconstruction{S<:AbstractVector,T<:AbstractFloat}
    signal::S
    mask::BitMatrix
    weights::Matrix{T}
end

"""Multiple conservative ridge-band reconstructions and their SST residual."""
struct ComponentSet{S<:AbstractVector,T<:AbstractFloat}
    components::Vector{S}
    masks::Vector{BitMatrix}
    weights::Vector{Matrix{T}}
    residual::S
    residual_weights::Matrix{T}
    ownership::Symbol
end
