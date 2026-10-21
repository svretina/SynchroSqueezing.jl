"""Envelope, unwrapped phase, and optional phase-derived frequency of a component."""
struct ComponentFeatures{T<:AbstractFloat}
    times::Vector{T}
    amplitude::Vector{T}
    phase::Vector{T}
    instantaneous_frequency::Union{Nothing,Vector{T}}
    valid::BitVector
end

"""Result and diagnostics from a constant exponential damping fit."""
struct DampingFit{T<:AbstractFloat}
    damping_time::T
    decay_rate::T
    reference_amplitude::T
    reference_time::T
    fitted_amplitude::Vector{T}
    residuals::Vector{T}
    log_residuals::Vector{T}
    selected::BitVector
    fit_window::Tuple{T,T}
    method::Symbol
    status::Symbol
    rmse::T
    log_rmse::T
    r_squared::T
    decay_rate_standard_error::T
    damping_time_standard_error::T
    iterations::Int
    converged::Bool
end

"""One reconstruction/ridge/window case in a damping sensitivity study."""
struct DampingSensitivityCase{T<:AbstractFloat,F}
    ridge_index::Int
    bandwidth_hz::T
    fit_window::Tuple{T,T}
    fit::F
end

"""Between-case damping sensitivity, kept separate from fit uncertainty."""
struct DampingSensitivity{T<:AbstractFloat,C}
    cases::Vector{C}
    minimum_damping_time::T
    maximum_damping_time::T
    damping_time_spread::T
end

"""Local WSST ridge-to-sideband energy prominence and its validity metadata."""
struct RidgeProminence{T<:AbstractFloat}
    ratio::Vector{T}
    ridge_energy::Vector{T}
    background_energy::Vector{T}
    valid::BitVector
    ridge_bandwidth_hz::T
    background_band_hz::Tuple{T,T}
end

"""Configuration recorded with an operational mode-termination result."""
struct ModeSupportConfig{T<:AbstractFloat}
    relative_amplitude_threshold::Union{Nothing,T}
    prominence_threshold::Union{Nothing,T}
    persistence_seconds::T
    recovery_factor::T
    reference_amplitude::T
end

"""Supported prefix, fit mask, diagnostics, and absorbing termination time."""
struct ModeSupportResult{T<:AbstractFloat,C}
    relative_amplitude::Vector{T}
    prominence::Vector{T}
    amplitude_pass::BitVector
    prominence_pass::BitVector
    raw_activity::BitVector
    supported::BitVector
    fit_mask::BitVector
    uncertain::BitVector
    last_supported_index::Union{Nothing,Int}
    last_supported_time::Union{Nothing,T}
    stop_time::Union{Nothing,T}
    stop_reason::Symbol
    config::C
end

"""A raw component reconstruction and its separately masked usable signal."""
struct SupportedComponentReconstruction{R,S,M}
    raw::R
    signal::S
    support::M
end

"""One amplitude-threshold case with its support decision and damping fit."""
struct ModeSupportSensitivityCase{T<:AbstractFloat,S,F}
    relative_amplitude_threshold::T
    support::S
    damping_fit::F
end

"""Damping variation across operational relative-amplitude thresholds."""
struct ModeSupportSensitivity{T<:AbstractFloat,C}
    cases::Vector{C}
    minimum_damping_time::T
    maximum_damping_time::T
    damping_time_spread::T
end
