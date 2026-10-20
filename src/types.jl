"""A backend capable of producing complex continuous wavelet coefficients."""
abstract type AbstractCWTBackend end

"""Marker type for the `ContinuousWavelets.jl` CWT backend."""
struct ContinuousWaveletsBackend <: AbstractCWTBackend end

"""Analytic-support diagnostics for the rows of a CWT."""
struct WaveletValidation{T<:AbstractFloat}
    dc_leakage::Vector{T}
    valid_rows::BitVector
    tolerance::T
end

"""Immutable, reusable CWT configuration for one signal length and precision."""
struct CWTPlan{B<:AbstractCWTBackend,C,D,S,F,I,T<:AbstractFloat}
    backend::B
    backend_config::C
    daughters::D
    scales::S
    frequencies::F
    frequency_weights::Vector{T}
    inverse_style::I
    validation::WaveletValidation{T}
    sample_rate::T
    signal_length::Int
end

"""Metadata needed to interpret and invert a CWT result."""
struct CWTMetadata{T<:AbstractFloat,C,I}
    sample_rate::T
    sample_interval::T
    signal_length::Int
    padded_length::Int
    crop::UnitRange{Int}
    boundary::Symbol
    normalization::Symbol
    inverse_style::I
    backend_config::C
end

"""
    CWTResult

Complex CWT coefficients with scale/frequency on the first axis and time on the
second. The coefficient and derivative arrays are non-copying transposed views
of the backend's native time × scale matrices.
"""
struct CWTResult{A<:AbstractMatrix,D<:AbstractMatrix,T<:AbstractFloat,M}
    coefficients::A
    derivative::D
    frequencies::Vector{T}
    scales::Vector{T}
    frequency_weights::Vector{T}
    valid_rows::BitVector
    reliability::BitMatrix
    metadata::M
end

"""Instantaneous-frequency estimates and their coefficient-level validity mask."""
struct FrequencyEstimate{T<:AbstractFloat}
    values::Matrix{T}
    valid::BitMatrix
end
