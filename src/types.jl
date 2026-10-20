"""A backend capable of producing complex continuous wavelet coefficients."""
abstract type AbstractCWTBackend end

"""A policy for accepting CWT coefficients during phase reassignment."""
abstract type AbstractMagnitudeThreshold end

"""Apply no positive magnitude cutoff; exact zeros remain phase-invalid."""
struct NoThreshold <: AbstractMagnitudeThreshold end

"""Accept coefficients whose magnitude is strictly greater than `value`."""
struct AbsoluteThreshold{T<:AbstractFloat} <: AbstractMagnitudeThreshold
    value::T
    function AbsoluteThreshold(value::T) where {T<:AbstractFloat}
        isfinite(value) || value == T(Inf) || throw(ArgumentError(
            "absolute threshold must not be NaN",
        ))
        value >= zero(T) || throw(ArgumentError("absolute threshold must be nonnegative"))
        return new{T}(value)
    end
end

"""Accept coefficients above a fraction of the global finite CWT maximum."""
struct RelativeThreshold{T<:AbstractFloat} <: AbstractMagnitudeThreshold
    fraction::T
    reference::Symbol
    function RelativeThreshold{T}(fraction::T, reference::Symbol) where {T<:AbstractFloat}
        isfinite(fraction) && fraction >= zero(T) || throw(ArgumentError(
            "relative threshold fraction must be finite and nonnegative",
        ))
        reference === :global || throw(ArgumentError(
            "only reference=:global is currently supported for relative thresholds",
        ))
        return new{T}(fraction, reference)
    end
end

RelativeThreshold(fraction::T; reference::Symbol = :global) where {T<:AbstractFloat} =
    RelativeThreshold{T}(fraction, reference)

"""Marker type for the `ContinuousWavelets.jl` CWT backend."""
struct ContinuousWaveletsBackend <: AbstractCWTBackend end

"""Analytic-support diagnostics for the rows of a CWT."""
struct WaveletValidation{T<:AbstractFloat}
    dc_leakage::Vector{T}
    valid_rows::BitVector
    tolerance::T
end

"""Immutable, reusable CWT configuration for one signal length and precision."""
struct CWTPlan{B<:AbstractCWTBackend,C,D,S,F,I,R,T<:AbstractFloat}
    backend::B
    backend_config::C
    daughters::D
    scales::S
    frequencies::F
    frequency_weights::Vector{T}
    inverse_style::I
    reconstruction_weights::R
    squeezing_weights::Vector{T}
    synthesis_response::Vector{T}
    admissibility::Complex{T}
    standard_compatible::Bool
    reliable_frequency_range::Tuple{T,T}
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
    signal_kind::Symbol
    inverse_style::I
    backend_config::C
end

"""
    CWTResult

Complex CWT coefficients with scale/frequency on the first axis and time on the
second. The coefficient and derivative arrays are non-copying transposed views
of the backend's native time × scale matrices.
"""
struct CWTResult{A<:AbstractMatrix,D<:AbstractMatrix,T<:AbstractFloat,M,S,R}
    coefficients::A
    derivative::D
    frequencies::Vector{T}
    scales::Vector{T}
    frequency_weights::Vector{T}
    reconstruction_weights::R
    squeezing_weights::Vector{T}
    synthesis_response::Vector{T}
    admissibility::Complex{T}
    standard_compatible::Bool
    reliable_frequency_range::Tuple{T,T}
    valid_rows::BitVector
    reliability::BitMatrix
    metadata::M
    dc_component::S
end

"""Instantaneous-frequency estimates and their coefficient-level validity mask."""
struct FrequencyEstimate{T<:AbstractFloat}
    values::Matrix{T}
    valid::BitMatrix
end

"""Backend-derived weights and synthesis rule used by inverse WSST."""
struct ReconstructionConvention{T<:AbstractFloat}
    scale_weights::Vector{T}
    admissibility::Complex{T}
    synthesis_response::Vector{T}
    response_floor::T
    reliable_frequency_range::Tuple{T,T}
    method::Symbol
    normalization::Symbol
    signal_kind::Symbol
end

"""Per-column accounting for a reference synchrosqueezing calculation."""
struct WSSTDiagnostics{T<:AbstractFloat,C<:Complex{T}}
    accepted::Vector{Int}
    rejected_magnitude::Vector{Int}
    rejected_unreliable::Vector{Int}
    rejected_nonfinite_input::Vector{Int}
    rejected_unstable_phase::Vector{Int}
    rejected_out_of_range::Vector{Int}
    unrecoverable::Vector{Int}
    accepted_sum::Vector{C}
    deposited_sum::Vector{C}
end

"""
    WSSTResult

Bin-integrated complex synchrosqueezed coefficients with frequency on the
first axis and time on the second. CWT, derivative, and phase arrays are
retained only when requested.
"""
struct WSSTResult{T<:AbstractFloat,C<:Complex{T},M,P,R,DC,W,D,PH,PV}
    coefficients::Matrix{C}
    frequencies::Vector{T}
    frequency_edges::Vector{T}
    frequency_spacing::Symbol
    scales::Vector{T}
    reassignment::Symbol
    threshold::P
    diagnostics::WSSTDiagnostics{T,C}
    reconstruction::ReconstructionConvention{T}
    residual::R
    dc_component::DC
    metadata::M
    cwt::W
    derivative::D
    phase::PH
    phase_valid::PV
end
