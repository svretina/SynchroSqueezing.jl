function _detect_frequency_spacing(frequencies::Vector{T}) where {T<:AbstractFloat}
    steps = diff(frequencies)
    linear_error = maximum(abs.(steps .- first(steps))) / first(steps)
    ratios = frequencies[2:end] ./ frequencies[1:end-1]
    log_error = maximum(abs.(ratios .- first(ratios))) / first(ratios)
    tolerance = sqrt(eps(T))
    linear_error <= tolerance && return :linear
    log_error <= tolerance && return :log
    return :arbitrary
end

function _resolve_frequency_grid(
    cwt_result::CWTResult{A,D,T};
    frequency_grid,
    frequency_edges,
    freqrange,
    nfrequencies,
    frequency_spacing,
) where {A,D,T<:AbstractFloat}
    frequency_grid === nothing || freqrange === nothing || throw(ArgumentError(
        "provide either frequency_grid or freqrange, not both",
    ))
    nfrequencies > 1 || throw(ArgumentError("nfrequencies must exceed one"))
    spacing = frequency_spacing
    centers = if frequency_grid !== nothing
        T.(collect(frequency_grid))
    elseif freqrange === nothing
        nfrequencies == length(cwt_result.frequencies) || throw(ArgumentError(
            "nfrequencies only changes a generated grid when freqrange is provided",
        ))
        T.(cwt_result.frequencies)
    else
        selected_range = freqrange
        length(selected_range) == 2 || throw(ArgumentError("freqrange must contain (fmin, fmax) centers"))
        spacing === :auto && (spacing = :linear)
        if spacing === :linear
            linear_frequency_grid(selected_range[1], selected_range[2], nfrequencies; precision = T)
        elseif spacing === :log
            log_frequency_grid(selected_range[1], selected_range[2], nfrequencies; precision = T)
        else
            throw(ArgumentError("generated grids require frequency_spacing=:linear or :log"))
        end
    end
    length(centers) > 1 || throw(ArgumentError("at least two output frequencies are required"))
    all(isfinite, centers) && all(>(zero(T)), centers) && all(>(zero(T)), diff(centers)) ||
        throw(ArgumentError("frequency centers must be finite, positive, and strictly ascending"))
    centers[end] < cwt_result.metadata.sample_rate / 2 || throw(ArgumentError(
        "the highest frequency center must be below Nyquist",
    ))
    spacing === :auto && (spacing = _detect_frequency_spacing(centers))
    spacing in (:linear, :log, :arbitrary) || throw(ArgumentError(
        "frequency_spacing must be :auto, :linear, :log, or :arbitrary",
    ))
    edges = frequency_edges === nothing ?
        _frequency_bin_edges(centers; spacing) :
        _validate_frequency_edges(T.(collect(frequency_edges)), centers)
    return centers, edges, spacing
end

"""
    wsst(cwt_result; kwargs...)

Compute a positive-frequency first-order WSST from a validated CWT. The output
stores backend-derived, quadrature-weighted complex contributions integrated
per frequency bin. `freqrange` denotes the first and last bin centers. Set
`implementation=:reference` to use the independent two-stage phase path.
"""
function wsst(
    cwt_result::CWTResult{A,D,T};
    frequency_grid = nothing,
    frequency_edges = nothing,
    freqrange = nothing,
    nfrequencies::Integer = length(cwt_result.frequencies),
    frequency_spacing::Symbol = :auto,
    reassignment::Symbol = :nearest,
    threshold::AbstractMagnitudeThreshold = RelativeThreshold(T(1e-6)),
    preserve_residual::Bool = true,
    keep_cwt::Bool = false,
    keep_derivative::Bool = false,
    keep_phase::Bool = false,
    implementation::Symbol = :optimized,
) where {A,D,T<:AbstractFloat}
    cwt_result.metadata.signal_kind === :two_sided_complex && throw(ArgumentError(
        "positive-frequency WSST supports real or positive-frequency analytic complex signals only",
    ))
    implementation in (:reference, :optimized) || throw(ArgumentError(
        "implementation must be :reference or :optimized",
    ))
    centers, edges, spacing = _resolve_frequency_grid(
        cwt_result;
        frequency_grid,
        frequency_edges,
        freqrange,
        nfrequencies,
        frequency_spacing,
    )
    C = eltype(cwt_result.coefficients)
    coefficients = zeros(C, length(centers), size(cwt_result.coefficients, 2))
    phase = if implementation === :reference || keep_phase
        instantaneous_frequency(cwt_result; threshold)
    else
        nothing
    end
    diagnostics, residual = if implementation === :reference
        synchrosqueeze_reference!(
            coefficients,
            cwt_result,
            phase,
            centers,
            edges,
            threshold,
            reassignment;
            preserve_residual,
        )
    else
        synchrosqueeze!(
            coefficients,
            cwt_result,
            centers,
            edges,
            threshold,
            reassignment;
            preserve_residual,
        )
    end
    reconstruction = ReconstructionConvention(
        copy(cwt_result.squeezing_weights),
        cwt_result.admissibility,
        copy(cwt_result.synthesis_response),
        eps(T) * maximum(abs, cwt_result.synthesis_response),
        cwt_result.reliable_frequency_range,
        cwt_result.standard_compatible ? :standard_analytic_l1 : :diagonal_wavelet_frame,
        cwt_result.metadata.normalization,
        cwt_result.metadata.signal_kind,
    )
    return WSSTResult(
        coefficients,
        centers,
        edges,
        spacing,
        copy(cwt_result.scales),
        reassignment,
        threshold,
        diagnostics,
        reconstruction,
        residual,
        cwt_result.dc_component,
        cwt_result.metadata,
        keep_cwt ? cwt_result.coefficients : nothing,
        keep_derivative ? cwt_result.derivative : nothing,
        keep_phase ? phase.values : nothing,
        keep_phase ? phase.valid : nothing,
    )
end

"""Compute a CWT with an existing plan and immediately synchrosqueeze it."""
function wsst(plan::CWTPlan, signal::AbstractVector; kwargs...)
    return wsst(cwt(plan, signal); kwargs...)
end

"""
    wsst(signal; fs, cwt_options=(;), kwargs...)

Plan, compute, and synchrosqueeze one signal. Backend planning options belong
in the `cwt_options` named tuple; remaining keywords configure WSST.
"""
function wsst(signal::AbstractVector; fs::Real, cwt_options::NamedTuple = (;), kwargs...)
    (:fs in keys(cwt_options) || :precision in keys(cwt_options)) && throw(ArgumentError(
        "cwt_options must not contain fs or precision",
    ))
    precision = _realtype(eltype(signal))
    plan = plan_cwt(length(signal); fs, precision, cwt_options...)
    return wsst(plan, signal; kwargs...)
end
