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

function _reconstruction_convention(cwt_result::CWTResult, ::Type{T}) where {T<:AbstractFloat}
    return ReconstructionConvention(
        copy(cwt_result.squeezing_weights),
        cwt_result.admissibility,
        copy(cwt_result.synthesis_response),
        eps(T) * maximum(abs, cwt_result.synthesis_response),
        cwt_result.reliable_frequency_range,
        cwt_result.standard_compatible ? :standard_analytic_l1 : :diagonal_wavelet_frame,
        cwt_result.metadata.normalization,
        cwt_result.metadata.signal_kind,
    )
end

function _first_order_metadata(cwt_result::CWTResult, implementation::Symbol, ::Type{T}) where {T}
    return ReassignmentMetadata(
        1,
        :phase_derivative,
        implementation,
        T(NaN),
        :none,
        T(Inf),
        zero(T),
        zero(T),
        zero(T),
        cwt_result.metadata.backend_config,
    )
end

function _second_order_metadata(
    cwt_result::CWTResult,
    estimate::SecondOrderEstimate{T},
    implementation::Symbol,
    denominator_threshold,
    fallback,
    max_abs_chirp_rate,
) where {T}
    candidates = count(isfinite, estimate.first_order)
    denominator = max(candidates, 1)
    second_order_valid = count(estimate.second_order_valid)
    fallback_count = count(estimate.used_fallback)
    return ReassignmentMetadata(
        2,
        :time_derivative_wsst2,
        implementation,
        T(denominator_threshold),
        fallback,
        T(max_abs_chirp_rate),
        T(second_order_valid / denominator),
        T(fallback_count / denominator),
        T((candidates - second_order_valid - fallback_count) / denominator),
        cwt_result.metadata.backend_config,
    )
end

function _wsst_from_phase(
    cwt_result::CWTResult{A,D,T},
    phase,
    reassignment_metadata;
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
    implementation::Symbol,
) where {A,D,T<:AbstractFloat}
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
    return WSSTResult(
        coefficients,
        centers,
        edges,
        spacing,
        copy(cwt_result.scales),
        reassignment,
        threshold,
        reassignment_metadata,
        diagnostics,
        _reconstruction_convention(cwt_result, T),
        residual,
        cwt_result.dc_component,
        cwt_result.metadata,
        keep_cwt ? cwt_result.coefficients : nothing,
        keep_derivative ? cwt_result.derivative : nothing,
        keep_phase ? phase.values : nothing,
        keep_phase ? phase.valid : nothing,
    )
end

"""
    wsst(cwt_result; order=1, kwargs...)

Compute a positive-frequency first-order WSST from a validated CWT. For
`order=2`, call `wsst(plan, signal; order=2, ...)` so the required auxiliary
wavelet transforms can be computed without changing first-order CWT cost.
"""
function wsst(
    cwt_result::CWTResult{A,D,T};
    order::Integer = 1,
    implementation::Symbol = :auto,
    threshold::AbstractMagnitudeThreshold = RelativeThreshold(T(1e-6)),
    kwargs...,
) where {A,D,T<:AbstractFloat}
    order == 1 || throw(ArgumentError(
        "order=2 requires wsst(plan, signal; order=2, ...) or a SecondOrderOperators object",
    ))
    cwt_result.metadata.signal_kind === :two_sided_complex && throw(ArgumentError(
        "positive-frequency WSST supports real or positive-frequency analytic complex signals only",
    ))
    resolved = implementation === :auto ? :optimized : implementation
    resolved in (:reference, :optimized) || throw(ArgumentError(
        "implementation must be :auto, :reference, or :optimized",
    ))
    phase = resolved === :reference || get(kwargs, :keep_phase, false) ?
        instantaneous_frequency(cwt_result; threshold) : nothing
    metadata = _first_order_metadata(cwt_result, resolved, T)
    return _wsst_from_phase(
        cwt_result,
        phase,
        metadata;
        threshold,
        implementation = resolved,
        kwargs...,
    )
end

"""Compute second-order WSST from an explicitly prepared operator bundle."""
function wsst(
    operators::SecondOrderOperators;
    order::Integer = 2,
    implementation::Symbol = :auto,
    threshold = nothing,
    denominator_threshold = nothing,
    fallback::Symbol = :first_order,
    max_abs_chirp_rate::Real = Inf,
    kwargs...,
)
    order == 2 || throw(ArgumentError("SecondOrderOperators require order=2"))
    operators.cwt.metadata.signal_kind === :two_sided_complex && throw(ArgumentError(
        "positive-frequency WSST supports real or positive-frequency analytic complex signals only",
    ))
    T = eltype(operators.cwt.frequencies)
    threshold === nothing && (threshold = RelativeThreshold(T(1e-6)))
    threshold isa AbstractMagnitudeThreshold || throw(ArgumentError(
        "threshold must be an AbstractMagnitudeThreshold",
    ))
    denominator_threshold === nothing && (denominator_threshold = sqrt(eps(T)))
    resolved = implementation === :auto ? :reference : implementation
    resolved === :reference || throw(ArgumentError(
        "second-order WSST currently provides only implementation=:reference",
    ))
    estimate = second_order_frequency(
        operators;
        threshold,
        denominator_threshold,
        fallback,
        max_abs_chirp_rate,
    )
    phase = FrequencyEstimate(estimate.values, estimate.valid)
    metadata = _second_order_metadata(
        operators.cwt,
        estimate,
        resolved,
        denominator_threshold,
        fallback,
        max_abs_chirp_rate,
    )
    return _wsst_from_phase(
        operators.cwt,
        phase,
        metadata;
        threshold,
        implementation = resolved,
        kwargs...,
    )
end

"""Compute a CWT with an existing plan and immediately synchrosqueeze it."""
function wsst(plan::CWTPlan, signal::AbstractVector; order::Integer = 1, kwargs...)
    order == 1 && return wsst(cwt(plan, signal); order, kwargs...)
    order == 2 && return wsst(second_order_operators(plan, signal); order, kwargs...)
    throw(ArgumentError("WSST order must be 1 or 2"))
end

"""
    wsst(signal; fs, cwt_options=(;), kwargs...)

Plan, compute, and synchrosqueeze one signal. Backend planning options belong
in the `cwt_options` named tuple; remaining keywords configure WSST.
"""
function wsst(
    signal::AbstractVector;
    fs::Real,
    cwt_options::NamedTuple = (;),
    order::Integer = 1,
    kwargs...,
)
    (:fs in keys(cwt_options) || :precision in keys(cwt_options)) && throw(ArgumentError(
        "cwt_options must not contain fs or precision",
    ))
    precision = _realtype(eltype(signal))
    plan = plan_cwt(length(signal); fs, precision, cwt_options...)
    return wsst(plan, signal; order, kwargs...)
end
