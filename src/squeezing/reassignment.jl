"""Return zero outside the grid, otherwise the nearest physical center index."""
function _nearest_bin(frequency::T, centers::Vector{T}, edges::Vector{T}) where {T<:AbstractFloat}
    edges[1] <= frequency <= edges[end] || return 0
    upper = searchsortedfirst(centers, frequency)
    upper == 1 && return 1
    upper > length(centers) && return length(centers)
    return frequency - centers[upper - 1] <= centers[upper] - frequency ? upper - 1 : upper
end

"""Deposit one contribution into neighboring physical centers."""
function _deposit_linear!(output, frequency::T, value, centers::Vector{T}, edges::Vector{T}, column::Int) where {T<:AbstractFloat}
    edges[1] <= frequency <= edges[end] || return false
    upper = searchsortedfirst(centers, frequency)
    if upper == 1
        output[1, column] += value
    elseif upper > length(centers)
        output[end, column] += value
    else
        lower = upper - 1
        fraction = (frequency - centers[lower]) / (centers[upper] - centers[lower])
        output[lower, column] += (one(T) - fraction) * value
        output[upper, column] += fraction * value
    end
    return true
end

function _accounting(::Type{C}, ntime::Int) where {T<:AbstractFloat,C<:Complex{T}}
    return (
        accepted = zeros(Int, ntime),
        rejected_magnitude = zeros(Int, ntime),
        rejected_unreliable = zeros(Int, ntime),
        rejected_nonfinite_input = zeros(Int, ntime),
        rejected_unstable_phase = zeros(Int, ntime),
        rejected_out_of_range = zeros(Int, ntime),
        unrecoverable = zeros(Int, ntime),
        accepted_sum = zeros(C, ntime),
        deposited_sum = zeros(C, ntime),
    )
end

function _diagnostics(accounting)
    return WSSTDiagnostics(
        accounting.accepted,
        accounting.rejected_magnitude,
        accounting.rejected_unreliable,
        accounting.rejected_nonfinite_input,
        accounting.rejected_unstable_phase,
        accounting.rejected_out_of_range,
        accounting.unrecoverable,
        accounting.accepted_sum,
        accounting.deposited_sum,
    )
end

function _deposit!(output, reassignment, frequency, contribution, centers, edges, column)
    if reassignment === :nearest
        output[_nearest_bin(frequency, centers, edges), column] += contribution
    elseif reassignment === :linear
        _deposit_linear!(output, frequency, contribution, centers, edges, column)
    else
        throw(ArgumentError("reassignment must be :nearest or :linear"))
    end
    return nothing
end

function _add_residual!(residual, contribution, column)
    residual === nothing || (residual[column] += contribution)
    return nothing
end

"""
    synchrosqueeze_reference!(output, cwt_result, phase, centers, edges, threshold,
                              reassignment; preserve_residual=true)

Readable two-stage reference reassignment. Frequencies are computed separately
in `phase`; deposited values use positive log-frequency quadrature and are
stored as bin-integrated complex contributions.
"""
function synchrosqueeze_reference!(
    output::Matrix{C},
    cwt_result::CWTResult,
    phase::FrequencyEstimate{T},
    centers::Vector{T},
    edges::Vector{T},
    threshold::AbstractMagnitudeThreshold,
    reassignment::Symbol;
    preserve_residual::Bool = true,
) where {T<:AbstractFloat,C<:Complex{T}}
    size(phase.values) == size(cwt_result.coefficients) || throw(DimensionMismatch(
        "phase estimates must match CWT dimensions",
    ))
    weights = cwt_result.squeezing_weights
    fill!(output, zero(C))
    ntime = size(output, 2)
    accounting = _accounting(C, ntime)
    residual = preserve_residual ? zeros(C, ntime) : nothing
    magnitude_threshold = T(_threshold_magnitude(threshold, cwt_result.coefficients))

    for column in axes(cwt_result.coefficients, 2), row in axes(cwt_result.coefficients, 1)
        coefficient = cwt_result.coefficients[row, column]
        derivative = cwt_result.derivative[row, column]
        finite_coefficient = isfinite(coefficient)
        contribution = finite_coefficient ? weights[row] * coefficient : zero(C)

        if !cwt_result.reliability[row, column]
            accounting.rejected_unreliable[column] += 1
            finite_coefficient ? _add_residual!(residual, contribution, column) :
                                 (accounting.unrecoverable[column] += 1)
        elseif !finite_coefficient || !isfinite(derivative)
            accounting.rejected_nonfinite_input[column] += 1
            finite_coefficient ? _add_residual!(residual, contribution, column) :
                                 (accounting.unrecoverable[column] += 1)
        elseif abs(coefficient) <= magnitude_threshold
            accounting.rejected_magnitude[column] += 1
            _add_residual!(residual, contribution, column)
        elseif !phase.valid[row, column]
            accounting.rejected_unstable_phase[column] += 1
            _add_residual!(residual, contribution, column)
        else
            frequency = phase.values[row, column]
            if frequency < edges[1] || frequency > edges[end]
                accounting.rejected_out_of_range[column] += 1
                _add_residual!(residual, contribution, column)
            else
                accounting.accepted[column] += 1
                accounting.accepted_sum[column] += contribution
                _deposit!(output, reassignment, frequency, contribution, centers, edges, column)
                accounting.deposited_sum[column] += contribution
            end
        end
    end
    return _diagnostics(accounting), residual
end

"""
    synchrosqueeze!(output, cwt_result, centers, edges, threshold, reassignment;
                    preserve_residual=true)

Fused phase-and-deposition kernel. This is allocation-free per coefficient and
is regression-tested against `synchrosqueeze_reference!`.
"""
function synchrosqueeze!(
    output::Matrix{C},
    cwt_result::CWTResult,
    centers::Vector{T},
    edges::Vector{T},
    threshold::AbstractMagnitudeThreshold,
    reassignment::Symbol;
    preserve_residual::Bool = true,
) where {T<:AbstractFloat,C<:Complex{T}}
    weights = cwt_result.squeezing_weights
    fill!(output, zero(C))
    ntime = size(output, 2)
    accounting = _accounting(C, ntime)
    residual = preserve_residual ? zeros(C, ntime) : nothing
    magnitude_threshold = T(_threshold_magnitude(threshold, cwt_result.coefficients))

    for column in axes(cwt_result.coefficients, 2), row in axes(cwt_result.coefficients, 1)
        coefficient = cwt_result.coefficients[row, column]
        derivative = cwt_result.derivative[row, column]
        finite_coefficient = isfinite(coefficient)
        contribution = finite_coefficient ? weights[row] * coefficient : zero(C)

        if !cwt_result.reliability[row, column]
            accounting.rejected_unreliable[column] += 1
            finite_coefficient ? _add_residual!(residual, contribution, column) :
                                 (accounting.unrecoverable[column] += 1)
        elseif !finite_coefficient || !isfinite(derivative)
            accounting.rejected_nonfinite_input[column] += 1
            finite_coefficient ? _add_residual!(residual, contribution, column) :
                                 (accounting.unrecoverable[column] += 1)
        elseif abs(coefficient) <= magnitude_threshold
            accounting.rejected_magnitude[column] += 1
            _add_residual!(residual, contribution, column)
        else
            frequency = _phase_frequency(coefficient, derivative)
            if !isfinite(frequency)
                accounting.rejected_unstable_phase[column] += 1
                _add_residual!(residual, contribution, column)
            elseif frequency < edges[1] || frequency > edges[end]
                accounting.rejected_out_of_range[column] += 1
                _add_residual!(residual, contribution, column)
            else
                accounting.accepted[column] += 1
                accounting.accepted_sum[column] += contribution
                _deposit!(output, reassignment, frequency, contribution, centers, edges, column)
                accounting.deposited_sum[column] += contribution
            end
        end
    end
    return _diagnostics(accounting), residual
end
