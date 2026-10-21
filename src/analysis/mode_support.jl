"""
    ridge_prominence(result, ridge; ridge_bandwidth_hz,
                     background_band_hz)

Compare mean squared WSST magnitude inside a ridge-centered band with the mean
in two local sidebands. `background_band_hz=(inner, outer)` is measured from
the ridge center. Both sidebands must contain bins, which independently marks
frequency-boundary cases invalid. This is local transform background, not a
detector-noise statistic.
"""
function ridge_prominence(
    result::WSSTResult,
    ridge::RidgeResult;
    ridge_bandwidth_hz::Real = 1,
    background_band_hz = (2, 6),
)
    _validate_ridge_for_result(result, ridge)
    T = eltype(result.frequencies)
    ridge_width = T(ridge_bandwidth_hz)
    inner = T(first(background_band_hz))
    outer = T(last(background_band_hz))
    isfinite(ridge_width) && ridge_width >= 0 || throw(ArgumentError(
        "ridge_bandwidth_hz must be finite and nonnegative",
    ))
    isfinite(inner) && isfinite(outer) && ridge_width < inner < outer ||
        throw(ArgumentError(
            "background_band_hz must be finite, ordered, and outside the ridge band",
        ))

    ntime = size(result.coefficients, 2)
    ratio = fill(T(NaN), ntime)
    ridge_energy = fill(T(NaN), ntime)
    background_energy = fill(T(NaN), ntime)
    valid = falses(ntime)
    for column in 1:ntime
        ridge.valid[column] || continue
        center = ridge.frequencies[column]
        ridge_rows = findall(frequency -> abs(frequency - center) <= ridge_width,
                             result.frequencies)
        lower_rows = findall(frequency -> inner <= center - frequency <= outer,
                             result.frequencies)
        upper_rows = findall(frequency -> inner <= frequency - center <= outer,
                             result.frequencies)
        isempty(ridge_rows) && continue
        (isempty(lower_rows) || isempty(upper_rows)) && continue
        ridge_values = @view result.coefficients[ridge_rows, column]
        lower_values = @view result.coefficients[lower_rows, column]
        upper_values = @view result.coefficients[upper_rows, column]
        all(isfinite, ridge_values) && all(isfinite, lower_values) &&
            all(isfinite, upper_values) || continue
        ridge_energy[column] = sum(abs2, ridge_values) / length(ridge_values)
        background_sum = sum(abs2, lower_values) + sum(abs2, upper_values)
        background_count = length(lower_values) + length(upper_values)
        background_energy[column] = background_sum / background_count
        ratio[column] = if iszero(background_energy[column])
            iszero(ridge_energy[column]) ? zero(T) : T(Inf)
        else
            ridge_energy[column] / background_energy[column]
        end
        valid[column] = true
    end
    return RidgeProminence(
        ratio, ridge_energy, background_energy, valid, ridge_width, (inner, outer),
    )
end

function _termination_reference(features, reference_amplitude, reference_mask)
    length(reference_mask) == length(features.amplitude) || throw(DimensionMismatch(
        "reference_mask must match the component length",
    ))
    if reference_amplitude !== nothing
        reference_amplitude > 0 && isfinite(reference_amplitude) || throw(ArgumentError(
            "reference_amplitude must be finite and positive",
        ))
        return eltype(features.amplitude)(reference_amplitude)
    end
    eligible = features.valid .& reference_mask .& isfinite.(features.amplitude)
    any(eligible) || throw(ArgumentError("no valid samples are available for the amplitude reference"))
    reference = maximum(features.amplitude[eligible])
    reference > 0 || throw(ArgumentError("the amplitude reference must be positive"))
    return reference
end

"""
    mode_support(features; prominence=nothing,
                 relative_amplitude_threshold=0.01,
                 prominence_threshold=nothing,
                 persistence_seconds=0, recovery_factor=1,
                 reference_amplitude=nothing, reference_mask=features.valid,
                 boundary_mask=trues(...))

Construct an operational, noise-model-independent support interval. Activity
requires the component's amplitude relative to its own reference to pass its
threshold and, when requested, local WSST prominence to pass independently.
A failed criterion must persist for `persistence_seconds` before an absorbing
termination is applied. The raw optimized ridge remains unchanged.

`fit_mask` excludes failed and boundary-unreliable samples even during a
temporary pending interval. `supported` is the stateful prefix and may include
short uncertain dips. A `recovery_factor > 1` supplies hysteresis: after a
failure, every enabled criterion must exceed its threshold times this factor
to cancel pending termination.
"""
function mode_support(
    features::ComponentFeatures;
    prominence::Union{Nothing,RidgeProminence} = nothing,
    relative_amplitude_threshold::Union{Nothing,Real} = 0.01,
    prominence_threshold::Union{Nothing,Real} = nothing,
    persistence_seconds::Real = 0,
    recovery_factor::Real = 1,
    reference_amplitude::Union{Nothing,Real} = nothing,
    reference_mask = features.valid,
    boundary_mask = trues(length(features.times)),
)
    n = length(features.times)
    length(boundary_mask) == n || throw(DimensionMismatch(
        "boundary_mask must match the component length",
    ))
    T = eltype(features.times)
    amplitude_threshold = relative_amplitude_threshold === nothing ? nothing :
        T(relative_amplitude_threshold)
    persistence = eltype(features.times)(persistence_seconds)
    recovery = eltype(features.times)(recovery_factor)
    if amplitude_threshold !== nothing
        isfinite(amplitude_threshold) && 0 < amplitude_threshold <= 1 ||
            throw(ArgumentError("relative_amplitude_threshold must lie in (0, 1]"))
    end
    isfinite(persistence) && persistence >= 0 || throw(ArgumentError(
        "persistence_seconds must be finite and nonnegative",
    ))
    isfinite(recovery) && recovery >= 1 || throw(ArgumentError(
        "recovery_factor must be finite and at least one",
    ))
    if prominence_threshold !== nothing
        prominence === nothing && throw(ArgumentError(
            "prominence data are required when prominence_threshold is set",
        ))
        isfinite(prominence_threshold) && prominence_threshold > 0 ||
            throw(ArgumentError("prominence_threshold must be finite and positive"))
    end
    amplitude_threshold === nothing && prominence_threshold === nothing &&
        throw(ArgumentError("at least one activity criterion must be enabled"))
    if prominence !== nothing
        length(prominence.ratio) == n || throw(DimensionMismatch(
            "prominence and component features must have equal length",
        ))
    end

    reference = _termination_reference(features, reference_amplitude, reference_mask)
    relative_amplitude = features.amplitude ./ reference
    amplitude_pass = amplitude_threshold === nothing ? trues(n) :
        BitVector(features.valid .& isfinite.(relative_amplitude) .&
                  (relative_amplitude .>= amplitude_threshold))
    prominence_values = prominence === nothing ? fill(T(NaN), n) : prominence.ratio
    prominence_pass = if prominence_threshold === nothing
        trues(n)
    else
        BitVector(prominence.valid .& isfinite.(prominence_values) .&
                  (prominence_values .>= T(prominence_threshold)))
    end
    base_valid = BitVector(features.valid .& boundary_mask)
    raw_activity = base_valid .& amplitude_pass .& prominence_pass
    supported = falses(n)
    uncertain = falses(n)

    start_index = findfirst(raw_activity)
    if start_index === nothing
        config = ModeSupportConfig(
            amplitude_threshold,
            prominence_threshold === nothing ? nothing : T(prominence_threshold),
            persistence,
            recovery,
            reference,
        )
        return ModeSupportResult(
            relative_amplitude, prominence_values, amplitude_pass, prominence_pass,
            raw_activity, supported, falses(n), uncertain, nothing, nothing, nothing,
            :no_initial_support, config,
        )
    end

    pending_start = nothing
    termination_index = nothing
    stop_reason = :end_of_data
    failure_reason(index) = if !boundary_mask[index]
        :boundary
    elseif !features.valid[index]
        :invalid_component
    elseif !amplitude_pass[index] && !prominence_pass[index]
        :relative_amplitude_and_insufficient_prominence
    elseif !amplitude_pass[index]
        :relative_amplitude
    else
        :insufficient_prominence
    end
    for index in start_index:n
        if pending_start === nothing
            if raw_activity[index]
                supported[index] = true
            else
                pending_start = index
                supported[index] = true
                uncertain[index] = true
                if iszero(persistence)
                    termination_index = index
                    stop_reason = failure_reason(index)
                    supported[index:end] .= false
                    break
                end
            end
        else
            amplitude_recovered = amplitude_threshold === nothing ||
                (base_valid[index] &&
                 relative_amplitude[index] >= amplitude_threshold * recovery)
            prominence_recovered = prominence_threshold === nothing ||
                (prominence.valid[index] &&
                 prominence_values[index] >= T(prominence_threshold) * recovery)
            if amplitude_recovered && prominence_recovered
                supported[pending_start:index] .= true
                uncertain[pending_start:index-1] .= true
                pending_start = nothing
            elseif features.times[index] - features.times[pending_start] >= persistence
                termination_index = pending_start
                stop_reason = failure_reason(pending_start)
                supported[pending_start:end] .= false
                uncertain[pending_start:index] .= true
                break
            else
                supported[index] = true
                uncertain[index] = true
            end
        end
    end
    fit_mask = supported .& raw_activity
    last_supported_index = findlast(fit_mask)
    last_supported_time = last_supported_index === nothing ? nothing :
        features.times[last_supported_index]
    termination_time = termination_index === nothing ? nothing :
        features.times[termination_index]
    config = ModeSupportConfig(
        amplitude_threshold,
        prominence_threshold === nothing ? nothing : T(prominence_threshold),
        persistence,
        recovery,
        reference,
    )
    return ModeSupportResult(
        relative_amplitude, prominence_values, amplitude_pass, prominence_pass,
        raw_activity, supported, fit_mask, uncertain, last_supported_index,
        last_supported_time, termination_time, stop_reason, config,
    )
end

"""
    supported_reconstruction(component_or_signal, support)

Preserve a raw `ComponentReconstruction` or component vector and return a
separate signal that is zero outside the operationally supported prefix. No
raw coefficients, mask, weights, or ridge trajectory are mutated.
"""
function supported_reconstruction(
    component::ComponentReconstruction,
    support::ModeSupportResult,
)
    length(component.signal) == length(support.supported) || throw(DimensionMismatch(
        "component and support result must have equal length",
    ))
    signal = ifelse.(support.supported, component.signal, zero(eltype(component.signal)))
    return SupportedComponentReconstruction(component, signal, support)
end

function supported_reconstruction(
    signal::AbstractVector,
    support::ModeSupportResult,
)
    length(signal) == length(support.supported) || throw(DimensionMismatch(
        "signal and support result must have equal length",
    ))
    supported_signal = ifelse.(support.supported, signal, zero(eltype(signal)))
    return SupportedComponentReconstruction(signal, supported_signal, support)
end

"""
    mode_support_sensitivity(features; relative_amplitude_thresholds, kwargs...)

Repeat operational support and downstream exponential damping fits over
several amplitude thresholds. The spread reports decision sensitivity rather
than conditional regression uncertainty.
"""
function mode_support_sensitivity(
    features::ComponentFeatures;
    relative_amplitude_thresholds,
    fit_kwargs = NamedTuple(),
    support_kwargs...,
)
    thresholds = collect(relative_amplitude_thresholds)
    isempty(thresholds) && throw(ArgumentError("at least one amplitude threshold is required"))
    T = eltype(features.times)
    cases = ModeSupportSensitivityCase[]
    for threshold in thresholds
        support = mode_support(
            features;
            relative_amplitude_threshold = threshold,
            support_kwargs...,
        )
        fit = fit_exponential_damping(
            features;
            valid_mask = support.fit_mask,
            fit_kwargs...,
        )
        push!(cases, ModeSupportSensitivityCase(T(threshold), support, fit))
    end
    damping_times = [case.damping_fit.damping_time for case in cases
                     if isfinite(case.damping_fit.damping_time)]
    isempty(damping_times) && return ModeSupportSensitivity(
        cases, T(NaN), T(NaN), T(NaN),
    )
    minimum_time, maximum_time = extrema(damping_times)
    return ModeSupportSensitivity(
        cases, minimum_time, maximum_time, maximum_time - minimum_time,
    )
end
