function _predictive_frequency_cost(previous, current, config::PredictiveRidgeConfig)
    return config.frequency_penalty * abs2(
        (current - previous) / config.frequency_scale,
    )
end

_physical_slope(previous, current, delta_time) = (current - previous) / delta_time

function _predictive_curvature(previous, current, following, previous_step, next_step,
                               config::PredictiveRidgeConfig)
    previous_slope = _physical_slope(previous, current, previous_step)
    next_slope = _physical_slope(current, following, next_step)
    slope_change = next_slope - previous_slope
    cost = config.curvature_penalty * abs2(slope_change / config.curvature_scale)
    return cost, slope_change, previous_slope
end

function _predictive_pair_allowed(previous, current, delta_time,
                                  config::PredictiveRidgeConfig)
    config.max_rate === nothing && return true
    return abs(current - previous) <= config.max_rate * delta_time
end

function _predictive_transition_allowed(
    previous,
    current,
    following,
    previous_step,
    next_step,
    config::PredictiveRidgeConfig,
)
    _predictive_pair_allowed(current, following, next_step, config) || return false
    _, slope_change, previous_slope = _predictive_curvature(
        previous,
        current,
        following,
        previous_step,
        next_step,
        config,
    )
    if config.max_chirp_rate_change !== nothing &&
       abs(slope_change) > config.max_chirp_rate_change
        return false
    end
    if config.method === :pruned
        prediction = current + previous_slope * next_step
        abs(following - prediction) <= config.prediction_bandwidth || return false
    end
    return true
end

function _predictive_initial_cost(first_frequency, initial_slope,
                                  config::PredictiveRidgeConfig)
    T = typeof(config.epsilon)
    cost = zero(T)
    if config.initial_frequency !== nothing
        cost += config.initial_frequency_penalty * abs2(
            (first_frequency - config.initial_frequency) / config.initial_frequency_scale,
        )
    end
    if config.initial_slope !== nothing
        cost += config.initial_slope_penalty * abs2(
            (initial_slope - config.initial_slope) / config.initial_slope_scale,
        )
    end
    return cost
end

function _predictive_fallback(
    coefficients,
    frequencies,
    times,
    freqrange,
    config;
    valid_mask,
    frequency_envelope,
)
    initial_prior = if config.initial_frequency === nothing
        nothing
    else
        config.initial_frequency_penalty .* abs2.(
            (frequencies .- config.initial_frequency) ./ config.initial_frequency_scale,
        )
    end
    return extract_ridge(
        coefficients,
        frequencies;
        times,
        freqrange,
        method = :dense,
        penalty = config.frequency_penalty,
        continuity_scale = config.frequency_scale,
        max_rate = config.max_rate,
        epsilon = config.epsilon,
        activity_threshold = config.activity_threshold,
        valid_mask,
        frequency_envelope,
        initial_prior,
    )
end

function _second_order_dynamic_ridge(
    coefficients,
    frequencies,
    times,
    freqrange,
    config;
    valid_mask,
    frequency_envelope,
)
    T = eltype(frequencies)
    ntime = length(times)
    ntime <= 2 && return _predictive_fallback(
        coefficients,
        frequencies,
        times,
        freqrange,
        config;
        valid_mask,
        frequency_envelope,
    )
    eligible = _eligible_states(
        coefficients,
        frequencies,
        times,
        freqrange,
        valid_mask,
        frequency_envelope,
    )
    costs, energies, global_reference = _observation_costs(coefficients, eligible, config)
    allowed_rows = findall(frequency -> freqrange[1] <= frequency <= freqrange[2], frequencies)
    local_frequencies = frequencies[allowed_rows]
    local_eligible = @view eligible[allowed_rows, :]
    local_costs = @view costs[allowed_rows, :]
    nstates = length(allowed_rows)

    previous = fill(T(Inf), nstates, nstates)
    current = similar(previous)
    predecessors = zeros(Int32, nstates, nstates, ntime)
    first_step = times[2] - times[1]
    for first_state in 1:nstates
        local_eligible[first_state, 1] || continue
        for second_state in 1:nstates
            local_eligible[second_state, 2] || continue
            first_frequency = local_frequencies[first_state]
            second_frequency = local_frequencies[second_state]
            _predictive_pair_allowed(
                first_frequency,
                second_frequency,
                first_step,
                config,
            ) || continue
            initial_slope = _physical_slope(first_frequency, second_frequency, first_step)
            previous[first_state, second_state] =
                local_costs[first_state, 1] + local_costs[second_state, 2] +
                _predictive_frequency_cost(first_frequency, second_frequency, config) +
                _predictive_initial_cost(first_frequency, initial_slope, config)
        end
    end
    any(isfinite, previous) || throw(ArgumentError(
        "no feasible predictive ridge reaches the second time sample",
    ))

    for column in 3:ntime
        fill!(current, T(Inf))
        previous_step = times[column - 1] - times[column - 2]
        next_step = times[column] - times[column - 1]
        for middle_state in 1:nstates
            local_eligible[middle_state, column - 1] || continue
            lower = 1
            upper = nstates
            if config.max_rate !== nothing
                maximum_jump = config.max_rate * next_step
                lower = searchsortedfirst(
                    local_frequencies,
                    local_frequencies[middle_state] - maximum_jump,
                )
                upper = searchsortedlast(
                    local_frequencies,
                    local_frequencies[middle_state] + maximum_jump,
                )
            end
            for next_state in lower:upper
                local_eligible[next_state, column] || continue
                best_cost = T(Inf)
                best_predecessor = 0
                predecessor_lower = 1
                predecessor_upper = nstates
                if config.method === :pruned
                    ratio = next_step / previous_step
                    predicted_lower = local_frequencies[next_state] -
                                      config.prediction_bandwidth
                    predicted_upper = local_frequencies[next_state] +
                                      config.prediction_bandwidth
                    first_lower_frequency = local_frequencies[middle_state] -
                        (predicted_upper - local_frequencies[middle_state]) / ratio
                    first_upper_frequency = local_frequencies[middle_state] -
                        (predicted_lower - local_frequencies[middle_state]) / ratio
                    predecessor_lower = searchsortedfirst(
                        local_frequencies,
                        first_lower_frequency,
                    )
                    predecessor_upper = searchsortedlast(
                        local_frequencies,
                        first_upper_frequency,
                    )
                    predecessor_lower > predecessor_upper && continue
                end
                for first_state in predecessor_lower:predecessor_upper
                    isfinite(previous[first_state, middle_state]) || continue
                    _predictive_transition_allowed(
                        local_frequencies[first_state],
                        local_frequencies[middle_state],
                        local_frequencies[next_state],
                        previous_step,
                        next_step,
                        config,
                    ) || continue
                    curvature_cost, _, _ = _predictive_curvature(
                        local_frequencies[first_state],
                        local_frequencies[middle_state],
                        local_frequencies[next_state],
                        previous_step,
                        next_step,
                        config,
                    )
                    candidate = previous[first_state, middle_state] +
                                _predictive_frequency_cost(
                                    local_frequencies[middle_state],
                                    local_frequencies[next_state],
                                    config,
                                ) + curvature_cost
                    if candidate < best_cost
                        best_cost = candidate
                        best_predecessor = first_state
                    end
                end
                if best_predecessor != 0
                    current[middle_state, next_state] =
                        best_cost + local_costs[next_state, column]
                    predecessors[middle_state, next_state, column] = best_predecessor
                end
            end
        end
        any(isfinite, current) || throw(ArgumentError(
            "no feasible predictive ridge reaches time index $column",
        ))
        previous, current = current, previous
    end

    best_previous = 0
    best_final = 0
    best_cost = T(Inf)
    for previous_state in 1:nstates, final_state in 1:nstates
        candidate = previous[previous_state, final_state]
        if candidate < best_cost
            best_cost = candidate
            best_previous = previous_state
            best_final = final_state
        end
    end
    best_final != 0 || throw(ArgumentError("no reachable predictive terminal state"))
    local_indices = zeros(Int, ntime)
    local_indices[end - 1] = best_previous
    local_indices[end] = best_final
    for column in ntime:-1:3
        local_indices[column - 2] = predecessors[
            local_indices[column - 1],
            local_indices[column],
            column,
        ]
    end
    indices = allowed_rows[local_indices]
    return _path_diagnostics(
        indices,
        frequencies,
        times,
        eligible,
        energies,
        global_reference,
        best_cost,
        freqrange,
        config,
    )
end

"""
    extract_predictive_ridge(coefficients, frequencies; times, freqrange, ...)

Track a frequency-constrained ridge with a constant-velocity second-order
model. `method=:exact` globally minimizes the full pair-state objective.
`method=:pruned` adds a fixed-Hz prediction window and is heuristic relative
to the unpruned objective. One- and two-sample inputs fall back to first-order
dense DP because curvature is not identifiable.
"""
function extract_predictive_ridge(
    coefficients::AbstractMatrix,
    frequencies::AbstractVector;
    times::AbstractVector,
    freqrange,
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
    epsilon::Real = eps(float(eltype(frequencies))),
    activity_threshold::Real = 0,
    valid_mask = nothing,
    frequency_envelope = nothing,
)
    T = promote_type(float(eltype(frequencies)), float(eltype(times)))
    converted_frequencies = T.(frequencies)
    converted_times = T.(times)
    _validate_ridge_inputs(
        coefficients,
        converted_frequencies,
        converted_times,
        freqrange,
        valid_mask,
    )
    config = PredictiveRidgeConfig(
        T;
        method,
        frequency_penalty,
        frequency_scale,
        curvature_penalty,
        curvature_scale,
        max_rate,
        max_chirp_rate_change,
        prediction_bandwidth,
        initial_frequency,
        initial_frequency_penalty,
        initial_frequency_scale,
        initial_slope,
        initial_slope_penalty,
        initial_slope_scale,
        epsilon,
        activity_threshold,
    )
    return _second_order_dynamic_ridge(
        coefficients,
        converted_frequencies,
        converted_times,
        (T(freqrange[1]), T(freqrange[2])),
        config;
        valid_mask,
        frequency_envelope,
    )
end

function extract_predictive_ridge(
    result::WSSTResult;
    freqrange,
    kwargs...,
)
    ntime = size(result.coefficients, 2)
    times = (0:ntime-1) .* result.metadata.sample_interval
    return extract_predictive_ridge(
        result.coefficients,
        result.frequencies;
        times,
        freqrange,
        kwargs...,
    )
end

function _predictive_config_keywords(config::PredictiveRidgeConfig)
    return (
        method = config.method,
        frequency_penalty = config.frequency_penalty,
        frequency_scale = config.frequency_scale,
        curvature_penalty = config.curvature_penalty,
        curvature_scale = config.curvature_scale,
        max_rate = config.max_rate,
        max_chirp_rate_change = config.max_chirp_rate_change,
        prediction_bandwidth = config.prediction_bandwidth,
        initial_frequency = config.initial_frequency,
        initial_frequency_penalty = config.initial_frequency_penalty,
        initial_frequency_scale = config.initial_frequency_scale,
        initial_slope = config.initial_slope,
        initial_slope_penalty = config.initial_slope_penalty,
        initial_slope_scale = config.initial_slope_scale,
        epsilon = config.epsilon,
        activity_threshold = config.activity_threshold,
    )
end

function extract_predictive_ridge(
    coefficients::AbstractMatrix,
    frequencies::AbstractVector,
    config::PredictiveRidgeConfig;
    times,
    freqrange,
    kwargs...,
)
    isempty(kwargs) || throw(ArgumentError(
        "a PredictiveRidgeConfig cannot be combined with configuration keywords",
    ))
    return extract_predictive_ridge(
        coefficients,
        frequencies;
        times,
        freqrange,
        _predictive_config_keywords(config)...,
    )
end

function extract_predictive_ridge(
    result::WSSTResult,
    config::PredictiveRidgeConfig;
    freqrange,
    kwargs...,
)
    isempty(kwargs) || throw(ArgumentError(
        "a PredictiveRidgeConfig cannot be combined with configuration keywords",
    ))
    return extract_predictive_ridge(
        result;
        freqrange,
        _predictive_config_keywords(config)...,
    )
end
