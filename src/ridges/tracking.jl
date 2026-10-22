function _validate_ridge_inputs(coefficients, frequencies, times, freqrange, valid_mask)
    nfrequency, ntime = size(coefficients)
    length(frequencies) == nfrequency || throw(DimensionMismatch(
        "frequency coordinates must match the first coefficient dimension",
    ))
    length(times) == ntime || throw(DimensionMismatch(
        "time coordinates must match the second coefficient dimension",
    ))
    nfrequency > 0 && ntime > 0 || throw(ArgumentError(
        "ridge extraction requires a nonempty frequency-by-time matrix",
    ))
    all(isfinite, frequencies) && all(>(0), frequencies) &&
        all(>(0), diff(frequencies)) || throw(ArgumentError(
        "ridge frequencies must be finite, positive, and strictly increasing",
    ))
    all(isfinite, times) && all(>(0), diff(times)) || throw(ArgumentError(
        "ridge times must be finite and strictly increasing",
    ))
    length(freqrange) == 2 || throw(ArgumentError(
        "freqrange must contain (fmin, fmax)",
    ))
    fmin, fmax = freqrange
    isfinite(fmin) && isfinite(fmax) && fmin <= fmax || throw(ArgumentError(
        "freqrange endpoints must be finite and ordered",
    ))
    any(frequency -> fmin <= frequency <= fmax, frequencies) || throw(ArgumentError(
        "freqrange contains no frequency-grid bins",
    ))
    if valid_mask !== nothing
        size(valid_mask) == size(coefficients) || throw(DimensionMismatch(
            "valid_mask must match the coefficient matrix",
        ))
    end
    return nothing
end

function _frequency_envelope(envelope, times, ::Type{T}) where {T<:AbstractFloat}
    envelope === nothing && return nothing
    values = envelope isa Function ? [envelope(time) for time in times] : collect(envelope)
    length(values) == length(times) || throw(DimensionMismatch(
        "frequency_envelope must provide one interval per time sample",
    ))
    return map(values) do bounds
        length(bounds) == 2 || throw(ArgumentError(
            "each frequency_envelope entry must contain (fmin, fmax)",
        ))
        lower, upper = T(bounds[1]), T(bounds[2])
        isfinite(lower) && isfinite(upper) && lower <= upper || throw(ArgumentError(
            "frequency_envelope intervals must be finite and ordered",
        ))
        (lower, upper)
    end
end

function _eligible_states(
    coefficients,
    frequencies,
    times,
    freqrange,
    valid_mask,
    frequency_envelope,
)
    eligible = falses(size(coefficients))
    envelope = _frequency_envelope(frequency_envelope, times, eltype(frequencies))
    lower, upper = freqrange
    for column in axes(coefficients, 2), row in axes(coefficients, 1)
        frequency = frequencies[row]
        inside = lower <= frequency <= upper
        if envelope !== nothing
            local_lower, local_upper = envelope[column]
            inside &= local_lower <= frequency <= local_upper
        end
        mask_valid = valid_mask === nothing || valid_mask[row, column]
        eligible[row, column] = inside && mask_valid && isfinite(coefficients[row, column])
    end
    return eligible
end

function _observation_costs(coefficients, eligible, config)
    T = typeof(config.epsilon)
    energies = Matrix{T}(undef, size(coefficients))
    costs = fill(T(Inf), size(coefficients))
    global_reference = zero(T)
    for index in eachindex(coefficients)
        energy = isfinite(coefficients[index]) ? T(abs2(coefficients[index])) : zero(T)
        energies[index] = energy
        global_reference = max(global_reference, energy)
    end
    for column in axes(coefficients, 2)
        reference = zero(T)
        for row in axes(coefficients, 1)
            eligible[row, column] && (reference = max(reference, energies[row, column]))
        end
        for row in axes(coefficients, 1)
            eligible[row, column] || continue
            costs[row, column] = iszero(reference) ? zero(T) :
                -log(energies[row, column] / reference + config.epsilon)
        end
    end
    return costs, energies, global_reference
end

_coordinate(frequency, ::Val{:physical}) = frequency
_coordinate(frequency, ::Val{:log}) = log(frequency)

function _transition_cost(previous, current, config::RidgeConfig)
    coordinate = Val(config.continuity)
    change = _coordinate(current, coordinate) - _coordinate(previous, coordinate)
    return config.penalty * abs2(change / config.continuity_scale)
end

function _transition_allowed(previous, current, delta_time, config::RidgeConfig)
    config.max_rate === nothing && return true
    return abs(current - previous) <= config.max_rate * delta_time
end

function _path_diagnostics(
    indices,
    frequencies,
    times,
    eligible,
    energies,
    global_reference,
    total_cost,
    freqrange,
    config,
)
    T = eltype(frequencies)
    ntime = length(times)
    path_frequencies = fill(T(NaN), ntime)
    valid = falses(ntime)
    active = falses(ntime)
    boundary = falses(ntime)
    local_energy = zeros(T, ntime)
    activity_floor = config.activity_threshold * global_reference
    for column in 1:ntime
        row = indices[column]
        row == 0 && continue
        valid[column] = true
        path_frequencies[column] = frequencies[row]
        local_energy[column] = energies[row, column]
        active[column] = local_energy[column] > zero(T) &&
                         local_energy[column] >= activity_floor
        first_admissible = findfirst(@view eligible[:, column])
        last_admissible = findlast(@view eligible[:, column])
        boundary[column] = row == first_admissible || row == last_admissible
    end
    return RidgeResult(
        indices,
        path_frequencies,
        collect(times),
        valid,
        active,
        boundary,
        local_energy,
        T(total_cost),
        (T(freqrange[1]), T(freqrange[2])),
        config,
    )
end

"""Return the first and last eligible row in each time column without allocating."""
function _eligible_bounds(eligible::AbstractMatrix)
    ntime = size(eligible, 2)
    lower = zeros(Int, ntime)
    upper = zeros(Int, ntime)
    for column in axes(eligible, 2)
        first_row = findfirst(@view eligible[:, column])
        first_row === nothing && continue
        lower[column] = first_row
        upper[column] = findlast(@view eligible[:, column])
    end
    return lower, upper
end

"""Precompute the time-independent continuity objective on the used row span."""
function _transition_costs(frequencies, first_row, last_row, config::RidgeConfig)
    T = eltype(frequencies)
    count = last_row - first_row + 1
    costs = Matrix{T}(undef, count, count)
    for current in first_row:last_row, previous in first_row:last_row
        costs[previous - first_row + 1, current - first_row + 1] =
            _transition_cost(frequencies[previous], frequencies[current], config)
    end
    return costs
end

const _MAX_TRANSITION_CACHE_BYTES = 64 * 1024^2

function _cache_transition_costs(frequencies, first_row, last_row, config)
    count = last_row - first_row + 1
    maximum_cached_rows = floor(Int, sqrt(_MAX_TRANSITION_CACHE_BYTES ÷ sizeof(eltype(frequencies))))
    return count <= maximum_cached_rows ?
        _transition_costs(frequencies, first_row, last_row, config) : nothing
end

function _maximum_ridge(costs, energies, frequencies, times, eligible, global_reference,
                        freqrange, config)
    indices = zeros(Int, length(times))
    total_cost = zero(eltype(frequencies))
    for column in axes(costs, 2)
        best_row = 0
        best_energy = -one(eltype(frequencies))
        for row in axes(costs, 1)
            eligible[row, column] || continue
            energy = energies[row, column]
            if energy > best_energy
                best_energy = energy
                best_row = row
            end
        end
        indices[column] = best_row
        best_row == 0 || (total_cost += costs[best_row, column])
    end
    return _path_diagnostics(indices, frequencies, times, eligible, energies,
                             global_reference, total_cost, freqrange, config)
end

function _prior_vector(prior, nfrequency, ::Type{T}, name) where {T}
    prior === nothing && return zeros(T, nfrequency)
    length(prior) == nfrequency || throw(DimensionMismatch(
        "$name must contain one value per frequency bin",
    ))
    converted = T.(prior)
    all(value -> isfinite(value) || value == T(Inf), converted) || throw(ArgumentError(
        "$name values must be finite or positive infinity",
    ))
    return converted
end

function _dynamic_ridge(
    costs,
    energies,
    frequencies,
    times,
    eligible,
    global_reference,
    freqrange,
    config;
    banded,
    initial_prior,
    end_prior,
)
    T = eltype(frequencies)
    nfrequency, ntime = size(costs)
    start_cost = _prior_vector(initial_prior, nfrequency, T, "initial_prior")
    terminal_cost = _prior_vector(end_prior, nfrequency, T, "end_prior")
    previous = fill(T(Inf), nfrequency)
    current = similar(previous)
    predecessors = zeros(Int, nfrequency, ntime)
    lower_eligible, upper_eligible = _eligible_bounds(eligible)
    any(iszero, lower_eligible) && throw(ArgumentError(
        "each time column must contain at least one eligible ridge state",
    ))
    first_used = minimum(lower_eligible)
    last_used = maximum(upper_eligible)
    transition_costs = _cache_transition_costs(
        frequencies, first_used, last_used, config,
    )
    for row in lower_eligible[1]:upper_eligible[1]
        eligible[row, 1] && (previous[row] = costs[row, 1] + start_cost[row])
    end
    any(isfinite, previous) || throw(ArgumentError("no valid initial ridge state"))

    for column in 2:ntime
        fill!(current, T(Inf))
        delta_time = times[column] - times[column - 1]
        for row in lower_eligible[column]:upper_eligible[column]
            eligible[row, column] || continue
            lower = lower_eligible[column - 1]
            upper = upper_eligible[column - 1]
            if banded && config.max_rate !== nothing
                maximum_jump = config.max_rate * delta_time
                lower = max(
                    lower,
                    searchsortedfirst(frequencies, frequencies[row] - maximum_jump),
                )
                upper = min(
                    upper,
                    searchsortedlast(frequencies, frequencies[row] + maximum_jump),
                )
            end
            lower <= upper || continue
            best_cost = T(Inf)
            best_predecessor = 0
            @inbounds for predecessor in lower:upper
                isfinite(previous[predecessor]) || continue
                _transition_allowed(
                    frequencies[predecessor], frequencies[row], delta_time, config,
                ) || continue
                transition = if transition_costs === nothing
                    _transition_cost(frequencies[predecessor], frequencies[row], config)
                else
                    transition_costs[
                        predecessor - first_used + 1,
                        row - first_used + 1,
                    ]
                end
                candidate = previous[predecessor] + transition
                if candidate < best_cost
                    best_cost = candidate
                    best_predecessor = predecessor
                end
            end
            if best_predecessor != 0
                current[row] = costs[row, column] + best_cost
                predecessors[row, column] = best_predecessor
            end
        end
        any(isfinite, current) || throw(ArgumentError(
            "no feasible ridge trajectory reaches time index $column",
        ))
        previous, current = current, previous
    end

    best_final = 0
    best_cost = T(Inf)
    for row in 1:nfrequency
        candidate = previous[row] + terminal_cost[row]
        if candidate < best_cost
            best_cost = candidate
            best_final = row
        end
    end
    best_final != 0 || throw(ArgumentError("no reachable terminal ridge state"))
    indices = zeros(Int, ntime)
    indices[end] = best_final
    for column in ntime:-1:2
        indices[column - 1] = predecessors[indices[column], column]
    end
    return _path_diagnostics(indices, frequencies, times, eligible, energies,
                             global_reference, best_cost, freqrange, config)
end

"""
    extract_ridge(coefficients, frequencies; times, freqrange, ...)

Track one ridge inside the mandatory inclusive `freqrange`. `method=:maximum`
is the column-wise baseline, `:dense` is the auditable global DP, and `:dp` is
the banded global DP. The objective uses per-column normalized negative
log-energy plus a quadratic continuity penalty in physical or log frequency.
"""
function extract_ridge(
    coefficients::AbstractMatrix,
    frequencies::AbstractVector;
    times::AbstractVector,
    freqrange,
    method::Symbol = :dp,
    continuity::Symbol = :physical,
    penalty::Real = 1,
    continuity_scale::Real = 1,
    max_rate::Union{Nothing,Real} = nothing,
    epsilon::Real = eps(float(eltype(frequencies))),
    activity_threshold::Real = 0,
    valid_mask = nothing,
    frequency_envelope = nothing,
    initial_prior = nothing,
    end_prior = nothing,
)
    T = promote_type(float(eltype(frequencies)), float(eltype(times)))
    converted_frequencies = T.(frequencies)
    converted_times = T.(times)
    _validate_ridge_inputs(
        coefficients, converted_frequencies, converted_times, freqrange, valid_mask,
    )
    config = RidgeConfig(
        T;
        method,
        continuity,
        penalty,
        continuity_scale,
        max_rate,
        epsilon,
        activity_threshold,
    )
    eligible = _eligible_states(
        coefficients,
        converted_frequencies,
        converted_times,
        (T(freqrange[1]), T(freqrange[2])),
        valid_mask,
        frequency_envelope,
    )
    costs, energies, global_reference = _observation_costs(coefficients, eligible, config)
    if method === :maximum
        return _maximum_ridge(
            costs, energies, converted_frequencies, converted_times, eligible,
            global_reference, freqrange, config,
        )
    end
    return _dynamic_ridge(
        costs,
        energies,
        converted_frequencies,
        converted_times,
        eligible,
        global_reference,
        freqrange,
        config;
        banded = method === :dp,
        initial_prior,
        end_prior,
    )
end

function extract_ridge(
    result::WSSTResult;
    freqrange,
    kwargs...,
)
    ntime = size(result.coefficients, 2)
    times = (0:ntime-1) .* result.metadata.sample_interval
    return extract_ridge(
        result.coefficients,
        result.frequencies;
        times,
        freqrange,
        kwargs...,
    )
end

function extract_ridge(
    coefficients::AbstractMatrix,
    frequencies::AbstractVector,
    config::RidgeConfig;
    times,
    freqrange,
    kwargs...,
)
    isempty(kwargs) || throw(ArgumentError(
        "a RidgeConfig cannot be combined with tracking-configuration keywords",
    ))
    return extract_ridge(
        coefficients,
        frequencies;
        times,
        freqrange,
        _config_keywords(config)...,
    )
end


function extract_ridge(
    result::WSSTResult,
    config::RidgeConfig;
    freqrange,
    kwargs...,
)
    isempty(kwargs) || throw(ArgumentError(
        "a RidgeConfig cannot be combined with tracking-configuration keywords",
    ))
    return extract_ridge(result; freqrange, _config_keywords(config)...)
end

function _config_keywords(config::RidgeConfig)
    return (
        method = config.method,
        continuity = config.continuity,
        penalty = config.penalty,
        continuity_scale = config.continuity_scale,
        max_rate = config.max_rate,
        epsilon = config.epsilon,
        activity_threshold = config.activity_threshold,
    )
end

"""
    extract_ridges(result; freqranges, configs=nothing, ...)

Track independent modes in the supplied frequency-range order. `configs` may
contain one `RidgeConfig` or named tuple per range; otherwise common keyword
settings are used.
"""
function extract_ridges(result::WSSTResult; freqranges, configs = nothing, kwargs...)
    ranges = collect(freqranges)
    isempty(ranges) && throw(ArgumentError("freqranges must not be empty"))
    if configs !== nothing
        length(configs) == length(ranges) || throw(DimensionMismatch(
            "configs must contain one configuration per frequency range",
        ))
        isempty(kwargs) || throw(ArgumentError(
            "common tracking keywords cannot be combined with per-range configs",
        ))
    end
    ridges = map(eachindex(ranges)) do index
        if configs === nothing
            extract_ridge(result; freqrange = ranges[index], kwargs...)
        else
            config = configs[index]
            settings = config isa RidgeConfig ? _config_keywords(config) : config
            settings isa NamedTuple || throw(ArgumentError(
                "each ridge config must be a RidgeConfig or named tuple",
            ))
            extract_ridge(result; freqrange = ranges[index], settings...)
        end
    end
    for first_index in eachindex(ridges), second_index in first_index+1:length(ridges)
        overlap = ridge_overlap(ridges[first_index], ridges[second_index])
        any(overlap) && @warn "independently tracked ridges overlap" first_index second_index
    end
    return ridges
end

"""Return samples where two active trajectories are within `tolerance_hz`."""
function ridge_overlap(first::RidgeResult, second::RidgeResult; tolerance_hz::Real = 0)
    length(first.frequencies) == length(second.frequencies) || throw(DimensionMismatch(
        "ridge trajectories must have the same number of samples",
    ))
    tolerance_hz >= 0 || throw(ArgumentError("tolerance_hz must be nonnegative"))
    return first.active .& second.active .&
           (abs.(first.frequencies .- second.frequencies) .<= tolerance_hz)
end
