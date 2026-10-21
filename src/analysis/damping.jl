function _analysis_float_type(values...)
    types = map(value -> float(real(eltype(value))), values)
    return promote_type(types...)
end

function _validate_times(times, n)
    length(times) == n || throw(DimensionMismatch("times and samples must have equal length"))
    n >= 3 || throw(ArgumentError("at least three samples are required"))
    all(isfinite, times) || throw(ArgumentError("times must be finite"))
    all(diff(times) .> 0) || throw(ArgumentError("times must be strictly increasing"))
    return nothing
end

"""
    analytic_signal(signal)

Return the FFT-based analytic representation of a real signal. Complex input
is already analytic by contract and is copied unchanged.
"""
analytic_signal(signal::AbstractVector{<:Complex}) = collect(signal)

function analytic_signal(signal::AbstractVector{<:Real})
    n = length(signal)
    n > 0 || throw(ArgumentError("signal must not be empty"))
    spectrum = fft(signal)
    multiplier = zeros(real(eltype(spectrum)), n)
    multiplier[1] = one(eltype(multiplier))
    if iseven(n)
        n > 1 && (multiplier[2:n÷2] .= 2)
        n > 1 && (multiplier[n÷2 + 1] = 1)
    else
        multiplier[2:(n+1)÷2] .= 2
    end
    return ifft(spectrum .* multiplier)
end

function _unwrap_valid_phase!(phase, raw, valid)
    previous = zero(eltype(phase))
    offset = zero(eltype(phase))
    continuing = false
    for index in eachindex(raw)
        if !valid[index]
            phase[index] = eltype(phase)(NaN)
            continuing = false
            continue
        end
        value = raw[index]
        if continuing
            jump = value - previous
            offset -= 2 * eltype(phase)(π) * round(jump / (2 * eltype(phase)(π)))
        else
            offset = zero(eltype(phase))
        end
        phase[index] = value + offset
        previous = value
        continuing = true
    end
    return phase
end

function _differentiate_valid(values, times, valid)
    T = eltype(times)
    derivative = fill(T(NaN), length(values))
    for index in eachindex(values)
        valid[index] || continue
        left = index > firstindex(values) && valid[index - 1]
        right = index < lastindex(values) && valid[index + 1]
        if left && right
            derivative[index] = (values[index + 1] - values[index - 1]) /
                                (times[index + 1] - times[index - 1])
        elseif right
            derivative[index] = (values[index + 1] - values[index]) /
                                (times[index + 1] - times[index])
        elseif left
            derivative[index] = (values[index] - values[index - 1]) /
                                (times[index] - times[index - 1])
        end
    end
    return derivative
end

"""
    component_features(component, times; valid_mask, amplitude_floor,
                       estimate_frequency=true)

Extract an envelope and unwrapped phase from a complex component. Real input
is converted with [`analytic_signal`](@ref). Samples at or below
`amplitude_floor`, nonfinite samples, and samples rejected by `valid_mask` are
explicitly marked invalid.
"""
function component_features(
    component::AbstractVector,
    times::AbstractVector;
    valid_mask = trues(length(component)),
    amplitude_floor::Real = 0,
    estimate_frequency::Bool = true,
)
    _validate_times(times, length(component))
    length(valid_mask) == length(component) || throw(DimensionMismatch(
        "valid_mask must match the component length",
    ))
    amplitude_floor >= 0 && isfinite(amplitude_floor) || throw(ArgumentError(
        "amplitude_floor must be finite and nonnegative",
    ))
    T = _analysis_float_type(component, times)
    t = T.(times)
    analytic = Complex{T}.(analytic_signal(component))
    amplitude = T.(abs.(analytic))
    valid = BitVector(valid_mask .& isfinite.(real.(analytic)) .&
                      isfinite.(imag.(analytic)) .& (amplitude .> T(amplitude_floor)))
    phase = Vector{T}(undef, length(component))
    _unwrap_valid_phase!(phase, T.(angle.(analytic)), valid)
    frequency = estimate_frequency ? _differentiate_valid(phase, t, valid) ./ (2T(π)) : nothing
    return ComponentFeatures(t, amplitude, phase, frequency, valid)
end

component_features(
    component::ComponentReconstruction,
    times::AbstractVector;
    kwargs...,
) = component_features(component.signal, times; kwargs...)

function _selection_mask(amplitude, times, valid_mask, window, weights, amplitude_floor)
    n = length(amplitude)
    length(valid_mask) == n || throw(DimensionMismatch("valid_mask must match amplitude"))
    length(weights) == n || throw(DimensionMismatch("weights must match amplitude"))
    selected = BitVector(valid_mask .& isfinite.(amplitude) .& (amplitude .> amplitude_floor) .&
                         isfinite.(weights) .& (weights .> 0))
    if window !== nothing
        if window isa AbstractVector{Bool}
            length(window) == n || throw(DimensionMismatch("window mask must match amplitude"))
            selected .&= window
        elseif window isa AbstractUnitRange{<:Integer}
            indices = falses(n)
            checkbounds(indices, window)
            indices[window] .= true
            selected .&= indices
        elseif window isa Tuple && length(window) == 2
            first(window) <= last(window) || throw(ArgumentError("fit window must be ordered"))
            selected .&= (times .>= first(window)) .& (times .<= last(window))
        else
            throw(ArgumentError("window must be nothing, a Boolean mask, an index range, or (tmin, tmax)"))
        end
    end
    count(selected) >= 3 || throw(ArgumentError("the fit requires at least three selected samples"))
    return selected
end

function _fit_arrays(amplitude, times, selected, weights)
    indices = findall(selected)
    T = promote_type(eltype(amplitude), eltype(times), float(eltype(weights)))
    t = T.(times[indices])
    a = T.(amplitude[indices])
    w = T.(weights[indices])
    w ./= sum(w) / length(w)
    reference_time = sum(w .* t) / sum(w)
    return indices, t, a, w, reference_time, t .- reference_time
end

function _diagnostics(amplitude, times, selected, weights, reference_time, reference_amplitude,
                      decay_rate, method, status, rate_se, iterations, converged)
    T = promote_type(eltype(amplitude), eltype(times), typeof(reference_amplitude))
    fitted = T.(reference_amplitude .* exp.(-decay_rate .* (times .- reference_time)))
    residuals = T.(amplitude .- fitted)
    log_residuals = fill(T(NaN), length(amplitude))
    log_residuals[selected] .= log.(amplitude[selected]) .- log.(fitted[selected])
    w = T.(weights[selected])
    w ./= sum(w) / length(w)
    amp_residual = residuals[selected]
    log_residual = log_residuals[selected]
    rmse = sqrt(sum(w .* abs2.(amp_residual)) / sum(w))
    log_rmse = sqrt(sum(w .* abs2.(log_residual)) / sum(w))
    observed = T.(amplitude[selected])
    mean_observed = sum(w .* observed) / sum(w)
    total = sum(w .* (observed .- mean_observed) .^ 2)
    r_squared = iszero(total) ? T(NaN) : one(T) - sum(w .* amp_residual .^ 2) / total
    damping_time = decay_rate > zero(T) ? inv(decay_rate) : T(NaN)
    tau_se = decay_rate > zero(T) ? rate_se / decay_rate^2 : T(NaN)
    selected_times = times[selected]
    fit_window = (T(first(selected_times)), T(last(selected_times)))
    return DampingFit(
        damping_time, T(decay_rate), T(reference_amplitude), T(reference_time), fitted,
        residuals, log_residuals, selected, fit_window, method, status, rmse, log_rmse,
        r_squared, T(rate_se), T(tau_se), iterations, converged,
    )
end

function _loglinear_fit(amplitude, times, selected, weights, require_decay)
    _, _, a, w, reference_time, u = _fit_arrays(amplitude, times, selected, weights)
    y = log.(a)
    ybar = sum(w .* y) / sum(w)
    denominator = sum(w .* u .^ 2)
    denominator > 0 || throw(ArgumentError("selected fitting times have zero spread"))
    slope = sum(w .* u .* (y .- ybar)) / denominator
    decay_rate = -slope
    reference_amplitude = exp(ybar)
    log_residual = y .- (ybar .+ slope .* u)
    dof = length(y) - 2
    variance = sum(w .* log_residual .^ 2) / dof
    rate_se = sqrt(variance / denominator)
    status = require_decay && !(decay_rate > 0) ? :nondecaying : :ok
    return _diagnostics(amplitude, times, selected, weights, reference_time,
                        reference_amplitude, decay_rate, :loglinear, status,
                        rate_se, 1, true)
end

function _nonlinear_fit(amplitude, times, selected, weights, require_decay, max_iterations,
                        tolerance)
    _, _, observed, w, reference_time, u = _fit_arrays(amplitude, times, selected, weights)
    initial = _loglinear_fit(amplitude, times, selected, weights, false)
    log_amplitude = log(initial.reference_amplitude)
    decay_rate = initial.decay_rate
    require_decay && (decay_rate = max(decay_rate, zero(decay_rate)))
    converged = false
    iteration = 0
    for current_iteration in 1:max_iterations
        iteration = current_iteration
        model = exp.(log_amplitude .- decay_rate .* u)
        residual = observed .- model
        j1 = model
        j2 = .-u .* model
        h11 = sum(w .* j1 .* j1)
        h12 = sum(w .* j1 .* j2)
        h22 = sum(w .* j2 .* j2)
        determinant = h11 * h22 - h12^2
        determinant > eps(typeof(determinant)) * h11 * h22 || break
        g1 = sum(w .* j1 .* residual)
        g2 = sum(w .* j2 .* residual)
        delta1 = (h22 * g1 - h12 * g2) / determinant
        delta2 = (-h12 * g1 + h11 * g2) / determinant
        old_loss = sum(w .* residual .^ 2)
        step = one(decay_rate)
        accepted = false
        while step >= inv(typeof(step)(128))
            candidate_log_amplitude = log_amplitude + step * delta1
            candidate_rate = decay_rate + step * delta2
            require_decay && (candidate_rate = max(candidate_rate, zero(candidate_rate)))
            candidate = exp.(candidate_log_amplitude .- candidate_rate .* u)
            if sum(w .* (observed .- candidate) .^ 2) <= old_loss
                log_amplitude = candidate_log_amplitude
                decay_rate = candidate_rate
                accepted = true
                break
            end
            step /= 2
        end
        !accepted && break
        if max(abs(step * delta1), abs(step * delta2) * maximum(abs, u)) <= tolerance
            converged = true
            break
        end
    end
    model = exp.(log_amplitude .- decay_rate .* u)
    residual = observed .- model
    j1 = model
    j2 = .-u .* model
    h11 = sum(w .* j1 .* j1)
    h12 = sum(w .* j1 .* j2)
    h22 = sum(w .* j2 .* j2)
    determinant = h11 * h22 - h12^2
    variance = sum(w .* residual .^ 2) / (length(observed) - 2)
    rate_se = determinant > 0 ? sqrt(variance * h11 / determinant) : typeof(decay_rate)(NaN)
    status = require_decay && !(decay_rate > 0) ? :nondecaying :
             (converged ? :ok : :not_converged)
    return _diagnostics(amplitude, times, selected, weights, reference_time,
                        exp(log_amplitude), decay_rate, :nonlinear, status,
                        rate_se, iteration, converged)
end

"""
    fit_exponential_damping(amplitude, times; method=:loglinear, window=nothing,
                            weights, valid_mask, amplitude_floor=0,
                            require_decay=true)

Fit `A(t) = Aᵣ exp(-(t-tᵣ)/τ)`. `:loglinear` performs weighted regression in
log-amplitude; `:nonlinear` independently minimizes amplitude-space residuals.
The returned mask records every included sample. A constant or growing fit is
flagged with `status=:nondecaying` and a nonfinite damping time when
`require_decay=true`.
"""
function fit_exponential_damping(
    amplitude::AbstractVector{<:Real},
    times::AbstractVector;
    method::Symbol = :loglinear,
    window = nothing,
    weights = ones(_analysis_float_type(amplitude, times), length(amplitude)),
    valid_mask = trues(length(amplitude)),
    amplitude_floor::Real = 0,
    require_decay::Bool = true,
    max_iterations::Integer = 50,
    tolerance::Real = 1e-10,
)
    _validate_times(times, length(amplitude))
    amplitude_floor >= 0 && isfinite(amplitude_floor) || throw(ArgumentError(
        "amplitude_floor must be finite and nonnegative",
    ))
    method in (:loglinear, :nonlinear) || throw(ArgumentError(
        "method must be :loglinear or :nonlinear",
    ))
    max_iterations > 0 || throw(ArgumentError("max_iterations must be positive"))
    tolerance > 0 || throw(ArgumentError("tolerance must be positive"))
    T = _analysis_float_type(amplitude, times, weights)
    a = T.(amplitude)
    t = T.(times)
    w = T.(weights)
    selected = _selection_mask(a, t, valid_mask, window, w, T(amplitude_floor))
    return method === :loglinear ?
        _loglinear_fit(a, t, selected, w, require_decay) :
        _nonlinear_fit(a, t, selected, w, require_decay, max_iterations, T(tolerance))
end

fit_exponential_damping(features::ComponentFeatures; valid_mask = features.valid, kwargs...) =
    fit_exponential_damping(features.amplitude, features.times;
                            valid_mask = valid_mask, kwargs...)

"""Estimate the time-dependent rate `-d(log(A))/dt`; optionally smooth log-amplitude."""
function local_damping_rate(
    amplitude::AbstractVector{<:Real},
    times::AbstractVector;
    valid_mask = trues(length(amplitude)),
    amplitude_floor::Real = 0,
    smoothing_window::Integer = 1,
)
    _validate_times(times, length(amplitude))
    isodd(smoothing_window) && smoothing_window >= 1 || throw(ArgumentError(
        "smoothing_window must be a positive odd integer",
    ))
    T = _analysis_float_type(amplitude, times)
    a = T.(amplitude)
    t = T.(times)
    valid = BitVector(valid_mask .& isfinite.(a) .& (a .> T(amplitude_floor)))
    logarithm = fill(T(NaN), length(a))
    logarithm[valid] .= log.(a[valid])
    if smoothing_window > 1
        raw = copy(logarithm)
        radius = smoothing_window ÷ 2
        for index in eachindex(logarithm)
            valid[index] || continue
            lower = max(firstindex(logarithm), index - radius)
            upper = min(lastindex(logarithm), index + radius)
            neighbours = [j for j in lower:upper if valid[j]]
            logarithm[index] = sum(raw[neighbours]) / length(neighbours)
        end
    end
    return .-_differentiate_valid(logarithm, t, valid)
end

local_damping_rate(features::ComponentFeatures; valid_mask = features.valid, kwargs...) =
    local_damping_rate(features.amplitude, features.times;
                       valid_mask = valid_mask, kwargs...)

"""
    damping_sensitivity(result, ridges; bandwidths_hz, windows, method=:loglinear)

Repeat analytic ridge reconstruction and damping fitting over ridge candidates,
band widths, and fitting windows. The returned spread is reconstruction/window
sensitivity and is deliberately separate from each fit's standard error.
"""
function damping_sensitivity(
    result::WSSTResult,
    ridges::AbstractVector{<:RidgeResult};
    bandwidths_hz,
    windows,
    method::Symbol = :loglinear,
    amplitude_floor::Real = 0,
    fit_kwargs...,
)
    isempty(ridges) && throw(ArgumentError("at least one ridge is required"))
    isempty(bandwidths_hz) && throw(ArgumentError("at least one bandwidth is required"))
    isempty(windows) && throw(ArgumentError("at least one fit window is required"))
    T = eltype(result.frequencies)
    times = T.((0:size(result.coefficients, 2)-1) ./ result.metadata.sample_rate)
    cases = DampingSensitivityCase{T,DampingFit{T}}[]
    for (ridge_index, ridge) in pairs(ridges), bandwidth in bandwidths_hz
        reconstruction = reconstruct_component(
            result, ridge; bandwidth_hz = bandwidth, analytic = true,
        )
        valid = ridge.valid .& ridge.active
        features = component_features(
            reconstruction.signal, times;
            valid_mask = valid,
            amplitude_floor,
            estimate_frequency = false,
        )
        for window in windows
            fit = fit_exponential_damping(
                features; method, window, amplitude_floor, fit_kwargs...,
            )
            push!(cases, DampingSensitivityCase(
                ridge_index, T(bandwidth), fit.fit_window, fit,
            ))
        end
    end
    finite_times = [case.fit.damping_time for case in cases if isfinite(case.fit.damping_time)]
    isempty(finite_times) && return DampingSensitivity(cases, T(NaN), T(NaN), T(NaN))
    minimum_time, maximum_time = extrema(finite_times)
    return DampingSensitivity(cases, minimum_time, maximum_time, maximum_time - minimum_time)
end

damping_sensitivity(result::WSSTResult, ridge::RidgeResult; kwargs...) =
    damping_sensitivity(result, [ridge]; kwargs...)
