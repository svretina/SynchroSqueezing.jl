#!/usr/bin/env julia

using CairoMakie
using SynchroSqueezing

function nonexponential_case()
    fs = 256.0
    samples = 768
    time = (0:samples-1) ./ fs
    tau1 = 0.70
    tau2 = 0.90
    threshold = 0.05
    amplitude = exp.(-time ./ tau1 .- time .^ 2 ./ (2tau2^2))
    phase = 2π .* (34.0 .* time .+ 2.5 .* time .^ 2)
    signal = amplitude .* cos.(phase)
    result = wsst(
        signal;
        fs,
        freqrange = (25.0, 55.0),
        nfrequencies = 241,
        threshold = RelativeThreshold(1e-9),
        preserve_residual = false,
    )
    ridge = extract_ridge(
        result;
        freqrange = (30.0, 52.0),
        penalty = 1.0,
        continuity_scale = 1.0,
    )
    raw = reconstruct_component(result, ridge; bandwidth_hz = 5.0, analytic = true)
    features = component_features(raw, time)
    boundary = BitVector((time .>= 0.10) .& (time .<= 2.85))
    reference = BitVector((time .>= 0.10) .& (time .<= 0.30))
    support = mode_support(
        features;
        relative_amplitude_threshold = threshold,
        persistence_seconds = 0.06,
        recovery_factor = 1.2,
        reference_mask = reference,
        boundary_mask = boundary,
    )
    supported = supported_reconstruction(raw, support)
    expected_stop = tau2^2 * (
        -inv(tau1) + sqrt(inv(tau1)^2 + 2log(inv(threshold)) / tau2^2)
    )
    local_rate = local_damping_rate(features; smoothing_window =  nine_or_one(fs))
    return (; time, amplitude, features, support, supported, expected_stop,
            local_rate, tau1, tau2)
end

# Roughly 35 ms, always an odd sample count.
function nine_or_one(fs)
    samples = max(1, round(Int, 0.035fs))
    return isodd(samples) ? samples : samples + 1
end

function leakage_case(separation)
    fs = 256.0
    samples = 768
    time = (0:samples-1) ./ fs
    short_tau = 0.40
    short_amplitude = 0.35
    long_tau = 3.0
    long_amplitude = 1.20
    short_frequency = 35.0
    long_frequency = short_frequency + separation
    short_component = short_amplitude .* exp.(-time ./ short_tau) .*
                      cos.(2π * short_frequency .* time)
    long_component = long_amplitude .* exp.(-time ./ long_tau) .*
                     cos.(2π * long_frequency .* time)
    result = wsst(
        short_component + long_component;
        fs,
        freqrange = (25.0, 65.0),
        nfrequencies = 321,
        threshold = RelativeThreshold(1e-10),
        preserve_residual = false,
    )
    ridge = extract_ridge(
        result;
        freqrange = (32.0, 38.0),
        penalty = 1.0,
        continuity_scale = 1.0,
    )
    raw = reconstruct_component(result, ridge; bandwidth_hz = 2.0, analytic = true)
    features = component_features(raw, time)
    prominence = ridge_prominence(
        result, ridge; ridge_bandwidth_hz = 0.75,
        background_band_hz = (2.25, 5.5),
    )
    boundary = BitVector((time .>= 0.10) .& (time .<= 2.85))
    support = mode_support(
        features;
        relative_amplitude_threshold = 0.05,
        persistence_seconds = 0.06,
        recovery_factor = 1.2,
        reference_amplitude = short_amplitude,
        boundary_mask = boundary,
    )
    fit = fit_exponential_damping(features; valid_mask = support.fit_mask)
    expected_stop = short_tau * log(20)
    return (; separation, stop_time = support.stop_time,
            stop_reason = support.stop_reason,
            damping_time = fit.damping_time,
            damping_error = abs(fit.damping_time / short_tau - 1),
            expected_stop,
            stopping_time_error = support.stop_time === nothing ?
                time[end] - expected_stop : support.stop_time - expected_stop,
            false_continuation = support.stop_time === nothing ?
                time[end] - expected_stop : max(support.stop_time - expected_stop, 0),
            premature_termination = support.stop_time === nothing ? 0.0 :
                max(expected_stop - support.stop_time, 0),
            prominence_at_expected = prominence.ratio[argmin(abs.(time .- expected_stop))])
end

function main(output_path)
    nonexponential = nonexponential_case()
    separations = [4.0, 6.0, 8.0, 12.0, 18.0]
    leakage = [leakage_case(separation) for separation in separations]

    println("nonexponential decay")
    println("  expected 5% crossing: ", nonexponential.expected_stop, " s")
    println("  recovered stop: ", nonexponential.support.stop_time, " s (",
            nonexponential.support.stop_reason, ")")
    println("  stopping-time error: ",
            nonexponential.support.stop_time - nonexponential.expected_stop, " s")
    println("strong-mode leakage sweep")
    for case in leakage
        println("  Δf=$(case.separation) Hz: stop=$(case.stop_time) s (",
                case.stop_reason, "), false continuation=$(case.false_continuation) s, " *
                "premature termination=$(case.premature_termination) s, " *
                "τ error=$(case.damping_error), " *
                "prominence at expected stop=$(case.prominence_at_expected)")
    end

    figure = Figure(size = (1700, 1250), fontsize = 18)
    envelope_axis = Axis(
        figure[1, 1], ylabel = "relative envelope", yscale = log10,
        title = "Nonexponential decay: support follows the observed envelope",
    )
    lines!(envelope_axis, nonexponential.time, nonexponential.amplitude;
           color = :black, linestyle = :dash, linewidth = 2,
           label = "true nonexponential envelope")
    lines!(envelope_axis, nonexponential.time,
           nonexponential.support.relative_amplitude;
           color = :dodgerblue3, linewidth = 1.5, label = "reconstructed / reference")
    hlines!(envelope_axis, [0.05]; color = :darkorange2, linestyle = :dash,
            label = "5% threshold")
    vlines!(envelope_axis, [nonexponential.expected_stop]; color = :black,
            linestyle = :dot, linewidth = 2, label = "analytical crossing")
    vlines!(envelope_axis, [nonexponential.support.stop_time]; color = :dodgerblue3,
            linestyle = :dot, linewidth = 2, label = "operational stop")
    axislegend(envelope_axis; position = :rt, nbanks = 2)
    hidexdecorations!(envelope_axis; grid = false)

    rate_axis = Axis(
        figure[2, 1], xlabel = "time (s)", ylabel = "local decay rate (s⁻¹)",
        title = "Increasing local rate confirms that termination does not assume constant τ",
    )
    true_rate = inv(nonexponential.tau1) .+
                nonexponential.time ./ nonexponential.tau2^2
    lines!(rate_axis, nonexponential.time, true_rate; color = :black,
           linestyle = :dash, linewidth = 2, label = "true local rate")
    lines!(rate_axis, nonexponential.time, nonexponential.local_rate;
           color = :dodgerblue3, linewidth = 1.2, label = "reconstructed local rate")
    vlines!(rate_axis, [nonexponential.support.stop_time]; color = :dodgerblue3,
            linestyle = :dot, linewidth = 2)
    axislegend(rate_axis; position = :lt)
    linkxaxes!(envelope_axis, rate_axis)

    continuation_axis = Axis(
        figure[1, 2], xlabel = "frequency separation (Hz)",
        ylabel = "stopping-time error (s)",
        title = "Leakage from a 3.4× stronger long-lived mode\n" *
                "negative values are premature termination",
    )
    scatterlines!(continuation_axis, separations,
                  [case.stopping_time_error for case in leakage];
                  color = :firebrick3, linewidth = 2.5, markersize = 11)
    hlines!(continuation_axis, [0.0]; color = :gray35, linestyle = :dash)

    error_axis = Axis(
        figure[2, 2], xlabel = "frequency separation (Hz)",
        ylabel = "relative damping-time error", yscale = log10,
        title = "Downstream damping sensitivity to mode separation",
    )
    scatterlines!(error_axis, separations, [case.damping_error for case in leakage];
                  color = :mediumpurple3, linewidth = 2.5, markersize = 11)

    colgap!(figure.layout, 28)
    rowgap!(figure.layout, 20)
    save(output_path, figure; px_per_unit = 2)
    println("saved ", output_path)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_path = isempty(ARGS) ?
        joinpath(@__DIR__, "mode_termination_validation.png") : abspath(first(ARGS))
    main(output_path)
end
