#!/usr/bin/env julia

using CairoMakie
using Random
using SynchroSqueezing

relative_decibels(coefficients) =
    20 .* log10.(max.(abs.(coefficients) ./ maximum(abs, coefficients), 1e-7))

function finite_median(values)
    finite = sort(filter(isfinite, collect(values)))
    isempty(finite) && return NaN
    middle = length(finite) ÷ 2
    return isodd(length(finite)) ? finite[middle + 1] :
           (finite[middle] + finite[middle + 1]) / 2
end

function main(output_path)
    Random.seed!(404)
    fs = 512.0
    samples = 2048
    time = (0:samples-1) ./ fs
    damping_times = (0.32, 2.60)
    initial_amplitudes = (1.0, 0.58)
    injected_background = 0.012
    frequencies = (
        32.0 .+ 7.0 .* time,
        98.0 .- 4.0 .* time .+ 0.45 .* time .^ 2,
    )
    phases = (
        2π .* (32.0 .* time .+ 3.5 .* time .^ 2),
        2π .* (98.0 .* time .- 2.0 .* time .^ 2 .+ 0.15 .* time .^ 3),
    )
    analytic_truths = (
        initial_amplitudes[1] .* exp.(-time ./ damping_times[1]) .* exp.(im .* phases[1]),
        initial_amplitudes[2] .* exp.(-time ./ damping_times[2]) .* exp.(im .* phases[2]),
    )
    real_truths = (real.(analytic_truths[1]), real.(analytic_truths[2]))
    signal = real_truths[1] + real_truths[2] .+
             injected_background .* randn(samples)

    result = wsst(
        signal;
        fs,
        freqrange = (20.0, 110.0),
        nfrequencies = 541,
        threshold = RelativeThreshold(2e-4),
        preserve_residual = false,
    )
    ridges = extract_ridges(
        result;
        freqranges = [(27.0, 65.0), (78.0, 104.0)],
        penalty = 1.0,
        continuity_scale = 1.0,
        activity_threshold = 0,
    )
    recovered = reconstruct_components(
        result,
        ridges;
        bandwidth_hz = 5.5,
        ownership = :exclusive,
        analytic = true,
    )
    features = [component_features(recovered.components[index], time)
                for index in eachindex(ridges)]
    prominences = [ridge_prominence(
        result,
        ridges[index];
        ridge_bandwidth_hz = 1.0,
        background_band_hz = (2.5, 7.0),
    ) for index in eachindex(ridges)]

    boundary_mask = BitVector((time .>= 0.10) .& (time .<= 3.85))
    reference_masks = (
        BitVector((time .>= 0.10) .& (time .<= 0.30)),
        BitVector((time .>= 0.20) .& (time .<= 0.60)),
    )
    relative_thresholds = (0.10, 0.05, 0.02)
    support_sensitivity = [mode_support_sensitivity(
        features[index];
        relative_amplitude_thresholds = relative_thresholds,
        persistence_seconds = 0.08,
        recovery_factor = 1.25,
        reference_mask = reference_masks[index],
        boundary_mask,
    ) for index in eachindex(ridges)]
    amplitude_support = [
        [case.support for case in support_sensitivity[index].cases]
        for index in eachindex(ridges)
    ]

    prominence_thresholds = (10.0, 100.0, 1000.0)
    prominence_only_support = [
        [mode_support(
            features[index];
            prominence = prominences[index],
            relative_amplitude_threshold = nothing,
            prominence_threshold = threshold,
            persistence_seconds = 0.08,
            recovery_factor = 1.25,
            reference_mask = reference_masks[index],
            boundary_mask,
        ) for threshold in prominence_thresholds]
        for index in eachindex(ridges)
    ]

    # The prominence threshold is deliberately an exposed experiment setting.
    # It is local WSST contrast, not a detector-noise likelihood.
    prominence_threshold = 100.0
    operational_amplitude_threshold = 0.05
    combined_support = [mode_support(
        features[index];
        prominence = prominences[index],
        relative_amplitude_threshold = operational_amplitude_threshold,
        prominence_threshold,
        persistence_seconds = 0.08,
        recovery_factor = 1.25,
        reference_mask = reference_masks[index],
        boundary_mask,
    ) for index in eachindex(ridges)]
    supported_components = [supported_reconstruction(
        recovered.components[index], combined_support[index],
    ) for index in eachindex(ridges)]

    fits = [fit_exponential_damping(
        features[index]; valid_mask = combined_support[index].fit_mask,
    ) for index in eachindex(ridges)]
    naive_fits = [fit_exponential_damping(
        features[index]; window = (0.12, 3.50),
    ) for index in eachindex(ridges)]

    # The injected-background crossing is evaluation truth only. It is not an
    # input to mode_support or to either activity criterion.
    evaluation_floor = 5injected_background
    evaluation_stop = [damping_times[index] *
        log(initial_amplitudes[index] / evaluation_floor) for index in eachindex(ridges)]

    for index in eachindex(ridges)
        println("mode $index")
        println("  component reference amplitude: ",
                combined_support[index].config.reference_amplitude)
        for (threshold_index, threshold) in pairs(relative_thresholds)
            support = amplitude_support[index][threshold_index]
            threshold_fit = support_sensitivity[index].cases[threshold_index].damping_fit
            ideal_stop = damping_times[index] * log(
                initial_amplitudes[index] /
                (threshold * support.config.reference_amplitude),
            )
            println("  amplitude-only ε=$(threshold): stop=$(support.stop_time) s, " *
                    "ideal reconstructed-reference crossing=$(ideal_stop) s, " *
                    "fitted τ=$(threshold_fit.damping_time) s")
        end
        println("  damping spread over 2%, 5%, 10% thresholds: ",
                support_sensitivity[index].damping_time_spread, " s")
        for (threshold_index, threshold) in pairs(prominence_thresholds)
            prominence_support = prominence_only_support[index][threshold_index]
            println("  prominence-only threshold=$(threshold): stop=",
                    prominence_support.stop_time, " s (",
                    prominence_support.stop_reason, ")")
        end
        println("  median prominence before 5σ evaluation crossing: ", finite_median(
            prominences[index].ratio[time .<= min(evaluation_stop[index], time[end])],
        ))
        println("  median prominence after 5σ evaluation crossing: ", finite_median(
            prominences[index].ratio[time .> min(evaluation_stop[index], time[end])],
        ))
        println("  combined 5% + prominence stop: ",
                combined_support[index].stop_time, " s (",
                combined_support[index].stop_reason, ")")
        println("  independent 5σ evaluation crossing: ", evaluation_stop[index], " s")
        println("  supported damping τ=$(fits[index].damping_time) s, " *
                "relative error=$(abs(fits[index].damping_time / damping_times[index] - 1))")
        println("  naive damping τ=$(naive_fits[index].damping_time) s, " *
                "relative error=$(abs(naive_fits[index].damping_time / damping_times[index] - 1))")
    end

    colors = Makie.wong_colors()[1:2]
    threshold_colors = [:seagreen3, :darkorange2, :mediumpurple3]
    threshold_markers = [:circle, :rect, :diamond]
    figure = Figure(size = (1800, 1900), fontsize = 18)

    ridge_axis = Axis(
        figure[1:2, 1], ylabel = "frequency (Hz)",
        title = "Raw ridges and operationally supported prefixes",
    )
    heatmap_plot = heatmap!(
        ridge_axis, time, result.frequencies,
        permutedims(relative_decibels(result.coefficients));
        colormap = :inferno, colorrange = (-70, 0),
    )
    for index in eachindex(ridges)
        lines!(ridge_axis, time, frequencies[index]; color = colors[index],
               linewidth = 2, linestyle = :dash, label = "mode $index truth")
        lines!(ridge_axis, time, ridges[index].frequencies; color = (colors[index], 0.45),
               linewidth = 1.2, label = "mode $index raw DP path")
        lines!(ridge_axis, time,
               ifelse.(combined_support[index].supported,
                       ridges[index].frequencies, NaN);
               color = colors[index], linewidth = 3,
               label = "mode $index supported")
    end
    axislegend(ridge_axis; position = :rb, nbanks = 2)
    hidexdecorations!(ridge_axis; grid = false)
    Colorbar(figure[1:2, 2], heatmap_plot; label = "relative magnitude (dB)")

    amplitude_axis = Axis(
        figure[3, 1], ylabel = "relative envelope", yscale = log10,
        title = "Per-component relative amplitude and amplitude-only stopping times",
    )
    for (threshold_index, threshold) in pairs(relative_thresholds)
        hlines!(amplitude_axis, [threshold]; color = threshold_colors[threshold_index],
                linestyle = :dash, linewidth = 1.5, label = "ε=$(threshold)")
    end
    for index in eachindex(features)
        lines!(amplitude_axis, time, combined_support[index].relative_amplitude;
               color = colors[index], linewidth = 1.5,
               label = "mode $index recovered / own reference")
        for (threshold_index, _) in pairs(relative_thresholds)
            stop = amplitude_support[index][threshold_index].stop_time
            stop === nothing || scatter!(amplitude_axis, [stop],
                [relative_thresholds[threshold_index]];
                color = colors[index], marker = threshold_markers[threshold_index],
                markersize = 13)
        end
    end
    axislegend(amplitude_axis; position = :rt, nbanks = 2)
    hidexdecorations!(amplitude_axis; grid = false)

    prominence_axis = Axis(
        figure[4, 1], ylabel = "local prominence", yscale = log10,
        title = "Local ridge energy / sideband energy (diagnostic, not detector SNR)",
    )
    hlines!(prominence_axis, [prominence_threshold]; color = :gray25,
            linewidth = 2, linestyle = :dash, label = "trial threshold")
    for index in eachindex(prominences)
        lines!(prominence_axis, time, prominences[index].ratio;
               color = colors[index], linewidth = 1.2, label = "mode $index")
        stop = combined_support[index].stop_time
        stop === nothing || vlines!(prominence_axis, [stop]; color = colors[index],
                                    linewidth = 2, linestyle = :dot,
                                    label = "mode $index combined stop")
    end
    axislegend(prominence_axis; position = :rt, nbanks = 2)
    hidexdecorations!(prominence_axis; grid = false)

    reconstruction_axis = Axis(
        figure[5, 1], ylabel = "component amplitude",
        title = "Reconstruction is retained only over the supported prefix",
    )
    for index in eachindex(ridges)
        lines!(reconstruction_axis, time, real_truths[index];
               color = (colors[index], 0.45), linestyle = :dash, linewidth = 2,
               label = "mode $index truth")
        lines!(reconstruction_axis, time, real.(recovered.components[index]);
               color = (colors[index], 0.25), linewidth = 0.8,
               label = "mode $index raw reconstruction")
        lines!(reconstruction_axis, time,
               ifelse.(combined_support[index].supported,
                       real.(supported_components[index].signal), NaN);
               color = colors[index], linewidth = 1.4,
               label = "mode $index supported reconstruction")
    end
    axislegend(reconstruction_axis; position = :rt, nbanks = 2)
    hidexdecorations!(reconstruction_axis; grid = false)

    damping_axis = Axis(
        figure[6, 1], xlabel = "time (s)", ylabel = "envelope", yscale = log10,
        title = "Damping fits restricted by operational support",
    )
    for index in eachindex(fits)
        lines!(damping_axis, time, abs.(analytic_truths[index]);
               color = (colors[index], 0.45), linestyle = :dash, linewidth = 2,
               label = "mode $index true envelope, τ=$(damping_times[index]) s")
        lines!(damping_axis, time[fits[index].selected],
               fits[index].fitted_amplitude[fits[index].selected];
               color = colors[index], linewidth = 3,
               label = "mode $index supported fit, τ=" *
                       "$(round(fits[index].damping_time; digits=3)) s")
    end
    axislegend(damping_axis; position = :rt, nbanks = 2)

    linkxaxes!(ridge_axis, amplitude_axis, prominence_axis,
               reconstruction_axis, damping_axis)
    rowgap!(figure.layout, 13)
    save(output_path, figure; px_per_unit = 2)
    println("saved ", output_path)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_path = isempty(ARGS) ?
        joinpath(@__DIR__, "mode_termination_two_components.png") : abspath(first(ARGS))
    main(output_path)
end
