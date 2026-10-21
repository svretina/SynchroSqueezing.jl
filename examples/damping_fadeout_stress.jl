#!/usr/bin/env julia

using CairoMakie
using Random
using SynchroSqueezing

normalized_similarity(reference, estimate) =
    sum(reference .* estimate) / sqrt(sum(abs2, reference) * sum(abs2, estimate))

relative_decibels(coefficients) =
    20 .* log10.(max.(abs.(coefficients) ./ maximum(abs, coefficients), 1e-7))

function main(output_path)
    Random.seed!(404)
    fs = 512.0
    samples = 2048
    time = (0:samples-1) ./ fs
    noise_standard_deviation = 0.012
    damping_times = (0.32, 2.60)
    initial_amplitudes = (1.0, 0.58)
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
    clean_signal = real_truths[1] + real_truths[2]
    signal = clean_signal .+ noise_standard_deviation .* randn(samples)

    # Five input-noise standard deviations define the known synthetic
    # detectability boundary. It is an evaluation criterion, not fit data.
    amplitude_floor = 5noise_standard_deviation
    expected_masks = (
        BitVector(abs.(analytic_truths[1]) .>= amplitude_floor),
        BitVector(abs.(analytic_truths[2]) .>= amplitude_floor),
    )
    detectable_until = (
        time[findlast(expected_masks[1])],
        time[findlast(expected_masks[2])],
    )

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
    support_threshold = amplitude_floor^2
    support_ridges = extract_ridges(
        result;
        freqranges = [(27.0, 65.0), (78.0, 104.0)],
        penalty = 1.0,
        continuity_scale = 1.0,
        activity_threshold = support_threshold,
    )
    recovered = reconstruct_components(
        result,
        ridges;
        bandwidth_hz = [5.5, 5.5],
        ownership = :exclusive,
        analytic = true,
    )
    supported_recovered = reconstruct_components(
        result,
        support_ridges;
        bandwidth_hz = [5.5, 5.5],
        ownership = :exclusive,
        analytic = true,
    )
    features = [component_features(recovered.components[index], time)
                for index in eachindex(ridges)]
    supported_features = [component_features(
        supported_recovered.components[index], time;
        valid_mask = support_ridges[index].active,
    ) for index in eachindex(support_ridges)]

    safe_windows = ((0.12, detectable_until[1]), (0.20, 3.50))
    safe_fits = [fit_exponential_damping(
        supported_features[index];
        window = safe_windows[index],
        amplitude_floor,
    ) for index in eachindex(ridges)]
    naive_fits = [fit_exponential_damping(
        features[index];
        window = (safe_windows[index][1], 3.50),
    ) for index in eachindex(ridges)]

    fit_end_times = collect(range(0.55, 3.5; length = 30))
    fit_end_estimates = [fill(NaN, length(fit_end_times)) for _ in ridges]
    for index in eachindex(ridges), endpoint_index in eachindex(fit_end_times)
        endpoint = fit_end_times[endpoint_index]
        endpoint > safe_windows[index][1] || continue
        fit = fit_exponential_damping(
            features[index];
            window = (safe_windows[index][1], endpoint),
        )
        fit_end_estimates[index][endpoint_index] = fit.damping_time
    end

    ridge_statistics = [ridge_metrics(
        ridges[index], frequencies[index]; expected = expected_masks[index],
    ) for index in eachindex(ridges)]
    false_support = [
        count(support_ridges[index].active .& .!expected_masks[index]) /
        max(count(.!expected_masks[index]), 1) for index in eachindex(ridges)
    ]
    supported_coverage = [
        count(support_ridges[index].active .& expected_masks[index]) /
        count(expected_masks[index]) for index in eachindex(ridges)
    ]
    component_statistics = [component_metrics(
        real.(recovered.components[index][expected_masks[index]]),
        real_truths[index][expected_masks[index]],
    ) for index in eachindex(ridges)]
    component_similarities = [normalized_similarity(
        real_truths[index][expected_masks[index]],
        real.(recovered.components[index][expected_masks[index]]),
    ) for index in eachindex(ridges)]

    for index in eachindex(ridges)
        safe_error = abs(safe_fits[index].damping_time / damping_times[index] - 1)
        naive_error = abs(naive_fits[index].damping_time / damping_times[index] - 1)
        println("mode $index")
        println("  truth detectable through: $(detectable_until[index]) s")
        println("  ridge metrics while detectable: ", ridge_statistics[index])
        println("  thresholded ridge coverage while detectable: ", supported_coverage[index])
        println("  thresholded ridge activity fraction after truth falls below floor: ",
                false_support[index])
        println("  reliability-aware τ=$(safe_fits[index].damping_time) s, " *
                "relative error=$(safe_error)")
        println("  naive full-window τ=$(naive_fits[index].damping_time) s, " *
                "relative error=$(naive_error)")
        println("  reconstruction metrics while detectable: ", component_statistics[index])
        println("  reconstruction similarity while detectable: ",
                component_similarities[index])
    end

    colors = Makie.wong_colors()[1:2]
    figure = Figure(size = (1800, 1850), fontsize = 18)

    signal_axis = Axis(
        figure[1, 1], ylabel = "amplitude",
        title = "Fade-out stress test: τ₁=$(damping_times[1]) s, " *
                "τ₂=$(damping_times[2]) s, input noise σ=$noise_standard_deviation",
    )
    lines!(signal_axis, time, signal; color = (:gray25, 0.75), linewidth = 0.8,
           label = "noisy mixture")
    lines!(signal_axis, time, clean_signal; color = :black, linewidth = 1.1,
           label = "clean mixture")
    axislegend(signal_axis; position = :rt)
    hidexdecorations!(signal_axis; grid = false)

    ridge_axis = Axis(
        figure[2, 1], ylabel = "frequency (Hz)",
        title = "The short-lived ridge becomes noise-driven after its mode fades",
    )
    heatmap_plot = heatmap!(
        ridge_axis, time, result.frequencies,
        permutedims(relative_decibels(result.coefficients));
        colormap = :inferno, colorrange = (-70, 0),
    )
    for index in eachindex(ridges)
        lines!(ridge_axis, time, frequencies[index]; color = colors[index],
               linewidth = 2, linestyle = :dash, label = "mode $index truth")
        lines!(ridge_axis, time, ridges[index].frequencies; color = colors[index],
               linewidth = 1, label = "mode $index numerical path")
        lines!(ridge_axis, time,
               ifelse.(support_ridges[index].active,
                       support_ridges[index].frequencies, NaN);
               color = colors[index], linewidth = 3,
               label = "mode $index supported path")
    end
    vlines!(ridge_axis, [detectable_until[1]]; color = :white,
            linewidth = 2, linestyle = :dot, label = "mode 1 detectability limit")
    axislegend(ridge_axis; position = :rb, nbanks = 2)
    hidexdecorations!(ridge_axis; grid = false)
    Colorbar(figure[2, 2], heatmap_plot; label = "relative magnitude (dB)")

    reconstruction_axis = Axis(
        figure[3, 1], ylabel = "component amplitude",
        title = "Reconstructions and known synthetic detectability interval",
    )
    for index in eachindex(ridges)
        lines!(reconstruction_axis, time, real_truths[index]; color = (colors[index], 0.5),
               linewidth = 2, linestyle = :dash, label = "mode $index truth")
        lines!(reconstruction_axis, time, real.(recovered.components[index]);
               color = colors[index], linewidth = 1,
               label = "mode $index reconstructed; similarity=" *
                       "$(round(component_similarities[index]; digits=4))")
    end
    vlines!(reconstruction_axis, [detectable_until[1]]; color = :gray35,
            linewidth = 2, linestyle = :dot)
    axislegend(reconstruction_axis; position = :rt, nbanks = 2)
    hidexdecorations!(reconstruction_axis; grid = false)

    envelope_axis = Axis(
        figure[4, 1], ylabel = "envelope", yscale = log10,
        title = "Reliability-aware fits stop before the short-lived mode reaches the floor",
    )
    hlines!(envelope_axis, [amplitude_floor]; color = :gray25, linewidth = 2,
            linestyle = :dash, label = "known 5σ amplitude floor")
    for index in eachindex(ridges)
        lines!(envelope_axis, time, abs.(analytic_truths[index]);
               color = (colors[index], 0.5), linewidth = 2, linestyle = :dash,
               label = "mode $index true envelope")
        lines!(envelope_axis, time, features[index].amplitude;
               color = colors[index], linewidth = 1, label = "mode $index recovered")
        lines!(envelope_axis, time[safe_fits[index].selected],
               safe_fits[index].fitted_amplitude[safe_fits[index].selected];
               color = colors[index], linewidth = 3, linestyle = :dot,
               label = "mode $index safe fit τ=" *
                       "$(round(safe_fits[index].damping_time; digits=3)) s")
    end
    axislegend(envelope_axis; position = :rt, nbanks = 2)
    hidexdecorations!(envelope_axis; grid = false)

    sensitivity_axis = Axis(
        figure[5, 1], xlabel = "fit-window end time (s)",
        ylabel = "estimated damping time (s)",
        title = "Window-end stress test: fitting noise after fade-out biases τ₁",
    )
    for index in eachindex(ridges)
        hlines!(sensitivity_axis, [damping_times[index]]; color = (colors[index], 0.55),
                linewidth = 2, linestyle = :dash,
                label = "mode $index true τ")
        scatterlines!(sensitivity_axis, fit_end_times, fit_end_estimates[index];
                      color = colors[index], linewidth = 3, marker = :circle,
                      markersize = 7, label = "mode $index estimate")
    end
    vlines!(sensitivity_axis, [detectable_until[1]]; color = :gray25,
            linewidth = 2, linestyle = :dot, label = "mode 1 detectability limit")
    axislegend(sensitivity_axis; position = :rt, nbanks = 2)

    linkxaxes!(signal_axis, ridge_axis, reconstruction_axis, envelope_axis)
    rowgap!(figure.layout, 14)
    save(output_path, figure; px_per_unit = 2)
    println("saved ", output_path)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_path = isempty(ARGS) ?
        joinpath(@__DIR__, "damping_fadeout_stress.png") : abspath(first(ARGS))
    main(output_path)
end
