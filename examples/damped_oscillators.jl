#!/usr/bin/env julia

using CairoMakie
using SynchroSqueezing

normalized_similarity(reference, estimate) =
    sum(reference .* estimate) / sqrt(sum(abs2, reference) * sum(abs2, estimate))

relative_decibels(coefficients) =
    20 .* log10.(max.(abs.(coefficients) ./ maximum(abs, coefficients), 1e-6))

function main(output_path)
    fs = 512.0
    samples = 1536
    time = (0:samples-1) ./ fs
    damping_times = (1.15, 2.20)
    frequencies = (
        34.0 .+ 8.0 .* time,
        101.0 .- 6.0 .* time .+ 0.8 .* time .^ 2,
    )
    phases = (
        2π .* (34.0 .* time .+ 4.0 .* time .^ 2),
        2π .* (101.0 .* time .- 3.0 .* time .^ 2 .+ (0.8 / 3) .* time .^ 3),
    )
    analytic_truths = (
        1.0 .* exp.(-time ./ damping_times[1]) .* exp.(im .* phases[1]),
        0.72 .* exp.(-time ./ damping_times[2]) .* exp.(im .* phases[2]),
    )
    real_truths = (real.(analytic_truths[1]), real.(analytic_truths[2]))
    signal = real_truths[1] + real_truths[2]

    result = wsst(
        signal;
        fs,
        freqrange = (20.0, 115.0),
        nfrequencies = 571,
        threshold = RelativeThreshold(1e-9),
        preserve_residual = false,
    )
    ridges = extract_ridges(
        result;
        freqranges = [(28.0, 64.0), (80.0, 108.0)],
        penalty = 1.0,
        continuity_scale = 1.0,
        activity_threshold = 1e-8,
    )
    recovered = reconstruct_components(
        result,
        ridges;
        bandwidth_hz = [6.0, 6.0],
        ownership = :exclusive,
        analytic = true,
    )
    fit_window = (0.30, 2.55)
    features = [component_features(
        recovered.components[index], time;
        valid_mask = ridges[index].valid .& ridges[index].active,
    ) for index in eachindex(ridges)]
    fits = [fit_exponential_damping(feature; window = fit_window)
            for feature in features]

    interior = findall(t -> fit_window[1] <= t <= fit_window[2], time)
    reconstructed = real.(sum(recovered.components))
    total_metrics = component_metrics(reconstructed[interior], signal[interior])
    total_similarity = normalized_similarity(signal[interior], reconstructed[interior])
    ridge_statistics = [ridge_metrics(ridges[index], frequencies[index])
                        for index in eachindex(ridges)]
    component_statistics = [component_metrics(
        real.(recovered.components[index][interior]), real_truths[index][interior],
    ) for index in eachindex(ridges)]
    component_similarities = [normalized_similarity(
        real_truths[index][interior], real.(recovered.components[index][interior]),
    ) for index in eachindex(ridges)]

    sensitivity = [damping_sensitivity(
        result,
        ridges[index];
        bandwidths_hz = [5.0, 6.0, 7.0],
        windows = [(0.25, 2.45), fit_window, (0.35, 2.65)],
    ) for index in eachindex(ridges)]

    for index in eachindex(ridges)
        damping_error = abs(fits[index].damping_time / damping_times[index] - 1)
        println("mode $index")
        println("  ridge metrics: ", ridge_statistics[index])
        println("  damping time: truth=$(damping_times[index]) s, " *
                "estimate=$(fits[index].damping_time) s, " *
                "relative error=$(damping_error), " *
                "fit SE=$(fits[index].damping_time_standard_error) s")
        println("  bandwidth/window sensitivity range: ",
                sensitivity[index].minimum_damping_time, "–",
                sensitivity[index].maximum_damping_time, " s")
        println("  reconstruction metrics: ", component_statistics[index])
        println("  reconstruction similarity: ", component_similarities[index])
    end
    println("summed reconstruction metrics: ", total_metrics)
    println("summed reconstruction similarity: ", total_similarity)

    colors = Makie.wong_colors()[1:2]
    figure = Figure(size = (1800, 1500), fontsize = 18)

    signal_axis = Axis(
        figure[1, 1],
        ylabel = "amplitude",
        title = "Synthetic signal: two exponentially damped chirp oscillators",
    )
    lines!(signal_axis, time, signal; color = :black, linewidth = 1.2)
    hidexdecorations!(signal_axis; grid = false)

    ridge_axis = Axis(
        figure[2, 1],
        ylabel = "frequency (Hz)",
        title = "WSST and extracted ridges",
    )
    heatmap_plot = heatmap!(
        ridge_axis,
        time,
        result.frequencies,
        permutedims(relative_decibels(result.coefficients));
        colormap = :inferno,
        colorrange = (-60, 0),
    )
    for index in eachindex(ridges)
        lines!(
            ridge_axis, time, frequencies[index];
            color = colors[index], linewidth = 2, linestyle = :dash,
            label = "mode $index truth",
        )
        lines!(
            ridge_axis, time,
            ifelse.(ridges[index].active, ridges[index].frequencies, NaN);
            color = colors[index], linewidth = 3,
            label = "mode $index ridge (RMSE=$(round(ridge_statistics[index].rmse; digits=3)) Hz)",
        )
    end
    axislegend(ridge_axis; position = :rb, nbanks = 2)
    hidexdecorations!(ridge_axis; grid = false)
    Colorbar(figure[2, 2], heatmap_plot; label = "relative magnitude (dB)")

    reconstruction_axis = Axis(
        figure[3, 1],
        ylabel = "amplitude",
        title = "Separated reconstructions inside the damping-fit interval",
    )
    for index in eachindex(ridges)
        lines!(
            reconstruction_axis, time, real_truths[index];
            color = (colors[index], 0.55), linewidth = 2, linestyle = :dash,
            label = "mode $index truth",
        )
        lines!(
            reconstruction_axis, time, real.(recovered.components[index]);
            color = colors[index], linewidth = 1.2,
            label = "mode $index reconstruction: similarity " *
                    "$(round(component_similarities[index]; digits=4)), " *
                    "rel. L² $(round(component_statistics[index].relative_error; digits=4))",
        )
    end
    vspan!(reconstruction_axis, fit_window...; color = (:gray70, 0.12))
    axislegend(reconstruction_axis; position = :rt, nbanks = 2)
    hidexdecorations!(reconstruction_axis; grid = false)

    envelope_axis = Axis(
        figure[4, 1],
        xlabel = "time (s)",
        ylabel = "amplitude envelope",
        title = "Damping fits; summed reconstruction similarity = " *
                "$(round(total_similarity; digits=5)), relative L² error = " *
                "$(round(total_metrics.relative_error; digits=4))",
    )
    for index in eachindex(fits)
        relative_tau_error = 100abs(fits[index].damping_time / damping_times[index] - 1)
        lines!(
            envelope_axis, time, abs.(analytic_truths[index]);
            color = (colors[index], 0.55), linewidth = 2, linestyle = :dash,
            label = "mode $index true envelope, τ=$(damping_times[index]) s",
        )
        lines!(
            envelope_axis, time, features[index].amplitude;
            color = colors[index], linewidth = 1,
            label = "mode $index reconstructed envelope",
        )
        lines!(
            envelope_axis, time[fits[index].selected],
            fits[index].fitted_amplitude[fits[index].selected];
            color = colors[index], linewidth = 3, linestyle = :dot,
            label = "mode $index fit τ=$(round(fits[index].damping_time; digits=3)) s " *
                    "($(round(relative_tau_error; digits=2))% error)",
        )
    end
    vlines!(envelope_axis, collect(fit_window); color = (:gray30, 0.5), linestyle = :dash)
    axislegend(envelope_axis; position = :rt, nbanks = 2)

    linkxaxes!(signal_axis, ridge_axis, reconstruction_axis, envelope_axis)
    rowgap!(figure.layout, 14)
    save(output_path, figure; px_per_unit = 2)
    println("saved ", output_path)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_path = isempty(ARGS) ?
        joinpath(@__DIR__, "damped_oscillators.png") : abspath(first(ARGS))
    main(output_path)
end
