#!/usr/bin/env julia

using CairoMakie
using LinearAlgebra: dot, norm
using SynchroSqueezing

function relative_decibels(coefficients)
    relative = abs.(coefficients) ./ maximum(abs, coefficients)
    return 20 .* log10.(max.(relative, 1e-6))
end

normalized_similarity(reference, estimate) =
    dot(reference, estimate) / (norm(reference) * norm(estimate))

function plot_results(
    output_path,
    time,
    signal,
    result,
    ridges,
    true_frequencies,
    reconstructed,
    interior_interval,
    full_metrics,
    interior_metrics,
    full_similarity,
    interior_similarity,
)
    colors = Makie.wong_colors()[1:2]
    decibels = relative_decibels(result.coefficients)
    error_signal = reconstructed .- signal
    figure = Figure(size = (1800, 1300), fontsize = 19)

    ridge_axis = Axis(
        figure[1, 1],
        ylabel = "frequency (Hz)",
        title = "Frequency-constrained DP ridges on the two-chirp WSST",
    )
    image = heatmap!(
        ridge_axis,
        time,
        result.frequencies,
        permutedims(decibels);
        colormap = :inferno,
        colorrange = (-60, 0),
    )
    search_ranges = ((30.0, 60.0), (90.0, 115.0))
    for index in eachindex(ridges)
        hlines!(
            ridge_axis,
            collect(search_ranges[index]);
            color = (colors[index], 0.55),
            linewidth = 1.2,
            linestyle = :dot,
        )
        lines!(
            ridge_axis,
            time,
            true_frequencies[index];
            color = colors[index],
            linewidth = 2.0,
            linestyle = :dash,
            label = "mode $index truth",
        )
        measured_frequency = ifelse.(
            ridges[index].active,
            ridges[index].frequencies,
            NaN,
        )
        lines!(
            ridge_axis,
            time,
            measured_frequency;
            color = colors[index],
            linewidth = 3.0,
            label = "mode $index extracted ridge",
        )
    end
    axislegend(ridge_axis; position = :rb, nbanks = 2)
    hidexdecorations!(ridge_axis; grid = false)
    Colorbar(figure[1, 2], image; label = "magnitude relative to maximum (dB)")

    comparison_axis = Axis(
        figure[2, 1],
        ylabel = "amplitude",
        title = "Original two-chirp signal versus sum of reconstructed components",
    )
    lines!(
        comparison_axis,
        time,
        signal;
        color = :black,
        linewidth = 1.4,
        label = "original signal",
    )
    lines!(
        comparison_axis,
        time,
        reconstructed;
        color = :dodgerblue3,
        linewidth = 1.2,
        linestyle = :dash,
        label = "reconstructed components",
    )
    axislegend(comparison_axis; position = :rt, orientation = :horizontal)
    hidexdecorations!(comparison_axis; grid = false)

    error_axis = Axis(
        figure[3, 1],
        xlabel = "time (s)",
        ylabel = "reconstruction error",
        title = "Full: similarity = $(round(full_similarity; digits = 6)), " *
                "relative L² error = $(round(full_metrics.relative_error; digits = 5)), " *
                "SRER = $(round(full_metrics.srer_db; digits = 2)) dB\n" *
                "Interior ($(round(interior_interval[1]; digits = 3))–" *
                "$(round(interior_interval[2]; digits = 3)) s): " *
                "similarity = $(round(interior_similarity; digits = 6)), " *
                "relative L² error = $(round(interior_metrics.relative_error; digits = 5)), " *
                "SRER = $(round(interior_metrics.srer_db; digits = 2)) dB",
    )
    lines!(error_axis, time, error_signal; color = :firebrick, linewidth = 1.1)

    linkxaxes!(ridge_axis, comparison_axis, error_axis)
    rowgap!(figure.layout, 16)
    save(output_path, figure; px_per_unit = 2.0)
    return nothing
end

function main(output_path)
    fs = 512.0
    samples = 1024
    time = (0:samples-1) ./ fs
    frequency1 = 38.0 .+ 7.0 .* time
    frequency2 = 105.0 .- 5.0 .* time
    component1 = cos.(2π .* (38.0 .* time .+ 3.5 .* time .^ 2))
    component2 = 0.55 .* cos.(2π .* (105.0 .* time .- 2.5 .* time .^ 2))
    signal = component1 + component2
    result = wsst(
        signal;
        fs,
        freqrange = (20.0, 125.0),
        nfrequencies = 421,
        threshold = RelativeThreshold(1e-9),
        preserve_residual = false,
    )
    ridges = extract_ridges(
        result;
        freqranges = [(30.0, 60.0), (90.0, 115.0)],
        penalty = 1.0,
        continuity_scale = 1.0,
        activity_threshold = 1e-8,
    )
    recovered = reconstruct_components(
        result,
        ridges;
        bandwidth_hz = 8.0,
        ownership = :exclusive,
    )
    truths = (frequency1, frequency2)
    components = (component1, component2)
    for index in eachindex(ridges)
        println("mode $index ridge metrics: ", ridge_metrics(ridges[index], truths[index]))
        println("mode $index component metrics: ",
                component_metrics(recovered.components[index], components[index]))
    end
    reconstructed = sum(recovered.components)
    interior = 129:896
    full_metrics = component_metrics(reconstructed, signal)
    interior_metrics = component_metrics(reconstructed[interior], signal[interior])
    full_similarity = normalized_similarity(signal, reconstructed)
    interior_similarity = normalized_similarity(signal[interior], reconstructed[interior])
    represented = inverse_wsst(result)
    conservation_error = norm(sum(recovered.components) + recovered.residual - represented) /
                         norm(represented)
    println("summed reconstruction full-record metrics: ", full_metrics)
    println("summed reconstruction full-record normalized similarity: ", full_similarity)
    println("summed reconstruction interior metrics: ", interior_metrics)
    println("summed reconstruction interior normalized similarity: ", interior_similarity)
    println("component + SST-residual relative conservation error: ", conservation_error)
    plot_results(
        output_path,
        time,
        signal,
        result,
        ridges,
        truths,
        reconstructed,
        (time[first(interior)], time[last(interior)]),
        full_metrics,
        interior_metrics,
        full_similarity,
        interior_similarity,
    )
    println("saved ", output_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_path = isempty(ARGS) ?
        joinpath(@__DIR__, "ridge_components.png") : abspath(first(ARGS))
    main(output_path)
end
