#!/usr/bin/env julia

using CairoMakie
using Random: MersenneTwister, randn
using SynchroSqueezing

function synthetic_evidence()
    times = collect(range(0.0, 2.0; length = 160))
    frequencies = collect(range(20.0, 100.0; length = 161))
    truth = 34.0 .+ 13.0 .* times .+ 4.0 .* times .^ 2
    competitor = truth .+ 7.0 .* sin.(7π .* times)
    random = MersenneTwister(404)
    coefficients = Matrix{ComplexF64}(undef, length(frequencies), length(times))
    for column in eachindex(times), row in eachindex(frequencies)
        frequency = frequencies[row]
        target = exp(-0.5 * ((frequency - truth[column]) / 0.8)^2)
        distraction = 1.15 * exp(-0.5 * ((frequency - competitor[column]) / 0.8)^2)
        noise = 0.025 * complex(randn(random), randn(random))
        coefficients[row, column] = target + distraction + noise
    end
    return times, frequencies, coefficients, truth
end

function relative_decibels(coefficients)
    relative = abs.(coefficients) ./ maximum(abs, coefficients)
    return 20 .* log10.(max.(relative, 1e-6))
end

function main(output_path)
    times, frequencies, coefficients, truth = synthetic_evidence()
    common = (
        times = times,
        freqrange = (28.0, 78.0),
        max_rate = 200.0,
        activity_threshold = 1e-6,
    )
    first_order = extract_ridge(
        coefficients,
        frequencies;
        common...,
        method = :dp,
        penalty = 0.1,
        continuity_scale = 1.0,
    )
    exact = extract_predictive_ridge(
        coefficients,
        frequencies;
        common...,
        method = :exact,
        frequency_penalty = 0.1,
        frequency_scale = 1.0,
        curvature_penalty = 0.05,
        curvature_scale = 20.0,
    )
    pruned = extract_predictive_ridge(
        coefficients,
        frequencies;
        common...,
        method = :pruned,
        prediction_bandwidth = 1.5,
        frequency_penalty = 0.1,
        frequency_scale = 1.0,
        curvature_penalty = 0.05,
        curvature_scale = 20.0,
    )
    strategies = (
        "first-order DP" => first_order,
        "exact second-order DP" => exact,
        "prediction-window DP" => pruned,
    )
    for (name, ridge) in strategies
        println(name, ": ", ridge_metrics(ridge, truth))
    end

    figure = Figure(size = (1700, 950), fontsize = 20)
    axis = Axis(
        figure[1, 1],
        xlabel = "time (s)",
        ylabel = "frequency (Hz)",
        title = "First- and second-order tracking on identical time-frequency evidence",
    )
    image = heatmap!(
        axis,
        times,
        frequencies,
        permutedims(relative_decibels(coefficients));
        colormap = :inferno,
        colorrange = (-45, 0),
    )
    lines!(axis, times, truth; color = :white, linewidth = 3.0,
           linestyle = :dash, label = "ground truth")
    colors = (:deepskyblue2, :chartreuse3, :orange)
    styles = (:dot, :solid, :dashdot)
    for (index, (name, ridge)) in enumerate(strategies)
        metrics = ridge_metrics(ridge, truth)
        lines!(
            axis,
            times,
            ridge.frequencies;
            color = colors[index],
            linewidth = 2.5,
            linestyle = styles[index],
            label = "$name — RMSE $(round(metrics.rmse; digits = 3)) Hz",
        )
    end
    axislegend(axis; position = :rb)
    Colorbar(figure[1, 2], image; label = "relative magnitude (dB)")
    save(output_path, figure; px_per_unit = 2.0)
    println("saved ", output_path)
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_path = isempty(ARGS) ?
        joinpath(@__DIR__, "predictive_ridge_strategies.png") : abspath(first(ARGS))
    main(output_path)
end
