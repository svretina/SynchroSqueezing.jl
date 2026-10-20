using CairoMakie
using DelimitedFiles: readdlm
using SynchroSqueezing
using TOML

const REPOSITORY_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_CONFIG = joinpath(@__DIR__, "three_oscillators.toml")
const PYTHON_SCRIPT = joinpath(
    REPOSITORY_ROOT,
    "reference",
    "python",
    "compare_three_oscillators.py",
)

function three_oscillator_signal(config)
    samples = Int(config["samples"])
    sample_rate = config["sample_rate_hz"]
    time = collect(0:samples-1) ./ sample_rate

    oscillator1 = config["oscillator_1"]
    frequency1 = fill(oscillator1["frequency_hz"], samples)
    envelope1 = oscillator1["amplitude"] .* (
        1 .+ oscillator1["modulation_depth"] .*
        cos.(2π * oscillator1["modulation_frequency_hz"] .* time)
    )
    component1 = envelope1 .* cos.(2π * oscillator1["frequency_hz"] .* time)

    oscillator2 = config["oscillator_2"]
    rate2 = oscillator2["chirp_rate_hz_per_second"]
    frequency2 = oscillator2["initial_frequency_hz"] .+ rate2 .* time
    component2 = oscillator2["amplitude"] .* cos.(2π .* (
        oscillator2["initial_frequency_hz"] .* time .+ 0.5 .* rate2 .* time .^ 2
    ))

    oscillator3 = config["oscillator_3"]
    rate3 = oscillator3["chirp_rate_hz_per_second"]
    frequency3 = oscillator3["initial_frequency_hz"] .+ rate3 .* time
    component3 = oscillator3["amplitude"] .* cos.(2π .* (
        oscillator3["initial_frequency_hz"] .* time .+ 0.5 .* rate3 .* time .^ 2
    ))

    return time, (component1, component2, component3),
           (frequency1, frequency2, frequency3)
end

function python_executable()
    configured = get(ENV, "PYTHON", "")
    !isempty(configured) && return configured

    local_environment = joinpath(
        REPOSITORY_ROOT,
        "reference",
        "python",
        ".venv",
        "bin",
        "python",
    )
    isfile(local_environment) && return local_environment

    system_python = Sys.which("python3")
    system_python === nothing && error(
        "Python was not found; set PYTHON to an environment containing ssqueezepy 0.6.6",
    )
    return system_python
end

function python_wsst(config_path, expected_signal)
    return mktempdir() do output_directory
        run(`$(python_executable()) $PYTHON_SCRIPT --config $config_path --output $output_directory`)
        signal = vec(readdlm(joinpath(output_directory, "signal.csv"), ',', Float64))
        frequencies = vec(readdlm(
            joinpath(output_directory, "frequencies.csv"),
            ',',
            Float64,
        ))
        decibels = Matrix{Float64}(readdlm(
            joinpath(output_directory, "wsst_decibels.csv"),
            ',',
            Float64,
        ))
        maximum(abs, signal - expected_signal) < 1e-12 || error(
            "Julia and Python generated different signals",
        )
        return frequencies, decibels
    end
end

function relative_decibels(coefficients)
    relative = abs.(coefficients) ./ maximum(abs, coefficients)
    return 20 .* log10.(max.(relative, 1e-6))
end

function local_ridge_rmse(decibels, frequencies, exact_frequency; half_width = 5.0)
    estimates = similar(exact_frequency)
    for time_index in eachindex(exact_frequency)
        candidates = findall(
            frequency -> abs(frequency - exact_frequency[time_index]) <= half_width,
            frequencies,
        )
        isempty(candidates) && error("the local ridge window does not intersect the frequency grid")
        local_index = argmax(@view decibels[candidates, time_index])
        estimates[time_index] = frequencies[candidates[local_index]]
    end
    return sqrt(sum(abs2, estimates - exact_frequency) / length(estimates))
end

function plot_comparison(
    output_path,
    time,
    components,
    exact_frequencies,
    julia_result,
    python_frequencies,
    python_decibels,
)
    signal = sum(components)
    julia_decibels = relative_decibels(julia_result.coefficients)
    julia_rmse = [
        local_ridge_rmse(julia_decibels, julia_result.frequencies, frequency)
        for frequency in exact_frequencies
    ]
    python_rmse = [
        local_ridge_rmse(python_decibels, python_frequencies, frequency)
        for frequency in exact_frequencies
    ]

    colors = Makie.wong_colors()[1:3]
    figure = Figure(size = (1800, 1080), fontsize = 19)
    signal_axis = Axis(
        figure[1, 1:2],
        ylabel = "amplitude",
        title = "Identical three-oscillator signal supplied to both implementations",
    )
    lines!(signal_axis, time, signal; color = :black, linewidth = 1.2, label = "sum")
    for (index, component) in enumerate(components)
        lines!(
            signal_axis,
            time,
            component;
            color = (colors[index], 0.65),
            linewidth = 1.0,
            label = "oscillator $index",
        )
    end
    axislegend(signal_axis; position = :rt, orientation = :horizontal, nbanks = 2)
    hidexdecorations!(signal_axis; grid = false)

    julia_axis = Axis(
        figure[2, 1],
        xlabel = "time (s)",
        ylabel = "frequency (Hz)",
        title = "SynchroSqueezing.jl  |  ridge RMSE = " *
                join(string.(round.(julia_rmse; digits = 3)), ", ") * " Hz",
    )
    python_axis = Axis(
        figure[2, 2],
        xlabel = "time (s)",
        title = "ssqueezepy 0.6.6  |  ridge RMSE = " *
                join(string.(round.(python_rmse; digits = 3)), ", ") * " Hz",
    )
    image = heatmap!(
        julia_axis,
        time,
        julia_result.frequencies,
        permutedims(julia_decibels);
        colormap = :inferno,
        colorrange = (-60, 0),
    )
    heatmap!(
        python_axis,
        time,
        python_frequencies,
        permutedims(python_decibels);
        colormap = :inferno,
        colorrange = (-60, 0),
    )
    for (index, frequency) in enumerate(exact_frequencies)
        for axis in (julia_axis, python_axis)
            lines!(
                axis,
                time,
                frequency;
                color = colors[index],
                linewidth = 2.0,
                linestyle = :dash,
            )
        end
    end
    Colorbar(figure[2, 3], image; label = "magnitude relative to panel maximum (dB)")
    truth_elements = [
        LineElement(color = colors[index], linewidth = 2.0, linestyle = :dash)
        for index in eachindex(exact_frequencies)
    ]
    Legend(
        figure[3, 1:2],
        truth_elements,
        ["oscillator $index analytical IF" for index in eachindex(exact_frequencies)];
        orientation = :horizontal,
        tellwidth = false,
        framevisible = false,
        title = "Dashed curves: ground-truth instantaneous frequency",
    )
    linkaxes!(julia_axis, python_axis)
    linkxaxes!(signal_axis, julia_axis)
    rowgap!(figure.layout, 12)
    colgap!(figure.layout, 15)
    save(output_path, figure; px_per_unit = 2.0)

    println("saved ", output_path)
    println("SynchroSqueezing.jl ridge RMSE (Hz): ", julia_rmse)
    println("ssqueezepy ridge RMSE (Hz): ", python_rmse)
    println("reliable Julia frequency range: ",
            julia_result.reconstruction.reliable_frequency_range, " Hz")
end

function main(output_path; config_path = DEFAULT_CONFIG)
    config = TOML.parsefile(config_path)
    time, components, exact_frequencies = three_oscillator_signal(config)
    signal = sum(components)
    wavelet = config["wavelet"]
    cwt_options = (;
        voices_per_octave = Int(config["voices_per_octave"]),
        wavelet = SynchroSqueezing.ContinuousWavelets.Morse(
            wavelet["gamma"],
            wavelet["beta"],
            wavelet["center_frequency"],
        ),
    )
    result = wsst(
        signal;
        fs = config["sample_rate_hz"],
        cwt_options,
        freqrange = (config["frequency_min_hz"], config["frequency_max_hz"]),
        nfrequencies = Int(config["frequency_bins"]),
        reassignment = :nearest,
        threshold = NoThreshold(),
        preserve_residual = false,
    )
    python_frequencies, python_decibels = python_wsst(config_path, signal)
    maximum(abs, python_frequencies - result.frequencies) < 1e-12 || error(
        "Julia and Python used different output frequency grids",
    )
    plot_comparison(
        output_path,
        time,
        components,
        exact_frequencies,
        result,
        python_frequencies,
        python_decibels,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_path = isempty(ARGS) ?
        joinpath(@__DIR__, "wsst_three_oscillators.png") : abspath(first(ARGS))
    main(output_path)
end
