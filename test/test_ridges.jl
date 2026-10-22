using LinearAlgebra: norm
using Random: MersenneTwister, randn

function exhaustive_ridge(coefficients, frequencies, times, freqrange;
                          penalty, continuity_scale, max_rate = nothing)
    eligible = findall(frequency -> freqrange[1] <= frequency <= freqrange[2], frequencies)
    references = [maximum(abs2, @view coefficients[eligible, column])
                  for column in axes(coefficients, 2)]
    epsilon = eps(Float64)
    best_path = nothing
    best_cost = Inf
    ranges = ntuple(_ -> eligible, size(coefficients, 2))
    for candidate in Iterators.product(ranges...)
        feasible = all(2:length(candidate)) do column
            max_rate === nothing ||
                abs(frequencies[candidate[column]] - frequencies[candidate[column - 1]]) <=
                max_rate * (times[column] - times[column - 1])
        end
        feasible || continue
        cost = 0.0
        for column in eachindex(candidate)
            reference = references[column]
            cost += iszero(reference) ? 0.0 :
                -log(abs2(coefficients[candidate[column], column]) / reference + epsilon)
            if column > 1
                change = frequencies[candidate[column]] - frequencies[candidate[column - 1]]
                cost += penalty * (change / continuity_scale)^2
            end
        end
        if cost < best_cost
            best_cost = cost
            best_path = collect(candidate)
        end
    end
    return best_path, best_cost
end

@testset "Ridge continuity cache is memory bounded" begin
    config = RidgeConfig(Float64)
    small = SynchroSqueezing._cache_transition_costs(
        collect(1.0:16.0), 1, 16, config,
    )
    @test size(small) == (16, 16)
    @test small[3, 7] == SynchroSqueezing._transition_cost(3.0, 7.0, config)

    large_frequencies = collect(1.0:3000.0)
    @test SynchroSqueezing._cache_transition_costs(
        large_frequencies, 1, length(large_frequencies), config,
    ) === nothing
end

@testset "Frequency-constrained ridge baseline" begin
    frequencies = [10.0, 20.0, 30.0, 40.0]
    times = [0.0, 1.0, 2.0]
    coefficients = ComplexF64[
        1 1 1
        3 3 3
        2 2 2
        100 100 100
    ]
    ridge = extract_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (15.0, 25.0),
        method = :maximum,
    )
    @test ridge.indices == [2, 2, 2]
    @test ridge.frequencies == fill(20.0, 3)
    @test all(ridge.valid)
    @test all(ridge.active)
    @test all(ridge.at_frequency_boundary)

    tied = extract_ridge(
        ones(ComplexF64, 4, 3),
        frequencies;
        times,
        freqrange = (20.0, 30.0),
        method = :maximum,
    )
    @test tied.indices == fill(2, 3)

    zeros_result = extract_ridge(
        zeros(ComplexF64, 4, 3),
        frequencies;
        times,
        freqrange = (20.0, 30.0),
        method = :maximum,
    )
    @test all(zeros_result.valid)
    @test !any(zeros_result.active)

    masked = trues(4, 3)
    masked[2, 2] = false
    masked[3, 2] = false
    invalid = extract_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (20.0, 30.0),
        method = :maximum,
        valid_mask = masked,
    )
    @test invalid.indices[2] == 0
    @test !invalid.valid[2]
    @test isnan(invalid.frequencies[2])

    @test_throws ArgumentError extract_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (21.0, 29.0),
    )
    @test_throws ArgumentError extract_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (30.0, 20.0),
    )
end

@testset "Dense DP is globally optimal and deterministic" begin
    frequencies = [10.0, 18.0, 31.0]
    times = [0.0, 0.5, 1.25, 2.0]
    coefficients = ComplexF64[
        4 1 4 4
        2 5 2 2
        1 2 1 1
    ]
    settings = (
        penalty = 2.0,
        continuity_scale = 8.0,
        max_rate = 20.0,
    )
    exact_path, exact_cost = exhaustive_ridge(
        coefficients,
        frequencies,
        times,
        (10.0, 31.0);
        settings...,
    )
    dense = extract_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (10.0, 31.0),
        method = :dense,
        settings...,
    )
    banded = extract_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (10.0, 31.0),
        method = :dp,
        settings...,
    )
    @test dense.indices == exact_path
    @test dense.total_cost ≈ exact_cost
    @test banded.indices == dense.indices
    @test banded.total_cost == dense.total_cost

    ties = extract_ridge(
        ones(ComplexF64, 3, 4),
        frequencies;
        times,
        freqrange = (10.0, 31.0),
        method = :dense,
        penalty = 0.0,
    )
    @test ties.indices == fill(1, 4)

    jump = ComplexF64[10 0; 0 10]
    constrained = extract_ridge(
        jump,
        [10.0, 30.0];
        times = [0.0, 1.0],
        freqrange = (10.0, 30.0),
        method = :dp,
        penalty = 0.0,
        max_rate = 5.0,
    )
    @test constrained.indices in ([1, 1], [2, 2])

    envelope = [(10.0, 18.0), (18.0, 31.0), (10.0, 18.0), (10.0, 18.0)]
    enveloped = extract_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (10.0, 31.0),
        frequency_envelope = envelope,
        method = :dense,
        penalty = 0.0,
    )
    @test enveloped.frequencies[1] <= 18.0
    @test 18.0 <= enveloped.frequencies[2] <= 31.0

    random = MersenneTwister(2404)
    for _ in 1:12
        random_coefficients = complex.(randn(random, 3, 4), randn(random, 3, 4))
        exact_path, exact_cost = exhaustive_ridge(
            random_coefficients,
            frequencies,
            times,
            (10.0, 31.0);
            penalty = 0.7,
            continuity_scale = 9.0,
            max_rate = 30.0,
        )
        estimate = extract_ridge(
            random_coefficients,
            frequencies;
            times,
            freqrange = (10.0, 31.0),
            method = :dense,
            penalty = 0.7,
            continuity_scale = 9.0,
            max_rate = 30.0,
        )
        @test estimate.indices == exact_path
        @test estimate.total_cost ≈ exact_cost
    end
end

@testset "Independent constrained ridges preserve range order" begin
    fs = 256.0
    n = 256
    time = sample_times(n, fs)
    signal = cos.(2π * 30.0 .* time) .+ 0.5cos.(2π * 80.0 .* time)
    result = wsst(
        signal;
        fs,
        freqrange = (15.0, 100.0),
        nfrequencies = 341,
        threshold = RelativeThreshold(1e-10),
        preserve_residual = false,
    )
    ridges = extract_ridges(
        result;
        freqranges = [(70.0, 90.0), (20.0, 40.0)],
        penalty = 1.0,
        continuity_scale = 1.0,
    )
    @test maximum(abs.(ridges[1].frequencies .- 80.0)) <= 0.25
    @test maximum(abs.(ridges[2].frequencies .- 30.0)) <= 0.25
    @test !any(ridge_overlap(ridges[1], ridges[2]))
end

@testset "Ridge activity and metrics" begin
    frequencies = [10.0, 20.0]
    coefficients = ComplexF64[1 0 1; 0 0 0]
    ridge = extract_ridge(
        coefficients,
        frequencies;
        times = [0.0, 1.0, 2.0],
        freqrange = (10.0, 20.0),
        method = :dense,
        activity_threshold = 0.1,
    )
    @test ridge.active == BitVector([true, false, true])
    metrics = ridge_metrics(ridge, fill(10.0, 3))
    @test metrics.rmse == 0
    @test metrics.coverage == 2 / 3
    @test component_metrics([1.0, 2.0], [1.0, 2.0]).relative_error == 0
    @test component_metrics([1.0, 2.0], [1.0, 2.0]).srer_db == Inf
end

@testset "Activity, range boundaries, and refinement" begin
    frequencies = [10.0, 20.0, 30.0, 40.0]
    times = [0.0, 1.0, 2.0]
    noise_interval = ComplexF64[
        100 100 100
        0.01 0.01 0.01
        0.005 0.005 0.005
        0 0 0
    ]
    noise_ridge = extract_ridge(
        noise_interval,
        frequencies;
        times,
        freqrange = (20.0, 30.0),
        activity_threshold = 1e-4,
    )
    @test !any(noise_ridge.active)

    leaving = ComplexF64[
        8 1 1
        1 8 1
        0 1 8
        0 0 1
    ]
    clipped = extract_ridge(
        leaving,
        frequencies;
        times,
        freqrange = (10.0, 30.0),
        method = :maximum,
    )
    @test clipped.frequencies == [10.0, 20.0, 30.0]
    @test clipped.at_frequency_boundary[end]

    refinement_errors = Float64[]
    for bins in (160, 320, 640)
        grid = collect(range(20.0, 80.0; length = bins))
        local_times = collect(range(0.0, 1.0; length = 80))
        truth = 36.0 .+ 18.0 .* local_times
        coefficients = ComplexF64[
            exp(-abs2(frequency - truth[column]) / 2)
            for frequency in grid, column in eachindex(local_times)
        ]
        ridge = extract_ridge(
            coefficients,
            grid;
            times = local_times,
            freqrange = (32.0, 58.0),
            method = :dp,
            penalty = 0.2,
            continuity_scale = 1.0,
        )
        push!(refinement_errors, ridge_metrics(ridge, truth).rmse)
    end
    @test refinement_errors[3] < refinement_errors[1]
    @test maximum(refinement_errors) < 0.2
end

@testset "Typed and per-range ridge configurations" begin
    frequencies = Float32[10, 20, 40]
    times = Float32[0, 1, 2]
    coefficients = ComplexF32[1 1 1; 4 4 4; 1 1 1]
    config = RidgeConfig(
        Float32;
        method = :dense,
        continuity = :log,
        penalty = 0.5,
        continuity_scale = 0.2,
    )
    ridge = extract_ridge(
        coefficients,
        frequencies,
        config;
        times,
        freqrange = (10.0f0, 40.0f0),
    )
    @test eltype(ridge.frequencies) == Float32
    @test ridge.indices == fill(2, 3)

    fs = 128.0
    n = 128
    time = sample_times(n, fs)
    result = wsst(
        cos.(2π * 20.0 .* time) .+ 0.5cos.(2π * 45.0 .* time);
        fs,
        freqrange = (10.0, 55.0),
        nfrequencies = 181,
        preserve_residual = false,
    )
    configs = [
        (method = :maximum,),
        RidgeConfig(Float64; method = :dense, penalty = 2.0),
    ]
    ridges = extract_ridges(
        result;
        freqranges = [(15.0, 25.0), (40.0, 50.0)],
        configs,
    )
    @test all((15.0 .<= ridges[1].frequencies) .& (ridges[1].frequencies .<= 25.0))
    @test all((40.0 .<= ridges[2].frequencies) .& (ridges[2].frequencies .<= 50.0))
end
