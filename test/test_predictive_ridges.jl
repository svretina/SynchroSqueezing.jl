using Random: MersenneTwister, randn

function exhaustive_predictive(
    coefficients,
    frequencies,
    times,
    freqrange;
    frequency_penalty,
    frequency_scale,
    curvature_penalty,
    curvature_scale,
    max_rate = nothing,
    max_chirp_rate_change = nothing,
    prediction_bandwidth = nothing,
)
    eligible = findall(frequency -> freqrange[1] <= frequency <= freqrange[2], frequencies)
    references = [maximum(abs2, @view coefficients[eligible, column])
                  for column in axes(coefficients, 2)]
    ranges = ntuple(_ -> eligible, size(coefficients, 2))
    best_path = nothing
    best_cost = Inf
    for candidate in Iterators.product(ranges...)
        cost = 0.0
        feasible = true
        for column in eachindex(candidate)
            reference = references[column]
            cost += iszero(reference) ? 0.0 :
                -log(abs2(coefficients[candidate[column], column]) / reference + eps(Float64))
            column == 1 && continue
            delta_time = times[column] - times[column - 1]
            frequency_change = frequencies[candidate[column]] -
                               frequencies[candidate[column - 1]]
            if max_rate !== nothing && abs(frequency_change) > max_rate * delta_time
                feasible = false
                break
            end
            cost += frequency_penalty * (frequency_change / frequency_scale)^2
            column == 2 && continue
            previous_step = times[column - 1] - times[column - 2]
            previous_slope = (
                frequencies[candidate[column - 1]] - frequencies[candidate[column - 2]]
            ) / previous_step
            next_slope = frequency_change / delta_time
            slope_change = next_slope - previous_slope
            if max_chirp_rate_change !== nothing &&
               abs(slope_change) > max_chirp_rate_change
                feasible = false
                break
            end
            prediction = frequencies[candidate[column - 1]] + previous_slope * delta_time
            if prediction_bandwidth !== nothing &&
               abs(frequencies[candidate[column]] - prediction) > prediction_bandwidth
                feasible = false
                break
            end
            cost += curvature_penalty * (slope_change / curvature_scale)^2
        end
        feasible || continue
        if cost < best_cost
            best_path = collect(candidate)
            best_cost = cost
        end
    end
    return best_path, best_cost
end

@testset "Exact predictive DP matches exhaustive enumeration" begin
    frequencies = [10.0, 18.0, 31.0]
    times = [0.0, 0.4, 1.1, 2.0]
    random = MersenneTwister(4402)
    settings = (
        frequency_penalty = 0.3,
        frequency_scale = 8.0,
        curvature_penalty = 0.7,
        curvature_scale = 20.0,
        max_rate = 40.0,
        max_chirp_rate_change = 50.0,
    )
    for _ in 1:12
        coefficients = complex.(randn(random, 3, 4), randn(random, 3, 4))
        exact_path, exact_cost = exhaustive_predictive(
            coefficients,
            frequencies,
            times,
            (10.0, 31.0);
            settings...,
        )
        estimate = extract_predictive_ridge(
            coefficients,
            frequencies;
            times,
            freqrange = (10.0, 31.0),
            method = :exact,
            settings...,
        )
        @test estimate.indices == exact_path
        @test estimate.total_cost ≈ exact_cost
    end

    tied = extract_predictive_ridge(
        ones(ComplexF64, 3, 4),
        frequencies;
        times,
        freqrange = (10.0, 31.0),
        method = :exact,
        frequency_penalty = 0.0,
        curvature_penalty = 0.0,
    )
    @test tied.indices == fill(1, 4)

    coefficients = complex.(randn(random, 3, 4), randn(random, 3, 4))
    pruned_path, pruned_cost = exhaustive_predictive(
        coefficients,
        frequencies,
        times,
        (10.0, 31.0);
        frequency_penalty = 0.2,
        frequency_scale = 8.0,
        curvature_penalty = 0.5,
        curvature_scale = 20.0,
        prediction_bandwidth = 10.0,
    )
    pruned = extract_predictive_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (10.0, 31.0),
        method = :pruned,
        frequency_penalty = 0.2,
        frequency_scale = 8.0,
        curvature_penalty = 0.5,
        curvature_scale = 20.0,
        prediction_bandwidth = 10.0,
    )
    @test pruned.indices == pruned_path
    @test pruned.total_cost ≈ pruned_cost
end

@testset "Physical constant-velocity prediction" begin
    uniform_times = [0.0, 1.0, 2.0, 3.0]
    uniform_path = [10.0, 20.0, 30.0, 40.0]
    uniform_slopes = diff(uniform_path) ./ diff(uniform_times)
    @test diff(uniform_slopes) == [0.0, 0.0]

    nonuniform_times = [0.0, 1.0, 3.0, 4.0]
    nonuniform_path = 10.0 .+ 5.0 .* nonuniform_times
    nonuniform_slopes = diff(nonuniform_path) ./ diff(nonuniform_times)
    @test diff(nonuniform_slopes) == [0.0, 0.0]

    frequencies = [10.0, 15.0, 20.0, 25.0, 30.0]
    coefficients = fill(0.01 + 0im, length(frequencies), length(nonuniform_times))
    truth_indices = [1, 2, 4, 5]
    for column in eachindex(truth_indices)
        coefficients[truth_indices[column], column] = 1
    end
    ridge = extract_predictive_ridge(
        coefficients,
        frequencies;
        times = nonuniform_times,
        freqrange = (10.0, 30.0),
        curvature_penalty = 5.0,
        curvature_scale = 1.0,
        frequency_penalty = 0.0,
    )
    @test ridge.indices == truth_indices
end

@testset "Predictive strategies on competing smooth and erratic ridges" begin
    frequencies = collect(10.0:10.0:50.0)
    times = collect(0.0:4.0)
    truth_indices = [1, 2, 3, 4, 5]
    competitor_indices = [1, 3, 2, 4, 3]
    coefficients = fill(0.05 + 0im, 5, 5)
    for column in eachindex(times)
        coefficients[truth_indices[column], column] = 3
        coefficients[competitor_indices[column], column] = 4
    end
    first_order = extract_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (10.0, 50.0),
        method = :dense,
        penalty = 0.0,
    )
    exact = extract_predictive_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (10.0, 50.0),
        method = :exact,
        frequency_penalty = 0.0,
        curvature_penalty = 2.0,
        curvature_scale = 10.0,
    )
    pruned = extract_predictive_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (10.0, 50.0),
        method = :pruned,
        prediction_bandwidth = 5.0,
        frequency_penalty = 0.0,
        curvature_penalty = 2.0,
        curvature_scale = 10.0,
    )
    truth = frequencies[truth_indices]
    @test ridge_metrics(exact, truth).rmse < ridge_metrics(first_order, truth).rmse
    @test exact.indices == truth_indices
    @test pruned.indices == truth_indices
    @test exact.config.method == :exact
    @test pruned.config.method == :pruned
end

@testset "Predictive constraints, priors, and fallback" begin
    frequencies = [10.0, 20.0, 30.0, 40.0]
    times = [0.0, 1.0, 2.0, 3.0]
    abrupt_indices = [1, 1, 4, 4]
    coefficients = fill(0.01 + 0im, 4, 4)
    for column in eachindex(abrupt_indices)
        coefficients[abrupt_indices[column], column] = 10
    end
    abrupt = extract_predictive_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (10.0, 40.0),
        frequency_penalty = 0.0,
        curvature_penalty = 0.0,
    )
    @test abrupt.indices == abrupt_indices
    narrow = extract_predictive_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (10.0, 40.0),
        method = :pruned,
        prediction_bandwidth = 1.0,
        frequency_penalty = 0.0,
        curvature_penalty = 0.0,
    )
    wide = extract_predictive_ridge(
        coefficients,
        frequencies;
        times,
        freqrange = (10.0, 40.0),
        method = :pruned,
        prediction_bandwidth = 30.0,
        frequency_penalty = 0.0,
        curvature_penalty = 0.0,
    )
    @test narrow.indices != abrupt_indices
    @test wide.indices == abrupt_indices

    restricted = extract_predictive_ridge(
        [100 100 100 100; coefficients[2:3, :]; 100 100 100 100],
        frequencies;
        times,
        freqrange = (20.0, 30.0),
        frequency_penalty = 0.0,
        curvature_penalty = 0.0,
        initial_frequency = 1000.0,
        initial_frequency_penalty = 10.0,
    )
    @test all((20.0 .<= restricted.frequencies) .& (restricted.frequencies .<= 30.0))

    short = extract_predictive_ridge(
        ones(ComplexF64, 4, 2),
        frequencies;
        times = [0.0, 1.0],
        freqrange = (10.0, 40.0),
    )
    @test short.config isa RidgeConfig
    @test short.indices == [1, 1]

    @test_throws ArgumentError PredictiveRidgeConfig(method = :pruned)
    @test_throws ArgumentError PredictiveRidgeConfig(
        method = :exact,
        prediction_bandwidth = 2.0,
    )
end
