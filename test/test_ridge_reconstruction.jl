using LinearAlgebra: norm

@testset "Ridge-band component reconstruction and ownership" begin
    fs = 512.0
    n = 512
    time = sample_times(n, fs)
    component1 = cos.(2π * 40.0 .* time)
    component2 = 0.6 .* cos.(2π * 88.0 .* time)
    signal = component1 + component2
    result = wsst(
        signal;
        fs,
        freqrange = (20.0, 110.0),
        nfrequencies = 361,
        threshold = RelativeThreshold(1e-10),
        preserve_residual = false,
    )
    ridges = extract_ridges(
        result;
        freqranges = [(34.0, 46.0), (82.0, 94.0)],
        continuity_scale = 1.0,
    )
    reconstruction1 = reconstruct_component(result, ridges[1]; bandwidth_hz = 4.0)
    reconstruction2 = reconstruct_component(result, ridges[2]; bandwidth_bins = 16)
    @test size(reconstruction1.mask) == size(result.coefficients)
    @test reconstruction1.weights == Float64.(reconstruction1.mask)
    @test norm(reconstruction1.signal - component1) / norm(component1) < 0.02
    @test norm(reconstruction2.signal - component2) / norm(component2) < 0.02

    exclusive = reconstruct_components(
        result,
        ridges;
        bandwidth_hz = [30.0, 30.0],
        ownership = :exclusive,
    )
    fractional = reconstruct_components(
        result,
        ridges;
        bandwidth_hz = 30.0,
        ownership = :fractional,
    )
    for collection in (exclusive, fractional)
        ownership_sum = copy(collection.residual_weights)
        for weights in collection.weights
            ownership_sum .+= weights
        end
        @test maximum(abs, ownership_sum .- 1) < 1e-14
        represented = inverse_wsst(result)
        @test norm(sum(collection.components) + collection.residual - represented) /
              norm(represented) < 1e-12
        @test all(weights -> all(>=(0), weights), collection.weights)
        @test all(>=(0), collection.residual_weights)
    end

    custom = falses(size(result.coefficients))
    custom[1:10, :] .= true
    custom_reconstruction = reconstruct_component(result, ridges[1]; custom_mask = custom)
    @test custom_reconstruction.mask == custom
    @test_throws ArgumentError ridge_mask(
        result,
        ridges[1];
        bandwidth_hz = 2.0,
        bandwidth_bins = 2,
    )
    @test_throws ArgumentError reconstruct_components(
        result,
        ridges;
        bandwidth_hz = 4.0,
        ownership = :invalid,
    )
end
