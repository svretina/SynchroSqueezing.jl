"""Compute frequency error and coverage without treating inactive samples as measurements."""
function ridge_metrics(
    ridge::RidgeResult,
    truth::AbstractVector;
    expected = trues(length(truth)),
)
    length(truth) == length(ridge.frequencies) || throw(DimensionMismatch(
        "ridge and truth must have the same number of samples",
    ))
    length(expected) == length(truth) || throw(DimensionMismatch(
        "expected mask must match the ridge length",
    ))
    measured = ridge.valid .& ridge.active .& expected
    expected_count = count(expected)
    expected_count > 0 || throw(ArgumentError("expected mask must select at least one sample"))
    count(measured) > 0 || return (
        rmse = NaN,
        bias = NaN,
        maximum_error = NaN,
        coverage = 0.0,
        boundary_fraction = 0.0,
    )
    errors = ridge.frequencies[measured] .- truth[measured]
    return (
        rmse = sqrt(sum(abs2, errors) / length(errors)),
        bias = sum(errors) / length(errors),
        maximum_error = maximum(abs, errors),
        coverage = count(measured) / expected_count,
        boundary_fraction = count(ridge.at_frequency_boundary .& measured) / count(measured),
    )
end

"""Return primary relative waveform error and signal-to-reconstruction-error ratio."""
function component_metrics(reconstructed::AbstractVector, truth::AbstractVector)
    length(reconstructed) == length(truth) || throw(DimensionMismatch(
        "reconstructed and true components must have the same length",
    ))
    signal_power = sum(abs2, truth)
    signal_power > 0 || throw(ArgumentError("component truth must have nonzero norm"))
    error_power = sum(abs2, reconstructed .- truth)
    return (
        relative_error = sqrt(error_power / signal_power),
        srer_db = iszero(error_power) ? Inf : 10log10(signal_power / error_power),
    )
end
