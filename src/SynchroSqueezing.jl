module SynchroSqueezing

using FFTW
using SpecialFunctions: gamma
import ContinuousWavelets

include("types.jl")
include("backends/interface.jl")
include("backends/continuous_wavelets.jl")
include("squeezing/grids.jl")
include("squeezing/phase.jl")
include("squeezing/reassignment.jl")
include("squeezing/wsst.jl")
include("ridges/types.jl")
include("ridges/tracking.jl")
include("ridges/predictive.jl")
include("reconstruction/inverse_wsst.jl")
include("ridges/reconstruction.jl")
include("ridges/metrics.jl")
include("analysis/types.jl")
include("analysis/damping.jl")
include("analysis/mode_support.jl")

"""
    package_version()

Return the package version declared for the current development checkout.

This small API exists only to verify the package skeleton during Milestone 00.
Transform, ridge, and reconstruction APIs will be added after the backend audit
and experimental validation are complete.
"""
package_version() = v"0.1.0-DEV"

export AbstractCWTBackend,
       AbstractMagnitudeThreshold,
       AbsoluteThreshold,
       ContinuousWaveletsBackend,
       CWTPlan,
       CWTMetadata,
       CWTResult,
       FrequencyEstimate,
       NoThreshold,
       ReconstructionConvention,
       RelativeThreshold,
       RidgeConfig,
       PredictiveRidgeConfig,
       RidgeResult,
       ComponentReconstruction,
       ComponentSet,
       ComponentFeatures,
       DampingFit,
       DampingSensitivityCase,
       DampingSensitivity,
       RidgeProminence,
       ModeSupportConfig,
       ModeSupportResult,
       SupportedComponentReconstruction,
       ModeSupportSensitivityCase,
       ModeSupportSensitivity,
       WSSTDiagnostics,
       WSSTResult,
       WaveletValidation,
       admissibility_constant,
       cwt,
       component_metrics,
       component_features,
       analytic_signal,
       fit_exponential_damping,
       local_damping_rate,
       damping_sensitivity,
       ridge_prominence,
       mode_support,
       supported_reconstruction,
       mode_support_sensitivity,
       extract_ridge,
       extract_ridges,
       extract_predictive_ridge,
       instantaneous_frequency,
       inverse_cwt,
       inverse_wsst,
       linear_frequency_grid,
       log_frequency_grid,
       package_version,
       plan_cwt,
       reconstruct_region,
       reconstruct_component,
       reconstruct_components,
       ridge_mask,
       ridge_metrics,
       ridge_overlap,
       synchrosqueeze!,
       synchrosqueeze_reference!,
       wsst

end # module SynchroSqueezing
