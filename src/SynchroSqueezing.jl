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
include("reconstruction/inverse_wsst.jl")

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
       WSSTDiagnostics,
       WSSTResult,
       WaveletValidation,
       admissibility_constant,
       cwt,
       instantaneous_frequency,
       inverse_cwt,
       inverse_wsst,
       linear_frequency_grid,
       log_frequency_grid,
       package_version,
       plan_cwt,
       reconstruct_region,
       synchrosqueeze!,
       synchrosqueeze_reference!,
       wsst

end # module SynchroSqueezing
