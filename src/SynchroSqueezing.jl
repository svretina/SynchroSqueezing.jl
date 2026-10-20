module SynchroSqueezing

using FFTW
import ContinuousWavelets

include("types.jl")
include("backends/interface.jl")
include("backends/continuous_wavelets.jl")

"""
    package_version()

Return the package version declared for the current development checkout.

This small API exists only to verify the package skeleton during Milestone 00.
Transform, ridge, and reconstruction APIs will be added after the backend audit
and experimental validation are complete.
"""
package_version() = v"0.1.0-DEV"

export AbstractCWTBackend,
       ContinuousWaveletsBackend,
       CWTPlan,
       CWTMetadata,
       CWTResult,
       FrequencyEstimate,
       WaveletValidation,
       cwt,
       instantaneous_frequency,
       inverse_cwt,
       package_version,
       plan_cwt

end # module SynchroSqueezing
