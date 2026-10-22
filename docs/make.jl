using Documenter
using SynchroSqueezing

makedocs(
    modules = [SynchroSqueezing],
    sitename = "SynchroSqueezing.jl",
    remotes = nothing,
    checkdocs = :exports,
    format = Documenter.HTML(repolink = "https://github.com/svretina/SynchroSqueezing.jl"),
    pages = [
        "Home" => "index.md",
        "Conventions" => "conventions.md",
        "Dependency audit" => "dependency_audit.md",
        "Backend decision" => "backend_decision.md",
        "CWT backend adapter" => "backend_adapter.md",
        "WSST theory" => "wsst_theory.md",
        "Second-order WSST" => "wsst2.md",
        "WSST2 validation" => "wsst2_validation.md",
        "WSST conventions" => "wsst_conventions.md",
        "WSST inverse" => "wsst_inverse.md",
        "Reconstruction audit" => "reconstruction_audit.md",
        "Synthesis theory" => "wsst_synthesis_theory.md",
        "Backend validation" => "wsst_backend_validation.md",
        "Component reconstruction" => "wsst_component_reconstruction.md",
        "Ridge tracking" => "ridge_tracking.md",
        "Predictive ridge tracking" => "predictive_ridge_tracking.md",
        "Damping estimation" => "damping_estimation.md",
        "Performance and regression" => "performance.md",
        "Julia and Python performance" => "julia_python_performance.md",
        "WSST numerics" => "wsst_numerics.md",
        "WSST validation" => "wsst_validation.md",
        "API status" => "api.md",
    ],
)
