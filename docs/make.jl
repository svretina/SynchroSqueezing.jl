using Documenter
using SynchroSqueezing

makedocs(
    modules = [SynchroSqueezing],
    sitename = "SynchroSqueezing.jl",
    remotes = nothing,
    pages = [
        "Home" => "index.md",
        "Conventions" => "conventions.md",
        "Dependency audit" => "dependency_audit.md",
        "Backend decision" => "backend_decision.md",
        "API status" => "api.md",
    ],
)
