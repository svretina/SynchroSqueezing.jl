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
        "API status" => "api.md",
    ],
)
