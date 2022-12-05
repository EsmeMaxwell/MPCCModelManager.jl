
# Hacky until the package is registered, then there should be a Project.toml in docs with dependencies listed
push!(LOAD_PATH, "../")

using Documenter, MPCCModelManager

makedocs(
    sitename = "MPCCModelManager.jl",
    authors = "Peter Maxwell",
    clean = true,
    doctest = false,
    format = Documenter.HTML(),
    modules = [MPCCModelManager],
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "manual/function_list.md"
        ]
    ]
)


deploydocs(
    repo = "github.com/PeterMaxwell/MPCCModelManager.jl.git",
)