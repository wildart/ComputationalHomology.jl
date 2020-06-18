using Documenter, ComputationalHomology

makedocs(
    modules = [ComputationalHomology],
    doctest = false,
    clean = true,
    sitename = "ComputationalHomology.jl",
    pages = [
        "Home" => "index.md",
        #"Tutorial" => "tutorial.md",
        "Chains" => "chains.md",
    ]
)

deploydocs(repo = "github.com/wildart/ComputationalHomology.jl.git")
