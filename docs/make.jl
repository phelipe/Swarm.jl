using Documenter, Swarm

makedocs(
    modules = [Swarm],
    format = :html,
    checkdocs = :exports,
    sitename = "Swarm.jl",
    pages = Any["index.md"]
)

deploydocs(
    repo = "github.com/phelipe/Swarm.jl.git",
)
