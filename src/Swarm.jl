module Swarm

include("particles.jl")
include("functions.jl")
include("PSO.jl")
include("HPSO.jl")
include("IHPSO.jl")

export Particles
export pso, hpso, ihpso

end # module
