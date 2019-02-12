
mutable struct Particles{T<: AbstractFloat, J<: Integer} 
    x::Vector{Vector{T}} # particle positions
    v::Vector{Vector{T}} # particle velocities
    p::Vector{Vector{T}} # best particle positions
    g::Vector{Vector{T}} # best neighborhood swarm position

    lb::Vector{T} #lower bound
    ub::Vector{T} #uper bound

    swarmsize::J #swarm size
    particlelength::J # particle length
 
    fx::Vector{T} # current particle function values
    fs::Vector{Bool} # feasibility of each particle
    fp::Vector{T} # best particle function values
    fg::Vector{T} # best neighborhood swarm position value

    order::Vector{J} # neighborhood order
    best_position::Vector{T} # best position
    best_value::T # best position value
    history::Vector{T} # best position value history
    
    function Particles(swarmsize::V, lb::J, ub::J) where {V<:Integer, J<:Vector}
        @assert length(ub) == length(lb) "ub and lb must have same dimension"
        @assert all(ub .> lb) "Each value of 'ub' must be greater than 'lb'"

        S = length(lb)
        x = [ rand(S) for i = 1:swarmsize]
        v = [ rand(S) for i = 1:swarmsize]

        vhigh = abs.(ub .- lb)
        vlow = -vhigh
        x = map(x-> lb .+ x .* (ub .- lb), x) # particle positions
        v = map(x-> vlow .+ x .* (vhigh .- vlow), v) # particle velocities
        p = copy(x)  # best particle positions
        g = copy(x)  # best neighborhood swarm position

        new{typeof(rand()), typeof(swarmsize) }(x, v, p, g, lb, ub, swarmsize, S, fill(Inf,swarmsize),
         fill(false,swarmsize), fill(Inf,swarmsize), fill(Inf,swarmsize), collect(1:swarmsize),
          [Inf for i=1:S], Inf, [] );
    end
end


