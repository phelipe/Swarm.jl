
function make_constraints(::Type{Val{nothing}}, args, kwargs, verbose)
    verbose && println("No constraints given.")
    return x -> [0.0]
end

function make_constraints(eqs::Vector, args, kwargs, verbose)
    verbose && println("Converting ieqcons to a single constraint function.")
    return x -> [f(x, args...; kwargs...) for f in eqs]
end

function make_constraints(eqs::Function, args, kwargs, verbose)
    verbose && println("Single constraint function given in f_ieqcons.")
    return x -> eqs(x, args...; kwargs...)
end

make_constraints(eqs, args, kwargs, verbose) = make_constraints(Val{eqs}, args, kwargs, verbose)

# atualiza melhores informações das partículas
function update_position!(particle::Particles)
    i_update = (particle.fx .< particle.fp) .& particle.fs
    particle.p[i_update] = copy(particle.x[i_update])
    particle.fp[i_update] = particle.fx[i_update]
end

# funçõa auxiliar para ajuda na busca local
verify(x, size) = (abs(x) % size == 0) ? size : mod(x,size)

# função que atualiza o valor dos melhores vizinhos com base nos dados 
# dos n vizinhos mais próximos de cada partícula
function findBestLocal!(particle::Particles, n::Integer) 
    sortperm!(particle.order, particle.fx, rev=true)
    for i in particle.order
        temp = verify.((i-n):(i+n), particle.swarmsize)
        val, ind = findmin(particle.fp[particle.order[temp]])
        particle.fg[i] = val
        particle.g[i] = particle.p[particle.order[temp]][ind]
    end    
end

# função que atualiza os valores dos melhores vizinhos com base na melhor partícula
function findBestGlobal!(particle::Particles)
    val, ind = findmin(particle.fp)
    particle.fg = fill(val, particle.swarmsize)
    particle.g = fill(particle.p[ind], particle.swarmsize)
end    