# ω0 => fator inercial inicial
# ωf => fator inercial final
# ϕp => coeficiente de aceleração partícula
# ϕg => coeficiente de aceleração grupo
# η => passo incremental para criação de nova solução
# λ => decaimento da temperatura

function hpso(func::Function, particles::Particles, constraints, args, kwargs, 
    ω0, ωf , ϕp, ϕg, η, λ , maxiter, maxiterSA, minstep, minfunc, verbose)

    obj = x -> func(x, args...; kwargs...)
    cons = make_constraints(constraints, args, kwargs, verbose)
    is_feasible = x -> all(cons(x) .>= 0)

    S = particles.swarmsize
    D = particles.particlelength 
    verbose && println("Initialize current particle function values.")
    particles.fx = [obj(particles.x[i]) for i = 1:S]  # current particle function values
    particles.fs = [is_feasible(particles.x[i]) for i = 1:S]  # feasibility of each particle
    
    particles.history = []
    particles.fp = copy(particles.fx)  # best particle function values
    particles.order = sortperm(particles.fx, rev=true) #particle neighborhood order

    particles.g = copy(particles.x)  # best neighborhood swarm position
    particles.fg = ones(S) * Inf  # best neighborhood swarm position starting value

    # procurar alguma solução feasible e depois
    # pegar a com o menor valor dentre estas
    if .|(particles.fs...)
        aux1, aux2 = findmin(particles.fx[particles.fs])
        particles.best_position = copy(particles.x[particles.fs][aux2]) # best position
        particles.best_value = aux1 # best position value    
    else
        aux1, aux2 = findmin(particles.fx)
        particles.best_position = copy(particles.x[aux2]) # best position
        particles.best_value = aux1 # best position value
    end


    # Store particle's best position (if constraints are satisfied)
    verbose && println("Store particle's best position.")
    update_position!(particles)
    
    # Update swarm's best position
    verbose && println("Update swarm's best position.")
    findBestGlobal!(particles)
    
    #temperature
    temperature = -(maximum(particles.fp) - minimum(particles.fp) )/ log(0.1)


    it = 1
    it_best = 1
    verbose && println("Start iteration.")
    while it <= maxiter
    
        # Update the particles' velocities and positions using equations 1 and 2
        # Eq(1): v = ω*v + ϕp . rp . (p -x) + ϕg . rg . (g - x)
        # Eq(2): x = x +v          
        ω = ωf + (ω0 - ωf) * ( (maxiter - it) / (maxiter - 1) )
        rp = rand(particles.swarmsize) # random value 
        rg = rand(particles.swarmsize) # random value
        particles.v = (ω .* particles.v) .+ (ϕp .* rp .* (particles.p .- particles.x)) .+
        (ϕg .* rg .* (particles.g .- particles.x)) 
        particles.x = particles.x .+ particles.v    

        # Correct for bound violations
        maskl = [particles.x[i] .< particles.lb for i = 1:S]
        masku = [particles.x[i] .> particles.ub for i = 1:S]
        mask = [ maskl[i] .| masku[i] for i = 1:S]
        particles.x = map((x,y) -> x.*(.~(y)), particles.x, mask)
        particles.x .+= [particles.lb.*maskl[i] for i = 1:S]
        particles.x .+= [particles.ub.*masku[i] for i=1:S]

        # Update objectives and constraints
        particles.fx = map(x-> obj(x), particles.x)
        particles.fs = map(x-> is_feasible(x), particles.x)

        # Store particle's best position (if constraints are satisfied)
        update_position!(particles)

        # Busca local pelo algoritm SA conforme artigo
        verbose && println("Local Search.")
        # TODO: lembrar de no IHPSO fazer essa parte aqui ser um loop com as melhores posições
        indextemp = argmin(particles.fp)
        gtemp = copy(particles.p[indextemp])

        p_feasible = [is_feasible(particles.p[i]) for i = 1:S]
        if .|(p_feasible...)
            indextemp = argmin(particles.fp[p_feasible])
            gtemp = copy(particles.p[p_feasible][indextemp])
            indextemp = findfirst(x -> x == gtemp, particles.p)    
        end

        patemp = 1.
        
        itSA = 1.
        # TODO: remover esse print aqui depois de apresentar os resultados
        verbose && println("Local Search best at start $(gtemp) $(obj(gtemp)).")
        while itSA <= maxiterSA
            xtemp = gtemp + η * (particles.ub - particles.lb) .* randn(D)

            # Correct for bound violations
            maskltemp = xtemp .< particles.lb
            maskutemp = xtemp .> particles.ub 
            mask = maskltemp .| maskutemp
            xtemp = xtemp.*(.~(mask))
            xtemp .+= particles.lb.*maskltemp
            xtemp .+= particles.ub.*maskutemp

            fxtemp = obj(xtemp)
            if (is_feasible(xtemp) && is_feasible(gtemp))
                patemp = min(1., (exp((obj(gtemp) - fxtemp)/temperature)) ) 
            elseif (!is_feasible(xtemp) && is_feasible(gtemp))
                patemp = 0.0
            elseif (is_feasible(xtemp) && !is_feasible(gtemp))     
                patemp = 1.
            end

            if patemp > rand()
                #verbose && println("$(itSA) $(patemp) e $(obj(gtemp)) e $(obj(xtemp)) $(temperature)")
                gtemp = xtemp
            end

            itSA += 1
        end    
        particles.p[indextemp] = gtemp
        particles.fp[indextemp] = obj(gtemp)


        verbose && println("Local Search best at end $(gtemp) $(obj(gtemp)).")

        # Update swarm's best position
        findBestGlobal!(particles)
    

        # Compare swarm's best position with global best position
        i_min = findmin(particles.fp) 
        position_temp = copy(particles.p[i_min[2]])

        if .|(p_feasible...)
            i_min = findmin(particles.fp[p_feasible]) 
            position_temp = particles.p[p_feasible][i_min[2]]
        end


        if i_min[1] < particles.best_value
            it_best = it
            stepsize = √(sum((particles.best_position .- position_temp).^2))
            if abs.(particles.best_value .- i_min[1]) <= minfunc
                verbose && println("\n Stopping search: Swarm best objective change 
                less than $(minfunc)")
                break
            end
            if stepsize <= minstep
                verbose && println("Stopping search: Swarm best position change 
                less than $(minstep)")
                break
            end
            # if feasible change value
            if is_feasible(position_temp)
                particles.best_position = copy(position_temp) # best position
                particles.best_value = i_min[1] # best position value
            elseif !is_feasible(position_temp) && !is_feasible(particles.best_position)
                particles.best_position = copy(position_temp) # best position
                particles.best_value = i_min[1] # best position value
            end    
        end
        push!(particles.history, particles.best_value)
        verbose && (print("\r Iteration : $(it)/$(maxiter) | best-iteration $(it_best) 
        | Best value: $(particles.best_value) "))
        it += 1
    end
    verbose && println("\n Stopping search: maximum iterations reached --> $(maxiter)")
    is_feasible(particles.best_position) || print("However, the optimization couldn't find a
     feasible design. Sorry")
    return nothing
end


function hpso(func::Function, lb::Vector, ub::Vector; constraints=nothing, args=(), kwargs=Dict(),
    swarmsize=100, omegastart=0.9, omegaend=0.4, phip=0.5, phig=0.5, eta = 0.001, lambda = 0.94,
     maxiter=100, maxiterSA = 30, minstep=1e-8, minfunc=1e-8, verbose=false)

@assert length(ub) == length(lb) "ub and lb must have same dimension"
@assert all(ub .> lb) "Each value of 'ub' must be greater than 'lb'"
particles = Particles(swarmsize, lb, ub) #create particles

hpso(func, particles, constraints, args, kwargs, omegastart, omegaend, phip,
 phig, eta, lambda, maxiter, maxiterSA, minstep, minfunc, verbose)
return particles
end

function hpso(particles:: Particles, func::Function; constraints=nothing, args=(), kwargs=Dict(),
    swarmsize=100, omegastart=0.9, omegaend=0.4, phip=0.5, phig=0.5, eta = 0.001, lambda = 0.94,
    maxiter=100, maxiterSA = 30, minstep=1e-8, minfunc=1e-8, verbose=false)

hpso(func, particles, constraints, args, kwargs, omegastart, omegaend, phip,
 phig, eta, lambda, maxiter, maxiterSA, minstep, minfunc, verbose)

end
