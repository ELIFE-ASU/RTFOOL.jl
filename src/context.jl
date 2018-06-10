struct Context{Ns, Nb}
    β :: Float64

    system_space :: StateSpace{Ns}
    system_state :: Vector{Float64}

    bath_space :: StateSpace{Nb}
    bath_state :: Vector{Float64}

    degeneracies :: Vector{NTuple{2, NTuple{2, Int}}}
    P :: Vector{BigFloat}
    Pid :: BigFloat

    function Context{Ns, Nb}(β::Float64, system_space::StateSpace{Ns}, system_state::AbstractVector,
                             bath_space::StateSpace{Nb}) where{Ns, Nb}
        bath_state = boltzmann(β, bath_space)
        deg, P, Pid = degeneracies(system_space, bath_space)
        new(β, system_space, system_state, bath_space, bath_state, deg, P, Pid)
    end
end

function Context(β::Float64, system_space::StateSpace{N}, system_state::AbstractVector) where {N}
    Context{N,N}(β, system_space, system_state, system_space)
end

function timestep(rng::AbstractRNG, ctx::Context)
    p = rand(rng)
    if p > ctx.Pid
        p -= ctx.Pid
        i = 1
        while i < length(ctx.P) && p > ctx.P[i]
            p -= ctx.P[i]
            i += 1
        end
        if i <= length(ctx.P)
            ((a,b), (u,v)) = ctx.degeneracies[i]
            Pa, Pb = ctx.system_state[a], ctx.bath_state[b]
            Pu, Pv = ctx.system_state[u], ctx.bath_state[v]
            dP = Pu*Pv - Pa*Pb
            ctx.system_state[a] += dP
            ctx.system_state[u] -= dP
        end
    end
    ctx.system_state
end

timestep(ctx::Context) = timestep(Base.GLOBAL_RNG, ctx)
