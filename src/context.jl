struct Context{Ns, Nb}
    β :: Float64

    system_space :: StateSpace{Ns}
    system_state :: Vector{Float64}

    bath_space :: StateSpace{Nb}
    bath_state :: Vector{Float64}

    function Context{Ns, Nb}(β::Float64, system_space::StateSpace{Ns}, system_state::AbstractVector,
                             bath_space::StateSpace{Nb}) where{Ns, Nb}
        new(β, system_space, system_state, bath_space, boltzmann(β, bath_space))
    end
end

function Context(β::Float64, system_space::StateSpace{N}, system_state::AbstractVector) where {N}
    Context{N,N}(β, system_space, system_state, system_space)
end
