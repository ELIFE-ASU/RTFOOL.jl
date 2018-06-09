__precompile__(true)

module RTFOOL

using Combinatorics, DataStructures

export Resource, StateSpace

struct StateSpace{N}
    basis :: Vector{NTuple{N, Int}}
    energy :: Vector{Float64}
    deg :: Vector{BigInt}
end

function StateSpace(Hm, Nm, Hw, Nw)
    const Lm, Lw = length(Hm), length(Hw)

    if Lm < 2
        throw(ArgumentError("monomer Hamiltonian is too short"))
    elseif Lw < 2
        throw(ArgumentError("water Hamiltonian is too short"))
    elseif Nm < 1
        throw(ArgumentError("space has too few monomers"))
    elseif Nw < Nm
        throw(ArgumentError("space has fewer water molecules than monomers"))
    end

    const Nmf, Nwf = factorial(BigInt(Nm)), factorial(BigInt(Nw))

    basis, energy, degeneracy = NTuple{Nm + Nw, Int}[], Float64[], BigInt[]

    state = Array{Int}(Nm + Nw)
    for partition in Combinatorics.IntegerPartitions(Nm)
        if maximum(partition) <= Lm
            i = 1
            for n in partition
                for _ in 1:n
                    state[i] = n
                    i += 1
                end
            end

            menergy = mapreduce(n -> n*Hm[n], +, partition)
            mdeg = Nmf
            for (n, a) in counter(partition)
                mdeg /= factorial(BigInt(a))
                if n != 1
                    mdeg /= 2^a
                end
            end

            t = length(partition)
            state[i:end-t] = Lw
            state[end-t+1:end] = 1
            push!(basis, tuple(state...))

            wenergy = sum(@view Hw[state[i:end]])
            push!(energy, menergy + wenergy)

            wdeg = Nwf
            for (n, a) in counter(@view state[i:end])
                wdeg ÷= factorial(BigInt(a))
            end
            push!(degeneracy, mdeg * wdeg)

            while state[i] > 2
                for j in (length(state)-t):-1:i
                    if state[j] > 2
                        state[j] -= 1
                        state[j+1:length(state)-t] = state[j]
                        break
                    end
                end
                push!(basis, tuple(state...))

                wenergy = sum(@view Hw[state[i:end]])
                push!(energy, menergy + wenergy)

                wdeg = Nwf
                for (n, a) in counter(@view state[i:end])
                    wdeg ÷= factorial(BigInt(a))
                end
                push!(degeneracy, mdeg * wdeg)
            end
        end
    end

    StateSpace{Nm + Nw}(basis, energy, degeneracy)
end

boltzmann(β, E::Number) = exp(-β*E)

function boltzmann(β, E)
    factors = boltzmann.(β, E)
    partition = sum(factors)
    factors ./ partition
end

function boltzmann(β, E, deg)
    factors = map(*, boltzmann.(β, E), deg)
    partition = sum(factors)
    convert(Array{Float64}, factors ./ partition)
end

boltzmann(β, space::StateSpace) = boltzmann(β, space.energy, space.deg)

end
