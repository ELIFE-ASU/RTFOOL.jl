__precompile__(true)

module RTFOOL

using Combinatorics, DataStructures

export Resource, subspace

"""
    boltzmann(β, E::Number)

Compute the Boltzmann factor for inverse temperature `β` and energy `E`.

    boltzmann(β, E)

Construct the Boltzmann distributions with inverse temperature `β` over the energy states
`E`.
"""
boltzmann(β, E::Number) = exp(-β*E)

function boltzmann(β, E)
    factors = boltzmann.(β, E)
    partition = sum(factors)
    factors ./ partition
end

"""
    Resource(p, H)

Construct a resource with state (probability array) `p` and Hamiltonian (energy array) `H`.

    Resource(β, H)

Construct a resource with a Boltzmann distribution (at inverse temperature `β`) over the
Hamiltonian `H`.

    Resource(i, H)

Construct a pure-state resource over the Hamiltonian `H`

```jldoctest
julia> Resource([0.1, 0.9], [1, 2])
RTFOOL.Resource([0.1, 0.9], [1.0, 2.0])

julia> Resource(0.5, [1, 2])
RTFOOL.Resource([0.622459, 0.377541], [1.0, 2.0])

julia> Resource(1, [1, 2])
RTFOOL.Resource([1.0, 0.0], [1.0, 2.0])

julia> Resource(3, [2, 3, 3, 4])
RTFOOL.Resource([0.0, 0.0, 1.0, 0.0], [2.0, 3.0, 3.0, 4.0])
```
"""
struct Resource
    p::Vector{Float64}
    H::Vector{Float64}

    Resource(p, H) = if !(sum(p) ≈ 1.0)
        throw(ArgumentError("probability array is not normalized"))
    elseif size(p) != size(H)
        throw(DimensionMismatch("size(p) != size(H)"))
    else
        new(p, H)
    end
end
Resource(β::Float64, H) = Resource(boltzmann(β, H), H)
Resource(i::Int, H) = let p = zeros(size(H))
    p[i] = 1.0
    Resource(p, H)
end

Base.size(r::Resource) = size(r.p)

Base.length(r::Resource) = length(r.p)

Base.isapprox(x::Resource, y::Resource) = x.p ≈ y.p && x.H ≈ y.H

function subspace(Hm, Nm, Hw, Nw)
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

    basis, energy, degeneracy
end

end
