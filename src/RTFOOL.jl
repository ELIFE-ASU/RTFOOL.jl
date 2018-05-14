module RTFOOL

export Resource

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

```jldoctest
julia> Resource([0.1, 0.9], [1, 2])
RTFOOL.Resource([0.1, 0.9], [1.0, 2.0])

julia> Resource(0.5, [1, 2])
RTFOOL.Resource([0.622459, 0.377541], [1.0, 2.0])
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

Base.size(r::Resource) = size(r.p)

Base.length(r::Resource) = length(r.p)

end
