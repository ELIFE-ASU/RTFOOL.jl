module RTFOOL

export Resource, tensor

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

"""
    tensor(x::Tuple{Resource,Int}, xs::Tuple{Resource,Int}...)

Construct the tensor product of resources.

```jldoctest
julia> tensor((Resource(0.25, [1,2]), 2))
RTFOOL.Resource([0.316042, 0.246134, 0.246134, 0.191689], [2.0, 3.0, 3.0, 4.0])

julia> tensor((Resource([0,1], [1,2]), 2), (Resource([0,1,0], [1,2,3]), 1))
RTFOOL.Resource([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
                [3.0, 4.0, 5.0, 4.0, 5.0, 6.0, 4.0, 5.0, 6.0, 5.0, 6.0, 7.0])
```
"""
function tensor(x::Tuple{Resource,Int}, xs::Tuple{Resource,Int}...)
    size = mapreduce(t -> length(t[1])^t[2], *, length(x[1])^x[2], xs)

    p = ones(size)
    H = zeros(size)

    stretch = size
    let (r, n) = x
        const m = length(r)
        for j in 1:n
            stretch ÷= m
            for i in 0:size-1
                k = (i ÷ stretch) % m
                p[i+1] *= r.p[k+1]
                H[i+1] += r.H[k+1]
            end
        end
    end

    for (r, n) in xs
        const m = length(r)
        for j in 1:n
            stretch ÷= m
            for i in 0:size-1
                k = (i ÷ stretch) % m
                p[i+1] *= r.p[k+1]
                H[i+1] += r.H[k+1]
            end
        end
    end

    Resource(p, H)
end

end
