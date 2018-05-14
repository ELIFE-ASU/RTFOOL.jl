module RTFOOL

export Resource, tensor, pure_system
export Context

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

"""
    pure_system(Hm, ms, Hw, ws)

Construct a pure system from a collection of monomers and water molecules with Hamiltonians
`Hm` and `Hw`, respectively. The number of free monomers, dimers, trimers, etc... are
provided in `ms`, and `ws` provides the number of dissociated water molecules, and higher
energy states.

The result is a tuple of the from `(r, Nm, Nw)` where `r` is the constructed resource, `Nm`
is the total number of monomers (bound or otherwise), and `Nw` is the total number of water
molecules.

```jldoctest
julia> r, Nm, Nw = pure_system([0.1, 1.0], [0,1], [0.1, 1.0], [0,2]);

julia> Nm, Nw
(2, 2)

julia> r.p
16-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 1.0
```
"""
function pure_system(Hm, ms, Hw, ws)
    if length(Hm) != length(ms)
        throw(DimensionMismatch("length(Hm) != length(ms)"))
    elseif length(Hw) != length(ws)
        throw(DimensionMismatch("length(Hw) != length(ws)"))
    elseif ws[1] < ms[1]
        throw(ArgumentError("too few dissocated water molecules given number of free monomers"))
    end 

    ms = map(*, 1:length(ms), ms)
    if @views sum(ws[2:end]) < sum(ms[2:end])
        warn("""It is suggested that the number of associated water molecules be at least
equal to the total number of bonds in the system""")
    end

    monomers = map(z->(Resource(z[1], Hm), z[2]), enumerate(ms))
    Nm = sum(ms)

    waters = map(z -> (Resource(z[1], Hw), z[2]), enumerate(ws))
    Nw = sum(ws)

    tensor(monomers..., waters...), Nm, Nw
end

Base.size(r::Resource) = size(r.p)

Base.length(r::Resource) = length(r.p)

Base.isapprox(x::Resource, y::Resource) = x.p ≈ y.p && x.H ≈ y.H

"""
    tensor(x::Tuple{Resource,Int}, xs::Tuple{Resource,Int}...)

Construct the tensor product of resources.

```jldoctest
julia> tensor((Resource(0.25, [1,2]), 2))
RTFOOL.Resource([0.316042, 0.246134, 0.246134, 0.191689], [2.0, 3.0, 3.0, 4.0])

julia> tensor((Resource([0,1], [1,2]), 2), (Resource([0,1,0], [1,2,3]), 1))
RTFOOL.Resource([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0], [3.0, 4.0, 5.0, 4.0, 5.0, 6.0, 4.0, 5.0, 6.0, 5.0, 6.0, 7.0])
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

"""
    Context(β, Hm, Nm, Hw, Nw)

Construct a Context with a given inverse temperature `β`, and monomer and water Hamiltonians
`Hm` and `Hw`, respectively. A bath is constructed with `Nm` monomers and `Nw` water
molecules in Gibb's states.

```jldoctest
julia> ctx = Context(0.25, [1,2], 2, [2,3], 1);

julia> ctx.β
0.25

julia> ctx.monomer
RTFOOL.Resource([0.562177, 0.437823], [1.0, 2.0])

julia> ctx.water
RTFOOL.Resource([0.562177, 0.437823], [2.0, 3.0])

julia> ctx.bath
RTFOOL.Resource([0.177672, 0.138371, 0.138371, 0.107763, 0.138371, 0.107763, 0.107763, 0.0839261], [4.0, 5.0, 5.0, 6.0, 5.0, 6.0, 6.0, 7.0])
```
"""
struct Context
    β::Float64

    monomer::Resource
    water::Resource

    bath::Resource
    Nm::Int
    Nw::Int

    function Context(β, Hm, Nm, Hw, Nw)
        monomer = Resource(β, Hm)
        water = Resource(β, Hw)
        bath = tensor((monomer, Nm), (water, Nw))
        new(β, monomer, water, bath, Nm, Nw)
    end
end

end
