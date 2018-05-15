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
    if @views sum(ws[2:end]) < sum(ms[2:end]) - length(ms) + 1
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
    degeneracies(H::AbstractArray)

Construct an dictionary associating each state with the indicies of the states with whom it
is degenerate. Each state is only associated with those who's indicies are greater than its
own. The number of admissible swaps (excluding the identity) is also returned.

    degeneracies(r::Resource)

Construct the degeneracy dictionary given a resource.
"""
function degeneracies(H::AbstractArray)
    num_swaps = 0
    deg = Dict{Int, Vector{Int}}()
    for i in 1:length(H)
        for j in i+1:length(H)
            if H[i] ≈ H[j]
                num_swaps += 1
                if haskey(deg, i)
                    push!(deg[i], j)
                else
                    deg[i] = [j]
                end
            end
        end
    end
    deg, num_swaps
end
degeneracies(r::Resource) = degeneracies(r.H)

"""
    Context(β, Hm, Nm, Hw, Nw, system)

    Context(β, Hm, ms, Ws, ws)
"""
struct Context
    β::Float64

    monomer::Resource
    water::Resource

    bath::Resource
    Nm::Int
    Nw::Int

    system::Resource

    H::Vector{Float64}

    deg::Dict{Int, Vector{Int}}
    num_swaps::Int

    function Context(β, Hm, Nm, Hw, Nw, system)
        monomer = Resource(β, Hm)
        water = Resource(β, Hw)

        bath = tensor((monomer, Nm), (water, Nw))

        H = kron(system.H, ones(bath.H)) + kron(ones(system.H), bath.H)

        deg, num_swaps = degeneracies(H)

        new(β, monomer, water, bath, Nm, Nw, system, H, deg, num_swaps)
    end
end
Context(β, Hm, ms, Hw, ws) = let (system, Nm, Nw) = pure_system(Hm, ms, Hw, ws)
    Context(β, Hm, Nm, Hw, Nw, system)
end

"""
    randswap([rng=GLOBAL_RNG], ctx)

Select a valid swap at random from the given context.
"""
function randswap(rng::AbstractRNG, ctx::Context)
    n = rand(rng, 0:ctx.num_swaps)
    if n == 0
        (0, 0)
    else
        a, b = 0, 0
        for (k, v) in ctx.deg
            m = length(v)
            if n > m
                n -= m
            else
                a = k
                b = v[n]
                break
            end
        end
        (a, b)
    end
end
randswap(ctx::Context) = randswap(Base.GLOBAL_RNG, ctx)

end
