__precompile__(true)

module RTFOOL

using Combinatorics, DataStructures

export StateSpace, Context, degeneracies, timestep
export Monotones

include("space.jl")
include("boltzmann.jl")
include("context.jl")
include("monotones.jl")

end
