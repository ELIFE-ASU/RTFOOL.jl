__precompile__(true)

module RTFOOL

using Combinatorics, DataStructures

export StateSpace, Context, degeneracies, timestep

include("space.jl")
include("boltzmann.jl")
include("context.jl")

end
