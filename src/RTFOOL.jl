__precompile__(true)

module RTFOOL

using Combinatorics, DataStructures

export StateSpace, Context

include("space.jl")
include("boltzmann.jl")
include("context.jl")

end
