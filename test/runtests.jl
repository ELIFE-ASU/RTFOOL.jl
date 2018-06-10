using RTFOOL
using RTFOOL.Monotones
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

include("space.jl")
include("boltzmann.jl")
include("context.jl")
include("monotones.jl")
