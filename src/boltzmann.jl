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
