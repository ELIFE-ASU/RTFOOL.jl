module RTFOOL

"""
    boltzmann(β::Float64, E::Float64)

Compute the Boltzmann factor for inverse temperature β and energy E.
"""
boltzmann(β::Float64, E::Float64) = exp(-β*E)

end
