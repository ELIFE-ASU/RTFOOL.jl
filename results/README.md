---
title: "Laboratory Notebook"
date: "May 2018"
author: "Nicole Y. Halpern, Nicolas Mathis, Douglas G. Moore, Sara I. Walker"
---

# Full Rebuild - 2018-06-10

DGM went about rebuilding the framework to better handle large(ish)
systems. RTFOOL is now built out of two primary components: `StateSpace`
and `Context`.

## The StateSpace

For the purposes of this work, a state space is built of three components:
```julia
type StateSpace
	basis
	energy
	deg
end
```

The `basis` field holds the pure configurations of the system as a
tuple of monomial and water configurations, i.e. a tuple filled with the length
of the polymer each monomer is in, and the state of each water molecule. The
`energy` field holds the energy of each configuration, and the `deg` field
holds the degeneracy, i.e. number of microstates, of each basis configuration.

The state space may correspond to a micro-scale description (if the degeneracy
of each of the basis configurations is unity) or a macro-scale otherwise.
The primary constructor for this type constructs a macro-scale description,
but we could develop micro-level constructors as well.

## The Context

The `Context` type represents a thermalization context, i.e. a system and a
bath together with all of the metadata necessary to perform thermalization
time steps. At construction, the bath's Gibbs state is built, the admissible
swaps and their respective probabilities are constructed. The `timestep`
method can them perform individual thermalization timesteps, computing the
next state of the system (stored as the `system_state` field).

## Monotones

We also implemented the `RTFOOL.Monotones` module which provides a suite of
possible monotones. At the moment, only `relative_entropy` is implemented.

# Implement Context - 2018-05-15

DGM created the `Resource` and `Context` types together with a small collection
of utility functions and constructors to implement the basic simulation of this
resource theory. The `Resource` type provides a basic type for encapsulating the
state and Hamiltonian for each resource, while `Context` combines a bath and a
system which can be stocastically transformed via `transform!` to implement
simulation.

A minimal example might look like
```julia
const β  = 1e-2         # inverse temperature
const Hm = [0.50, 0.75] # monomer Hamiltonian
const Hw = [0.00, 0.10] # water Hamiltonian
const Nm = [0,1]        # start with a pure-state dimer
const Nw = [0,1]        # start with a pure-state liquid water

ctx = Context(1e-2, Hm, Nm, Hw, Nw)
for _ in 1:100
    transform!(ctx)
end
println(ctx.system.p)
```

There are several changes that need to be made before we can call this initial
development phase complete:
1. *Remove unphysical states from the tensor product spaces.* MOST of the states
   we are dealing with are unphysical and those states end up with non-zero
   probabilities in the bath. We should only consider the span of the physical
   states instead of the full tensor space.
2. *Do not construct the composite state of the system-bath.* The permutations
   can be applied directly to the system's state if we are careful.

# Setup the project - 2018-05-13

DGM setup the repository, continuous integration, API documentation and started
creating some basic functions, e.g. `RTFOOL.boltzmann` which computes
Boltzmann's factor for a given inverse temperature and energy.
