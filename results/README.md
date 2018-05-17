---
title: "Laboratory Notebook"
date: "May 2018"
author: "Nicole Y. Halpern, Nicolas Mathis, Douglas G. Moore, Sara I. Walker"
---

# Setup the project - 2018-05-13

DGM setup the repository, continuous integration, API documentation and started
creating some basic functions, e.g. `RTFOOL.boltzmann` which computes
Boltzmann's factor for a given inverse temperature and energy.

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
