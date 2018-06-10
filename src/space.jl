struct StateSpace{N}
    basis :: Vector{NTuple{N, Int}}
    energy :: Vector{Float64}
    deg :: Vector{BigInt}
end

function StateSpace(Hm, Nm, Hw, Nw)
    const Lm, Lw = length(Hm), length(Hw)

    if Lm < 2
        throw(ArgumentError("monomer Hamiltonian is too short"))
    elseif Lw < 2
        throw(ArgumentError("water Hamiltonian is too short"))
    elseif Nm < 1
        throw(ArgumentError("space has too few monomers"))
    elseif Nw < Nm
        throw(ArgumentError("space has fewer water molecules than monomers"))
    end

    const Nmf, Nwf = factorial(BigInt(Nm)), factorial(BigInt(Nw))

    basis, energy, degeneracy = NTuple{Nm + Nw, Int}[], Float64[], BigInt[]

    state = Array{Int}(Nm + Nw)
    for partition in Combinatorics.IntegerPartitions(Nm)
        if maximum(partition) <= Lm
            i = 1
            for n in partition
                for _ in 1:n
                    state[i] = n
                    i += 1
                end
            end

            menergy = mapreduce(n -> n*Hm[n], +, partition)
            mdeg = Nmf
            for (n, a) in counter(partition)
                mdeg /= factorial(BigInt(a))
                if n != 1
                    mdeg /= 2^a
                end
            end

            t = length(partition)
            state[i:end-t] = Lw
            state[end-t+1:end] = 1
            push!(basis, tuple(state...))

            wenergy = sum(@view Hw[state[i:end]])
            push!(energy, menergy + wenergy)

            wdeg = Nwf
            for (n, a) in counter(@view state[i:end])
                wdeg ÷= factorial(BigInt(a))
            end
            push!(degeneracy, mdeg * wdeg)

            while state[i] > 2
                for j in (length(state)-t):-1:i
                    if state[j] > 2
                        state[j] -= 1
                        state[j+1:length(state)-t] = state[j]
                        break
                    end
                end
                push!(basis, tuple(state...))

                wenergy = sum(@view Hw[state[i:end]])
                push!(energy, menergy + wenergy)

                wdeg = Nwf
                for (n, a) in counter(@view state[i:end])
                    wdeg ÷= factorial(BigInt(a))
                end
                push!(degeneracy, mdeg * wdeg)
            end
        end
    end

    StateSpace{Nm + Nw}(basis, energy, degeneracy)
end

Base.length(s::StateSpace) = length(s.basis)

function degeneracies(s::StateSpace, t::StateSpace)
    deg_states = NTuple{2, NTuple{2,Int}}[]
    N = BigInt(0)
    A = BigInt(0)
    counts = BigInt[]
    for i in 1:length(s)
        for j in 1:length(t)
            e1 = s.energy[i] + t.energy[j]
            d1 = s.deg[i] * t.deg[j]
            A += 0.5 * d1 * (d1 - 1)
            for l in (j+1):length(t)
                e2 = s.energy[i] + t.energy[l]
                d2 = s.deg[i] * t.deg[l]
                if e1 == e2
                    push!(deg_states, ((i,j), (i,l)))
                    push!(counts, d1*d2)
                    N += counts[end]
                end
            end
            for k in (i+1):length(s)
                for l in 1:length(t)
                    e2 = s.energy[k] + t.energy[l]
                    d2 = s.deg[k] * t.deg[l]
                    if e1 == e2
                        push!(deg_states, ((i,j), (k,l)))
                        push!(counts, d1*d2)
                        N += counts[end]
                    end
                end
            end
        end
    end
    N += A
    deg_states, counts/N, A/N
end
