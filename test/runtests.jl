using RTFOOL
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

re(p, q) = sum(map((p,q) -> (p != zero(p)) ? p*log(p/q) : zero(p), p, q))
entropy(p) = -dot(p, map(p -> (p != zero(p)) ? log(p) : 0.0, p))

@testset "Boltzmann" begin
    @test RTFOOL.boltzmann(1.0, 1.0) == e^-1.00
    @test RTFOOL.boltzmann(1.0, 2.0) == e^-2.00
    @test RTFOOL.boltzmann(0.5, 1.5) == e^-0.75

    @test RTFOOL.boltzmann(1.0, [1.0]) == [1.0]
    @test RTFOOL.boltzmann(1.0, [2.0]) == [1.0]
    @test RTFOOL.boltzmann(1.0, [1.0, 1.0]) == [0.5, 0.5]
    @test RTFOOL.boltzmann(1.0, [1.0, 2.0]) == [e/(1+e), 1/(1+e)]
    @test RTFOOL.boltzmann(0.5, [0.5, 1.5]) ≈ [e/(e + √e), √e/(e + √e)]
    @test RTFOOL.boltzmann(1/3, [1, 2, 3]) ==
        [e^-(1/3), e^-(2/3), e^-1] / (e^-(1/3) + e^-(2/3) + e^-1)

    @test RTFOOL.boltzmann(1.0, [1.0], [2]) == [1.0]
    @test RTFOOL.boltzmann(1.0, [2.0], [3]) == [1.0]
    @test RTFOOL.boltzmann(2.0, [1.0, 1.0], [1, 2]) == [1.0/3, 2.0/3]
    @test RTFOOL.boltzmann(1.0, [1.0, 2.0], [2, 1]) == [2e/(1+2e), 1/(1+2e)]
    @test RTFOOL.boltzmann(0.5, [0.5, 1.5], [2, 3]) ≈ [2e/(2e + 3*√e), 3*√e/(2e + 3*√e)]
    @test RTFOOL.boltzmann(1/3, [1, 2, 3], [1, 2, 3]) ==
        [e^-(1/3), 2e^-(2/3), 3e^-1] / (e^-(1/3) + 2e^-(2/3) + 3e^-1)

    let space = StateSpace{1}([(1,)], [1.0], BigInt[2])
        @test RTFOOL.boltzmann(1.0, space) == [1.0]
    end
    let space = StateSpace{1}([(1,)], [2.0], BigInt[3])
        @test RTFOOL.boltzmann(1.0, [2.0], [3]) == [1.0]
    end
    let space = StateSpace{2}([(1,1), (2,2)], [1.0, 1.0], BigInt[1,2])
        @test RTFOOL.boltzmann(2.0, space) == [1.0/3, 2.0/3]
    end
    let space = StateSpace{2}([(1,1),(2,2)], [1.0, 2.0], BigInt[2,1])
        @test RTFOOL.boltzmann(1.0, space) == [2e/(1+2e), 1/(1+2e)]
    end
    let space = StateSpace{2}([(1,1),(2,2)], [0.5, 1.5], BigInt[2,3])
        @test RTFOOL.boltzmann(0.5, space) ≈ [2e/(2e + 3*√e), 3*√e/(2e + 3*√e)]
    end
    let space = StateSpace{3}([(1,1,1), (1,1,2), (1,2,2)], [1,2,3], BigInt[1,2,3])
        @test RTFOOL.boltzmann(1/3, space) ≈
            [e^-(1/3), 2e^-(2/3), 3e^-1] / (e^-(1/3) + 2e^-(2/3) + 3e^-1)
    end
end

@testset "StateSpace" begin
    @test_throws ArgumentError StateSpace([0], 1, [0.1, 1], 1)
    @test_throws ArgumentError StateSpace([1,2], -1, [0.1, 1], 1)
    @test_throws ArgumentError StateSpace([1,2], 0, [0.1, 1], 1)

    @test_throws ArgumentError StateSpace([1,2], 1, [1], 1)
    @test_throws ArgumentError StateSpace([1,2], 1, [0.1, 1], -1)
    @test_throws ArgumentError StateSpace([1,2], 1, [0.1, 1], 0)
    @test_throws ArgumentError StateSpace([1,2], 2, [0.1, 1], 1)

    let space = StateSpace([1,2], 1, [0.1, 1.0], 1)
        @test space.basis == [(1,1)]
        @test space.energy ≈ [1.1]
        @test space.deg == [1]
    end

    let space = StateSpace([1,2], 1, [0.1, 1.0], 2)
        @test space.basis == [(1,2,1)]
        @test space.energy ≈ [2.1]
        @test space.deg == [2]
    end

    let space = StateSpace([1,2], 1, [0.1, 1.0, 10.0], 2)
        @test space.basis == [(1,3,1), (1,2,1)]
        @test space.energy ≈ [11.1, 2.1]
        @test space.deg == [2,2]
    end

    let space = StateSpace([1,2], 2, [0.1, 1.0], 2)
        @test space.basis == [(2,2,2,1), (1,1,1,1)]
        @test space.energy ≈ [5.1, 2.2]
        @test space.deg == [2,1]
    end

    let space = StateSpace([1,2], 2, [0.1, 1.0, 10.0], 2)
        @test space.basis == [(2,2,3,1), (2,2,2,1), (1,1,1,1)]
        @test space.energy ≈ [14.1, 5.1, 2.2]
        @test space.deg == [2,2,1]
    end

    let space = StateSpace([1,2], 2, [0.1, 1.0], 3)
        @test space.basis == [(2,2,2,2,1), (1,1,2,1,1)]
        @test space.energy ≈ [6.1, 3.2] 
        @test space.deg == [3,3]
    end

    let space = StateSpace([1,2], 2, [0.1, 1.0, 10.0], 3)
        @test space.basis == [(2,2,3,3,1), (2,2,3,2,1), (2,2,2,2,1), (1,1,3,1,1), (1,1,2,1,1)]
        @test space.energy ≈ [24.1, 15.1, 6.1, 12.2, 3.2]
        @test space.deg ==   [3,6,3,3,3]
    end

    let space = StateSpace([1,2], 3, [0.1, 1.0], 3)
        @test space.basis == [(2,2,1,2,1,1), (1,1,1,1,1,1)]
        @test space.energy ≈ [6.2, 3.3]
        @test space.deg == [9,1]
    end

    let space = StateSpace([1,2], 3, [0.1, 1.0, 10.0], 3)
        @test space.basis == [(2,2,1,3,1,1), (2,2,1,2,1,1), (1,1,1,1,1,1)]
        @test space.energy ≈ [15.2, 6.2, 3.3]
        @test space.deg == [9, 9, 1]
    end

    let space = StateSpace([1,2,3], 3, [0.1, 1.0], 3)
        @test space.basis == [(3,3,3,2,2,1), (2,2,1,2,1,1), (1,1,1,1,1,1)]
        @test space.energy ≈ [11.1, 6.2, 3.3]
        @test space.deg == [9, 9, 1]
    end

    let space = StateSpace([1,2,3], 3, [0.1, 1.0, 10.0], 3)
        @test space.basis == [(3,3,3,3,3,1), (3,3,3,3,2,1), (3,3,3,2,2,1), (2,2,1,3,1,1),
                              (2,2,1,2,1,1), (1,1,1,1,1,1)]
        @test space.energy ≈ [29.1, 20.1, 11.1, 15.2, 6.2, 3.3]
        @test space.deg == [9, 18, 9, 9, 9, 1]
    end

    let space = StateSpace([1,2,3,4], 4, [0.1, 1.0], 4)
        @test space.basis == [(4,4,4,4,2,2,2,1), (3,3,3,1,2,2,1,1),
                              (2,2,2,2,2,2,1,1), (2,2,1,1,2,1,1,1),
                              (1,1,1,1,1,1,1,1)]
        @test space.energy ≈ [19.1, 12.2, 10.2, 7.3, 4.4]
        @test space.deg == [48, 72, 18, 24, 1]
    end

    let space = StateSpace([1,2,3,4], 4, [0.1, 1.0, 10.0], 4)
        @test space.basis == [(4,4,4,4,3,3,3,1), (4,4,4,4,3,3,2,1), (4,4,4,4,3,2,2,1),
                              (4,4,4,4,2,2,2,1), (3,3,3,1,3,3,1,1), (3,3,3,1,3,2,1,1),
                              (3,3,3,1,2,2,1,1), (2,2,2,2,3,3,1,1), (2,2,2,2,3,2,1,1),
                              (2,2,2,2,2,2,1,1), (2,2,1,1,3,1,1,1), (2,2,1,1,2,1,1,1),
                              (1,1,1,1,1,1,1,1)]
        @test space.energy ≈ [46.1, 37.1, 28.1, 19.1, 30.2, 21.2, 12.2, 28.2, 19.2, 10.2, 16.3,
                              7.3, 4.4]
        @test space.deg == [48, 144, 144, 48, 72, 144, 72, 18, 36, 18, 24, 24, 1]
    end

    let space = StateSpace([1,2,3,4,5], 5, [0.1, 1.0], 5)
        @test space.basis == [(5,5,5,5,5,2,2,2,2,1), (4,4,4,4,1,2,2,2,1,1), (3,3,3,2,2,2,2,2,1,1),
                              (3,3,3,1,1,2,2,1,1,1), (2,2,2,2,1,2,2,1,1,1), (2,2,1,1,1,2,1,1,1,1),
                              (1,1,1,1,1,1,1,1,1,1)]
         @test space.energy ≈ [29.1, 20.2, 16.2, 13.3, 11.3, 8.4, 5.5]
        @test space.deg == [300, 600, 300, 300, 150, 50, 1]
    end
end

@testset "Context" begin
    let space = StateSpace([1,2], 2, [0.1, 1.0], 2)
        ctx = Context(0.5, space, [0,1])
        @test ctx.β == 0.5
        @test ctx.system_space === space
        @test ctx.system_state == [0,1]
        @test ctx.bath_space == space
        let a = e^-2.55, b = e^-1.1
            @test ctx.bath_state ≈ [2a/(2a + b), b/(2a + b)]
        end
    end
end

@testset "Degeneracies" begin
    let sys = StateSpace([1,2], 2, [0.1, 1.0], 2)
        deg, probs, pid = degeneracies(sys, sys)
        @test deg == [((1,2), (2,1))]
        @test probs ≈ [1/3]
        @test sum(probs) + pid ≈ 1.0
    end
    let sys = StateSpace([1,2], 4, [0.1, 1.0], 4)
        deg, probs, pid = degeneracies(sys, sys)
        @test deg == [((1,2), (2,1)), ((1,3), (3,1)), ((2,3), (3,2))]
        @test probs ≈ [186624, 324, 576] / 592500
        @test sum(probs) + pid ≈ 1.0
    end
    let sys = StateSpace(1:5, 5, [0.1, 1.0], 5)
        deg, probs, pid = degeneracies(sys, sys)
        @test deg == [((1, 2), (2, 1)), ((1, 3), (3, 1)), ((1, 4), (4, 1)), ((1, 5),
        (5, 1)), ((1, 6), (6, 1)), ((1, 7), (7, 1)), ((2, 3), (3, 2)), ((2, 4), (4, 2)),
        ((2, 5), (5, 2)), ((2, 6), (6, 2)), ((2, 7), (7, 2)), ((3, 4), (4, 3)), ((3, 5),
        (5, 3)), ((3, 6), (4, 5)), ((3, 6), (5, 4)), ((3, 6), (6, 3)), ((3, 7), (7, 3)),
        ((4, 5), (5, 4)), ((4, 5), (6, 3)), ((4, 6), (6, 4)), ((4, 7), (7, 4)), ((5, 4),
        (6, 3)), ((5, 6), (6, 5)), ((5, 7), (6, 6)), ((5, 7), (7, 5)), ((6, 6), (7, 5)),
        ((6, 7), (7, 6))]
        @test isapprox(probs, [0.0913913, 0.0228478, 0.0228478, 0.00571196, 0.000634662,
        2.53865e-7, 0.0913913, 0.0913913, 0.0228478, 0.00253865, 1.01546e-6, 0.0228478,
        0.00571196, 0.00190399, 0.00190399, 0.000634662, 2.53865e-7, 0.00571196,
        0.00190399, 0.000634662, 2.53865e-7, 0.00190399, 0.000158666, 1.05777e-6,
        6.34662e-8, 1.05777e-6, 7.0518e-9]; atol=1e-7)
        @test sum(probs) + pid ≈ 1.0
    end
end

let rng = MersenneTwister(2018)
    @testset "Timestep" begin
        let ctx = Context(0.5, StateSpace([1,2], 4, [0.1, 1.0], 4), [1, 0, 0])
            r = re(ctx.system_state, ctx.bath_state)
            for _ in 1:10000
                timestep(rng, ctx)
                rnext = re(ctx.system_state, ctx.bath_state)
                @test isapprox(rnext, r; atol=1e-9) || rnext < r
                r = rnext
            end
            @test ctx.system_state ≈ ctx.bath_state
        end
        let ctx = Context(0.5, StateSpace(1:4, 4, [0.1, 1.0], 4), [1,0,0,0,0])
            r = re(ctx.system_state, ctx.bath_state)
            for _ in 1:1000000
                timestep(rng, ctx)
                rnext = re(ctx.system_state, ctx.bath_state)
                @test isapprox(rnext, r; atol=1e-9) || rnext < r
                r = rnext
            end
            @test ctx.system_state ≈ ctx.bath_state
        end
    end
end
