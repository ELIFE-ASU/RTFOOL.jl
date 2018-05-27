using RTFOOL
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

re(p,q) = sum(map((p,q) -> (p != zero(p)) ? p*log(p/q) : 0.0, p, q))
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
end

@testset "Resources" begin
    @test_throws MethodError Resource([], [])
    @test_throws ArgumentError Resource([0.5], [0.0])
    @test_throws ArgumentError Resource([0.5, 0.25], [0.0, 1.0])
    @test_throws DimensionMismatch Resource([0.5, 0.5], [0.0])

    @test length(Resource([1], [0])) == 1
    @test length(Resource([0.5, 0.25, 0.25], [0, 1, 2])) == 3

    @test size(Resource([1], [0])) == (1,)
    @test size(Resource([0.5, 0.25, 0.25], [0, 1, 2])) == (3,)

    let r = Resource([0.1, 0.9], [1, 1e-1])
        @test r.p ≈ [0.1, 0.9]
        @test r.H ≈ [1, 1e-1]
    end

    let r = Resource([0.1, 0.5, 0.4], [1, 1e-1, 1e-2])
        @test r.p ≈ [0.1, 0.5, 0.4]
        @test r.H ≈ [1, 1e-1, 1e-2]
    end

    let r = Resource(2.0, [1, 10, 100])
        @test r.p ≈ [e^-2, e^-20, e^-200] / (e^-2 + e^-20 + e^-200)
        @test r.H ≈ [1, 10, 100]
    end

    let r = Resource(2.0, [1, 1, 2])
        @test r.p ≈ [e^2, e^2, 1] / (1 + 2e^2)
        @test r.H ≈ [1, 1, 2]
    end

    let r = Resource(1, [1,2])
        @test r.p ≈ [1.0, 0.0]
        @test r.H ≈ [1, 2]
    end

    let r = Resource(3, [2, 3, 3, 4])
        @test r.p ≈ [0.0, 0.0, 1.0, 0.0]
        @test r.H ≈ [2, 3, 3, 4]
    end
end

@testset "Subspace" begin
    @test_throws ArgumentError subspace([0], 1, [0.1, 1], 1)
    @test_throws ArgumentError subspace([1,2], -1, [0.1, 1], 1)
    @test_throws ArgumentError subspace([1,2], 0, [0.1, 1], 1)

    @test_throws ArgumentError subspace([1,2], 1, [1], 1)
    @test_throws ArgumentError subspace([1,2], 1, [0.1, 1], -1)
    @test_throws ArgumentError subspace([1,2], 1, [0.1, 1], 0)
    @test_throws ArgumentError subspace([1,2], 2, [0.1, 1], 1)

    #  let (basis, energy, deg) = subspace([1,2], 1)
    #      @test basis == [(1,)]
    #      @test energy ≈ [1]
    #      @test deg == [1]
    #  end
    #
    #  let (basis, energy, deg) = subspace([1,2], 3)
    #      @test basis == [(2,2,1), (1,1,1)]
    #      @test energy ≈ [5, 3]
    #      @test deg == [3, 1]
    #  end
    #
    #  let (basis, energy, deg) = subspace([1,2,3], 3)
    #      @test basis == [(3,3,3), (2,2,1), (1,1,1)]
    #      @test energy ≈ [9, 5, 3]
    #      @test deg == [3, 3, 1]
    #  end
    #
    #  let (basis, energy, deg) = subspace([1,2,3,4], 4)
    #      @test basis == [(4,4,4,4), (3,3,3,1), (2,2,2,2), (2,2,1,1), (1,1,1,1)]
    #      @test energy ≈ [16, 10, 8, 6, 4]
    #      @test deg == [12, 12, 3, 6, 1]
    #  end
    #
    #  let (basis, energy, deg) = subspace([1,2,3,4,5], 5)
    #      @test basis == [(5,5,5,5,5), (4,4,4,4,1), (3,3,3,2,2), (3,3,3,1,1), (2,2,2,2,1),
    #                      (2,2,1,1,1), (1,1,1,1,1)]
    #      @test energy ≈ [25, 17, 13, 11, 9, 7, 5]
    #      @test deg == [60, 60, 30, 30, 15, 10, 1]
    #  end
end
