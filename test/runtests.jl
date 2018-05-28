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

    let (basis, energy, deg) = subspace([1,2], 1, [0.1, 1.0], 1)
        @test basis == [(1,1)]
        @test energy ≈ [1.1]
        #  @test deg == [1]
    end

    let (basis, energy, deg) = subspace([1,2], 1, [0.1, 1.0], 2)
        @test basis == [(1,2,1), (1,1,1)]
        @test energy ≈ [2.1, 1.2]
    end

    let (basis, energy, deg) = subspace([1,2], 1, [0.1, 1.0, 10.0], 2)
        @test basis == [(1,3,1), (1,2,1), (1,1,1)]
        @test energy ≈ [11.1, 2.1, 1.2]
    end

    let (basis, energy, deg) = subspace([1,2], 2, [0.1, 1.0], 2)
        @test basis == [(2,2,2,1), (2,2,1,1), (1,1,1,1)]
        @test energy ≈ [5.1, 4.2, 2.2]
    end

    let (basis, energy, deg) = subspace([1,2], 2, [0.1, 1.0, 10.0], 2)
        @test basis == [(2,2,3,1), (2,2,2,1), (2,2,1,1), (1,1,1,1)]
        @test energy ≈ [14.1, 5.1, 4.2, 2.2]
    end

    let (basis, energy, deg) = subspace([1,2], 2, [0.1, 1.0], 3)
        @test basis == [(2,2,2,2,1), (2,2,2,1,1), (2,2,1,1,1), (1,1,2,1,1), (1,1,1,1,1)]
        @test energy ≈ [6.1, 5.2, 4.3, 3.2, 2.3] 
    end

    let (basis, energy, deg) = subspace([1,2], 2, [0.1, 1.0, 10.0], 3)
        @test basis == [(2,2,3,3,1), (2,2,3,2,1), (2,2,3,1,1), (2,2,2,2,1), (2,2,2,1,1),
                        (2,2,1,1,1), (1,1,3,1,1), (1,1,2,1,1), (1,1,1,1,1)]
        @test energy ≈ [24.1, 15.1, 14.2, 6.1, 5.2, 4.3, 12.2, 3.2, 2.3]
    end

    let (basis, energy, deg) = subspace([1,2], 3, [0.1, 1.0], 3)
        @test basis == [(2,2,1,2,1,1), (2,2,1,1,1,1), (1,1,1,1,1,1)]
        @test energy ≈ [6.2, 5.3, 3.3]
        #  @test deg == [3, 1]
    end

    let (basis, energy, deg) = subspace([1,2], 3, [0.1, 1.0, 10.0], 3)
        @test basis == [(2,2,1,3,1,1), (2,2,1,2,1,1), (2,2,1,1,1,1), (1,1,1,1,1,1)]
        @test energy ≈ [15.2, 6.2, 5.3, 3.3]
        #  @test deg == [3, 1]
    end

    let (basis, energy, deg) = subspace([1,2,3], 3, [0.1, 1.0], 3)
        @test basis == [(3,3,3,2,2,1), (3,3,3,2,1,1), (3,3,3,1,1,1), (2,2,1,2,1,1),
                        (2,2,1,1,1,1), (1,1,1,1,1,1)]
        @test energy ≈ [11.1, 10.2, 9.3, 6.2, 5.3, 3.3]
    #      @test deg == [3, 3, 1]
    end

    let (basis, energy, deg) = subspace([1,2,3], 3, [0.1, 1.0, 10.0], 3)
        @test basis == [(3,3,3,3,3,1), (3,3,3,3,2,1), (3,3,3,3,1,1), (3,3,3,2,2,1),
                        (3,3,3,2,1,1), (3,3,3,1,1,1), (2,2,1,3,1,1), (2,2,1,2,1,1),
                        (2,2,1,1,1,1), (1,1,1,1,1,1)]
        @test energy ≈ [29.1, 20.1, 19.2, 11.1, 10.2, 9.3, 15.2, 6.2, 5.3, 3.3]
    #      @test deg == [3, 3, 1]
    end

    let (basis, energy, deg) = subspace([1,2,3,4], 4, [0.1, 1.0], 4)
        @test basis == [(4,4,4,4,2,2,2,1), (4,4,4,4,2,2,1,1), (4,4,4,4,2,1,1,1),
                        (4,4,4,4,1,1,1,1), (3,3,3,1,2,2,1,1), (3,3,3,1,2,1,1,1),
                        (3,3,3,1,1,1,1,1), (2,2,2,2,2,2,1,1), (2,2,2,2,2,1,1,1),
                        (2,2,2,2,1,1,1,1), (2,2,1,1,2,1,1,1), (2,2,1,1,1,1,1,1),
                        (1,1,1,1,1,1,1,1)]
        @test energy ≈ [19.1, 18.2, 17.3, 16.4, 12.2, 11.3, 10.4, 10.2, 9.3, 8.4, 7.3, 6.4,
                        4.4]
        #  @test deg == [12, 12, 3, 6, 1]
    end

    let (basis, energy, deg) = subspace([1,2,3,4], 4, [0.1, 1.0, 10.0], 4)
        @test basis == [(4,4,4,4,3,3,3,1), (4,4,4,4,3,3,2,1), (4,4,4,4,3,3,1,1),
                        (4,4,4,4,3,2,2,1), (4,4,4,4,3,2,1,1), (4,4,4,4,3,1,1,1),
                        (4,4,4,4,2,2,2,1), (4,4,4,4,2,2,1,1), (4,4,4,4,2,1,1,1),
                        (4,4,4,4,1,1,1,1), (3,3,3,1,3,3,1,1), (3,3,3,1,3,2,1,1),
                        (3,3,3,1,3,1,1,1), (3,3,3,1,2,2,1,1), (3,3,3,1,2,1,1,1),
                        (3,3,3,1,1,1,1,1), (2,2,2,2,3,3,1,1), (2,2,2,2,3,2,1,1),
                        (2,2,2,2,3,1,1,1), (2,2,2,2,2,2,1,1), (2,2,2,2,2,1,1,1),
                        (2,2,2,2,1,1,1,1), (2,2,1,1,3,1,1,1), (2,2,1,1,2,1,1,1),
                        (2,2,1,1,1,1,1,1), (1,1,1,1,1,1,1,1)]
        @test energy ≈ [46.1, 37.1, 36.2, 28.1, 27.2, 26.3, 19.1, 18.2, 17.3, 16.4,
                        30.2, 21.2, 20.3, 12.2, 11.3, 10.4, 28.2, 19.2, 18.3, 10.2,
                        9.3, 8.4, 16.3, 7.3, 6.4, 4.4]
        #  @test deg == [12, 12, 3, 6, 1]
    end

    let (basis, energy, deg) = subspace([1,2,3,4,5], 5, [0.1, 1.0], 5)
        @test basis == [(5,5,5,5,5,2,2,2,2,1), (5,5,5,5,5,2,2,2,1,1), (5,5,5,5,5,2,2,1,1,1),
                        (5,5,5,5,5,2,1,1,1,1), (5,5,5,5,5,1,1,1,1,1), (4,4,4,4,1,2,2,2,1,1),
                        (4,4,4,4,1,2,2,1,1,1), (4,4,4,4,1,2,1,1,1,1), (4,4,4,4,1,1,1,1,1,1),
                        (3,3,3,2,2,2,2,2,1,1), (3,3,3,2,2,2,2,1,1,1), (3,3,3,2,2,2,1,1,1,1),
                        (3,3,3,2,2,1,1,1,1,1), (3,3,3,1,1,2,2,1,1,1), (3,3,3,1,1,2,1,1,1,1),
                        (3,3,3,1,1,1,1,1,1,1), (2,2,2,2,1,2,2,1,1,1), (2,2,2,2,1,2,1,1,1,1),
                        (2,2,2,2,1,1,1,1,1,1), (2,2,1,1,1,2,1,1,1,1), (2,2,1,1,1,1,1,1,1,1),
                        (1,1,1,1,1,1,1,1,1,1)]
        @test energy ≈ [29.1, 28.2, 27.3, 26.4, 25.5, 20.2, 19.3, 18.4, 17.5, 16.2, 15.3,
                        14.4, 13.5, 13.3, 12.4, 11.5, 11.3, 10.4, 9.5, 8.4, 7.5, 5.5]
        #  @test deg == [60, 60, 30, 30, 15, 10, 1]
    end
end
