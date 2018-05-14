using RTFOOL
@static if VERSION < v"0.7.0-DEV.2005"
    using Base.Test
else
    using Test
end

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

    let a = Resource(0.25, [1,2]), b = Resource(0.25, [1,10]), c = Resource(0.25, [1,1,2])
        let aa = tensor((a,2)), H = [2, 3, 3, 4]
            @test length(aa) == 4
            @test size(aa) == (4,)
            @test aa.H == H
            @test aa.p ≈ RTFOOL.boltzmann(0.25, H)
        end

        let ab = tensor((a,1), (b,2)), H = [3,12,12,21,4,13,13,22]
            @test length(ab) == 8
            @test size(ab) == (8,)
            @test ab.H == H
            @test ab.p ≈ RTFOOL.boltzmann(0.25, H)
        end

        let abc = tensor((a,1), (b,2), (c,1)),
            H = [4, 4, 5, 13, 13, 14, 13, 13, 14, 22, 22, 23,
                 5, 5, 6, 14, 14, 15, 14, 14, 15, 23, 23, 24]
            @test length(abc) == 24
            @test size(abc) == (24,)
            @test abc.H == H
            @test abc.p ≈ RTFOOL.boltzmann(0.25, H)
        end
    end
end
