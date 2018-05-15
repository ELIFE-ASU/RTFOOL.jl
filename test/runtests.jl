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

    @testset "Pure System" begin
        let Hm = [1, 10], Hw = [10, 20]
            @test_throws DimensionMismatch pure_system(Hm, [1], Hw, [1,0])
            @test_throws DimensionMismatch pure_system(Hm, [1,0], Hw, [1])
            @test_throws ArgumentError pure_system(Hm, [2,0], Hw, [1,0])

            let (r, Nm, Nw) = pure_system(Hm, [2,0], Hw, [2,0])
                @test Nm == 2
                @test Nw == 2
                @test r.p[1] == 1
            end

            let (r, Nm, Nw) = pure_system(Hm, [1,1], Hw, [1,1])
                @test Nm == 3
                @test Nw == 2
                @test r.p[14] == 1
            end
        end
    end
end

@testset "Degeneracies" begin
    let (deg, num_swaps) = RTFOOL.degeneracies([2,3,3,4])
        @test deg == Dict(2 => [3])
        @test num_swaps == 1
    end

    let (deg, num_swaps) = RTFOOL.degeneracies([3,4,4,5,4,5,5,6])
        @test deg == Dict(2 => [3,5], 3 => [5], 4 => [6,7], 6=>[7])
        @test num_swaps == 6
    end

    let (deg, num_swaps) = RTFOOL.degeneracies(tensor((Resource([1,0], [1,2]),3)))
        @test deg == Dict(2 => [3,5], 3 => [5], 4 => [6,7], 6=>[7])
        @test num_swaps == 6
    end
end

@testset "Context" begin
    let β = 0.5, Hm = [1,2], Hw = [1,2,3],
        system = pure_system(Hm, [0,1], Hw, [0,1,0])[1],
        ctx = Context(β, Hm, 2, Hw, 1, system)

        @test ctx.β == β

        @test ctx.monomer ≈ Resource(β, [1,2])
        @test ctx.water ≈ Resource(β, [1,2,3])

        @test ctx.bath ≈ Resource(β, [3,4,5,4,5,6,4,5,6,5,6,7])
        @test ctx.Nm == 2
        @test ctx.Nw == 1

        @test ctx.system == ctx.system
        @test length(ctx.H) == 144
        @test ctx.H[1] == 6
        @test ctx.H[end] == 14

        @test ctx.num_swaps == 1802
    end

    let β = 0.5, Hm = [1,2], Hw = [1,2,3], ms = [0,1], ws = [0,1,0],
        (system, Nm, Nw) = pure_system(Hm, ms, Hw, ws),
        ctx = Context(β, Hm, ms, Hw, ws)

        @test ctx.β == β

        @test ctx.monomer ≈ Resource(β, [1,2])
        @test ctx.water ≈ Resource(β, [1,2,3])

        @test ctx.bath ≈ Resource(β, [3,4,5,4,5,6,4,5,6,5,6,7])
        @test ctx.Nm == Nm
        @test ctx.Nw == Nw

        @test ctx.system == ctx.system
        @test length(ctx.H) == 144
        @test ctx.H[1] == 6
        @test ctx.H[end] == 14

        @test ctx.num_swaps == 1802
    end
end
