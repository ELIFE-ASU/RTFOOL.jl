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
end
