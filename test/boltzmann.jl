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
