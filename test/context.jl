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

let rng = MersenneTwister(2018)
    @testset "Timestep" begin
        let ctx = Context(0.5, StateSpace([1,2], 4, [0.1, 1.0], 4), [1, 0, 0])
            re = relative_entropy(ctx.system_state, ctx.bath_state)
            for _ in 1:10000
                timestep(rng, ctx)
                r = relative_entropy(ctx.system_state, ctx.bath_state)
                @test isapprox(r, re; atol=1e-9) || r < re
                re = r
            end
            @test ctx.system_state ≈ ctx.bath_state
        end
        let ctx = Context(0.5, StateSpace(1:4, 4, [0.1, 1.0], 4), [1,0,0,0,0])
            re = relative_entropy(ctx.system_state, ctx.bath_state)
            for _ in 1:1000000
                timestep(rng, ctx)
                r = relative_entropy(ctx.system_state, ctx.bath_state)
                @test isapprox(r, re; atol=1e-9) || r < re
                re = r
            end
            @test ctx.system_state ≈ ctx.bath_state
        end
    end
end
