@testset "Monotones" begin
    @test_throws MethodError relative_entropy([], [])
    @test_throws ErrorException relative_entropy([1.0], [0.5, 0.5])
    @test_throws ErrorException relative_entropy([0.5, 0.5], [1.0])
    @test_throws ErrorException relative_entropy([0.5], [1.0])
    @test_throws ErrorException relative_entropy([1.0], [0.5])

    @test relative_entropy([1.0, 0.0], [0.5, 0.5]) ≈ log(2.0)
    @test relative_entropy([0.5, 0.5], [0.5, 0.5]) ≈ 0.0
    @test isinf(relative_entropy([0.5, 0.5], [1.0, 0.0]))
    @test relative_entropy([0.25, 0.75], [0.5, 0.5]) ≈ 0.75*log(3) - log(2)
end
