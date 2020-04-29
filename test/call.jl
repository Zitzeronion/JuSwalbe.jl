@testset "call" begin
    @test call(1) == "Hello world"
    @test call(2) == "Hello world"
    @test_throws MethodError call(2) == 1
end