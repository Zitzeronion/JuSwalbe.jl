@testset "call" begin
    @test call(1) == "Hello world"
    @test call(2) == "Hello world"
    @test_throws MethodError call("hello") == "Hello world"
end