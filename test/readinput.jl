@testset "Read input" begin
    testinput = JuSwalbe.readinput("test.txt")
    @test isa(testinput, JuSwalbe.Inputconstants)
    @test testinput.lx == 10
    @test testinput.ly == 5
    @test testinput.maxruntime == 1000
    @test testinput.dumping == 100
    @test testinput.τ == 1.0
    @test testinput.gravity == 0.0
    @test testinput.γ == 0.01
    @test testinput.δ == 1.0
end

@testset "Find argument" begin
    arr = ["hello" 15; "my" 0.077; "wonderful" 10534; "world" 0.4231]
    @test findargument(arr, "wonderful") == 10534;
    @test findargument(arr, "hello") == 15;
    @test findargument(arr, "world") == 0.4231;
    @test findargument(arr, "my") == 0.077;
    # Returns nothing if the string can't be found
    @test findargument(arr, "lol") == nothing
    # Has to be a string
    @test_throws MethodError findargument(arr, 2)
end