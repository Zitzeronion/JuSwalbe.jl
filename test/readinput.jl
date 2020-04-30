@testset "Read input" begin
    testinput = JuSwalbe.readinput("test.txt")
    @test isa(testinput, JuSwalbe.inputconstants)
    @test testinput.lx == 10
    @test testinput.ly == 5
    @test testinput.maxruntime == 1000
    @test testinput.dumping == 100
    @test testinput.gravity == 0.0
    @test testinput.γ == 0.01
    @test testinput.δ == 1.0
end