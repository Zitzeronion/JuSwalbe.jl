@testset "Write/Read Data and Checkpoints" begin
    @testset "Checkpoints" begin
        typeset = [Float16 Float32 Float64]
        for T in typeset
            testdist = JuSwalbe.DistributionD1Q3(f0=ones(T, 10), f1=fill(T(0.1),10), f2=fill(T(-0.1),10))
            savecheckpoint(testdist)
            
            @test isfile("data/tmp/checkpoint_0.bson")
            readdist = loadcheckpoint("data/tmp/checkpoint_0.bson")
            @test testdist.f0 == readdist.f0
            @test testdist.f1 == readdist.f1
            @test testdist.f2 == readdist.f2

            rm("data", recursive=true)
        end
    end
end
