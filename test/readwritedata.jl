@testset "Write/Read Data and Checkpoints" begin
    @testset "Checkpoints" begin
        typeset = [Float16 Float32 Float64]
        for T in typeset
            testdist = JuSwalbe.DistributionD1Q3(f0=ones(10), f1=fill(T(0.1),10), f2=fill(T(-0.1),10))
            savecheckpoint(testdist)
            @test isfile("")