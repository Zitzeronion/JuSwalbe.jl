@testset "Write/Read Data and Checkpoints" begin
    @testset "Checkpoints D1Q3" begin
        typeset = [Float16 Float32 Float64]
        for T in typeset
            testdist = JuSwalbe.DistributionD1Q3(f0=ones(T, 10), f1=fill(T(0.1),10), f2=fill(T(-0.1),10))
            savecheckpoints(testdist)
            
            @test isfile("data/tmp/checkpoint_0.bson")
            readdist = loadcheckpoint("data/tmp/checkpoint_0.bson")
            @test testdist.f0 == readdist.f0
            @test testdist.f1 == readdist.f1
            @test testdist.f2 == readdist.f2

            rm("data", recursive=true)
        end
    end
    @testset "Checkpoints D2Q9" begin
        typeset = [Float16 Float32 Float64]
        for T in typeset
            testdist = simpledistD2Q9(5,5, T=T)
            @test isa(testdist, JuSwalbe.DistributionD2Q9{Matrix{T}}) 
            savecheckpoints(testdist)
            
            @test isfile("data/tmp/checkpoint_0.bson")
            readdist = loadcheckpoint("data/tmp/checkpoint_0.bson")
            @test testdist.f0 == readdist.f0
            @test testdist.f1 == readdist.f1
            @test testdist.f2 == readdist.f2
            @test testdist.f3 == readdist.f3
            @test testdist.f4 == readdist.f4
            @test testdist.f5 == readdist.f5
            @test testdist.f6 == readdist.f6
            @test testdist.f7 == readdist.f7
            @test testdist.f8 == readdist.f8

            rm("data", recursive=true)
        end
    end
    @testset "Checkpoint D1Q3" begin
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
    @testset "Checkpoints D2Q9" begin
        typeset = [Float16 Float32 Float64]
        for T in typeset
            testdist = simpledistD2Q9(5,5, T=T)
            @test isa(testdist, JuSwalbe.DistributionD2Q9{Matrix{T}}) 
            savecheckpoint(testdist)
            
            @test isfile("data/tmp/checkpoint_0.bson")
            readdist = loadcheckpoint("data/tmp/checkpoint_0.bson")
            @test testdist.f0 == readdist.f0
            @test testdist.f1 == readdist.f1
            @test testdist.f2 == readdist.f2
            @test testdist.f3 == readdist.f3
            @test testdist.f4 == readdist.f4
            @test testdist.f5 == readdist.f5
            @test testdist.f6 == readdist.f6
            @test testdist.f7 == readdist.f7
            @test testdist.f8 == readdist.f8

            rm("data", recursive=true)
        end
    end
    @testset "Save height data one dimension" begin
        for type in [Float16 Float32 Float64]
            mom = JuSwalbe.Macroquant(height = ones(type, 10), velocity = zeros(type, 10), pressure = zeros(type, 10), energy = zeros(type, 10))
            @test isa(mom, JuSwalbe.Macroquant{Vector{type}, Vector{type}})
            height2file(mom)
            @test isfile("data/tmp/height_0.bson")
            h = BSON.load("data/tmp/height_0.bson")
            @test length(h[:height]) == 10
            @test mom.height == h[:height]
        end
        rm("data", recursive=true)
    end
    @testset "Save height data two dimensions" begin
        for type in [Float16 Float32 Float64]
            mom = JuSwalbe.Macroquant(height = ones(type, (5,5)), velocity = simpleTwovector(5,5, T=type), pressure = zeros(type, (5,5)), energy = zeros(type, (5,5)))
            @test isa(mom, JuSwalbe.Macroquant{Matrix{type}, JuSwalbe.Twovector{Matrix{type}}})
            height2file(mom)
            @test isfile("data/tmp/height_0.bson")
            h = BSON.load("data/tmp/height_0.bson")
            @test size(h[:height]) == (5,5)
            @test mom.height == h[:height]
        end
        rm("data", recursive=true)
    end
    @testset "Save velocity data one dimension" begin
        for type in [Float16 Float32 Float64]
            mom = JuSwalbe.Macroquant(height = ones(type, 10), velocity = fill(type(0.1), 10), pressure = zeros(type, 10), energy = zeros(type, 10))
            @test isa(mom, JuSwalbe.Macroquant{Vector{type}, Vector{type}})
            velocity2file(mom)
            @test isfile("data/tmp/velocity_0.bson")
            v = BSON.load("data/tmp/velocity_0.bson")
            @test length(v[:velocity]) == 10
            @test mom.velocity == v[:velocity]
        end
        rm("data", recursive=true)
    end
    @testset "Save velocity data two dimensions" begin
        for type in [Float16 Float32 Float64]
            mom = JuSwalbe.Macroquant(height = ones(type, (5,5)), velocity = simpleTwovector(5,5, T=type), pressure = zeros(type, (5,5)), energy = zeros(type, (5,5)))
            @test isa(mom, JuSwalbe.Macroquant{Matrix{type}, JuSwalbe.Twovector{Matrix{type}}})
            velocity2file(mom)
            @test isfile("data/tmp/velocity_0.bson")
            v = BSON.load("data/tmp/velocity_0.bson")
            @test size(v[:velocity_x]) == (5,5)
            @test size(v[:velocity_y]) == (5,5)
            @test mom.velocity.x == v[:velocity_x]
            @test mom.velocity.y == v[:velocity_y]
        end
        rm("data", recursive=true)
    end
    @testset "Save velocity and height data one dimension" begin
        for type in [Float16 Float32 Float64]
            mom = JuSwalbe.Macroquant(height = ones(type, 10), velocity = fill(type(0.1), 10), pressure = zeros(type, 10), energy = zeros(type, 10))
            @test isa(mom, JuSwalbe.Macroquant{Vector{type}, Vector{type}})
            velocityandheight2file(mom)
            @test isfile("data/tmp/height_velocity_0.bson")
            hv = BSON.load("data/tmp/height_velocity_0.bson")
            @test length(hv[:velocity]) == 10
            @test length(hv[:height]) == 10
            @test mom.velocity == hv[:velocity]
            @test mom.height == hv[:height]
        end
        rm("data", recursive=true)
    end
    @testset "Save velocity and height data two dimensions" begin
        for type in [Float16 Float32 Float64]
            mom = JuSwalbe.Macroquant(height = ones(type, (5,5)), velocity = simpleTwovector(5,5, T=type), pressure = zeros(type, (5,5)), energy = zeros(type, (5,5)))
            @test isa(mom, JuSwalbe.Macroquant{Matrix{type}, JuSwalbe.Twovector{Matrix{type}}})
            velocityandheight2file(mom)
            @test isfile("data/tmp/height_velocity_0.bson")
            hv = BSON.load("data/tmp/height_velocity_0.bson")
            @test size(hv[:height]) == (5,5)
            @test size(hv[:velocity_x]) == (5,5)
            @test size(hv[:velocity_y]) == (5,5)
            @test mom.height == hv[:height]
            @test mom.velocity.x == hv[:velocity_x]
            @test mom.velocity.y == hv[:velocity_y]
        end
        rm("data", recursive=true)
    end
end
