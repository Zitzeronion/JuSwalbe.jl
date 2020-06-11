@testset "Minimal Example" begin
    @testset "One dimension" begin
        input, mom, force, dist = minimalsetup1d(10)
        @test isa(input, JuSwalbe.Inputconstants)
        @test isa(mom, JuSwalbe.Macroquant{Vector{Float64},Vector{Float64}})
        @test isa(force, JuSwalbe.Forces{Vector{Float64}})
        @test isa(dist, JuSwalbe.DistributionD1Q3)
        @test length(mom.height) == 10
        @test length(force.slip) == 10
        @test length(dist.f1) == 10
    end 
    @testset "Two dimension" begin
        input, mom, force, dist = minimalsetup2d(10,5)
        @test isa(input, JuSwalbe.Inputconstants)
        @test isa(mom, JuSwalbe.Macroquant{Matrix{Float64},JuSwalbe.Twovector{Matrix{Float64}}})
        @test isa(force, JuSwalbe.Forces{JuSwalbe.Twovector{Matrix{Float64}}})
        @test isa(dist, JuSwalbe.DistributionD2Q9)
        @test size(mom.height) == (10,5)
        @test size(force.slip.x) == (10,5)
        @test size(dist.f1) == (10,5)
    end 
end

@testset "Simple constructions" begin
    @testset "Simple Moment one dimension" begin
        for type in [Float16 Float32 Float64]
            mom = simplemoment1d(5, T=type)
            @test isa(mom, JuSwalbe.Macroquant{Vector{type}, Vector{type}})
            @test length(mom.height) == 5
            for i in mom.height
                @test i == type(1.0)
            end
            @test isa(mom.velocity, Vector{type})
            @test length(mom.velocity) == 5
            for i in mom.velocity
                @test i == type(0.1)
            end
            @test length(mom.pressure) == 5
            @test length(mom.energy) == 5
        end
    end
    @testset "Simple Twovector" begin
        for type in [Float16 Float32 Float64]
            xy = simpleTwovector(5,5; T=type)
            @test isa(xy, JuSwalbe.Twovector{Matrix{type}})
            @test size(xy.x) == (5,5)
            @test size(xy.y) == (5,5)
            for i in xy.x
                @test i == type(0.1)
            end
            for i in xy.y
                @test i == type(-0.1)
            end 
        end
    end
    @testset "Simple Moment two dimensions" begin
        for type in [Float16 Float32 Float64]
            mom = simplemoment2d(5,5, T=type)
            @test isa(mom, JuSwalbe.Macroquant{Matrix{type}, JuSwalbe.Twovector{Matrix{type}}})
            @test size(mom.height) == (5,5)
            for i in mom.height
                @test i == type(1.0)
            end
            @test isa(mom.velocity, JuSwalbe.Twovector{Matrix{type}})
            @test size(mom.velocity.x) == (5,5)
            for i in mom.velocity.x
                @test i == type(0.1)
            end
            @test size(mom.velocity.y) == (5,5)
            for i in mom.velocity.y
                @test i == type(-0.1)
            end
            @test size(mom.pressure) == (5,5)
            @test size(mom.energy) == (5,5)
        end
    end
    @testset "Simple D2Q9 Distribution" begin
        for type in [Float16 Float32 Float64]
            soldict = Dict(:f0 => type(1.0), :f1 => type(0.1), :f2 => type(-0.1), :f3 => type(0.01), :f4 => type(-0.01),
                           :f5 => type(0.2), :f6 => type(0.02), :f7 => type(0.02), :f8 => type(0.5))
            
                           xy = simpledistD2Q9(5,5; T=type)
            @test isa(xy, JuSwalbe.DistributionD2Q9{Matrix{type}})
            xydict = struct2dict(xy)
            for key in keys(xydict)
                @test size(xydict[key]) == (5,5)
                for element in xydict[key]
                    @test element == soldict[key]
                end
            end
            
        end
    end
end