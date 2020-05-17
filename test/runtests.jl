using JuSwalbe
using Test, DrWatson, DelimitedFiles, Images

@testset "JuSwalbe.jl" begin
    include("call.jl")
    include("macroscopic_typs.jl")
    include("distribution_types.jl")
    include("readinput.jl")
    include("equilibriumcalculation.jl")
    include("calculatepressure.jl")
    include("calcmoments.jl")
    include("forcing.jl")
    include("collide.jl")
end

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
