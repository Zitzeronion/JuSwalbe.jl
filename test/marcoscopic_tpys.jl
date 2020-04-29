@testset "Moments" begin
    # Test the types and the sizes of the moments
    # For 64 bit floating numbers, one dimensional
    @testset "float64_1d" begin
        N = 10 
        moments64 = JuSwalbe.macroquant64_1d(zeros(Float64,N),zeros(Float64,N),zeros(Float64,N))
        @test typeof(moments64) == JuSwalbe.macroquant64_1d
        @test typeof(moments64.height) == Array{Float64,1}
        @test typeof(moments64.velocity) == Array{Float64,1}
        @test typeof(moments64.energy) == Array{Float64,1} 
        @test size(moments64.height, 1) == N
        @test size(moments64.velocity, 1) == N
        @test size(moments64.energy, 1) == N
    end
    # For 32 bit floating numbers, one dimensional
    @testset "float32_1d" begin
        M = 10
        moments32 = JuSwalbe.macroquant32_1d(zeros(Float32,M),zeros(Float32,M),zeros(Float32,M))
        @test typeof(moments32) == JuSwalbe.macroquant32_1d
        @test typeof(moments32.height) == Array{Float64,1}
        @test typeof(moments32.velocity) == Array{Float64,1}
        @test typeof(moments32.energy) == Array{Float64,1} 
        @test size(moments32.height, 1) == M
        @test size(moments32.velocity, 1) == M
        @test size(moments32.energy, 1) == M
    end
    # For 64 bit floating numbers, two dimensional
    @testset "float64_2d" begin
        N = 10
        M = 5
        Moments64 = JuSwalbe.macroquant64_2d(zeros(Float64,N,M),zeros(Float64,N,M),zeros(Float64,N,M))
        @test typeof(Moments64) == JuSwalbe.macroquant64_2d
        @test typeof(Moments64.height) == Array{Float64,2}
        @test typeof(moments64.velocity) == Array{Float64,2}
        @test typeof(moments64.energy) == Array{Float64,2} 
        @test size(moments64.height) == (N,M)
        @test size(moments64.height, 1) == N
        @test size(moments64.height, 2) == M
        @test size(moments64.velocity) == (N,M)
        @test size(moments64.velocity, 1) == N
        @test size(moments64.velocity, 1) == M
        @test size(moments64.energy) == (N,M)
        @test size(moments64.energy, 1) == N
        @test size(moments64.energy, 1) == M
    end
    # For 32 bit floating numbers, two dimensional
    @testset "float32_2d" begin
        N = 10
        M = 5
        Moments32 = JuSwalbe.macroquant32_2d(zeros(Float32,N,M),zeros(Float32,N,M),zeros(Float32,N,M))
        @test typeof(Moments32) == JuSwalbe.macroquant32_2d
        @test typeof(Moments32.height) == Array{Float32,2}
        @test typeof(moments32.velocity) == Array{Float32,2}
        @test typeof(moments32.energy) == Array{Float32,2} 
        @test size(moments32.height) == (N,M)
        @test size(moments32.height, 1) == N
        @test size(moments32.height, 2) == M
        @test size(moments32.velocity) == (N,M)
        @test size(moments32.velocity, 1) == N
        @test size(moments32.velocity, 1) == M
        @test size(moments32.energy) == (N,M)
        @test size(moments32.energy, 1) == N
        @test size(moments32.energy, 1) == M
    end
end