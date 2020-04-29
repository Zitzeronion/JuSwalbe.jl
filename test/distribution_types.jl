@testset "Distributions" begin
    # Test the types and the sizes of the distributions
    # For 64 bit floating numbers, one dimensional
    @testset "float64_1d" begin
        N = 10 
        dist64 = JuSwalbe.dist64_1d(zeros(Float64,N),zeros(Float64,N),zeros(Float64,N))
        @test typeof(dist64) == JuSwalbe.dist64_1d
        @test typeof(dist64.f0) == Array{Float64,1}
        @test typeof(dist64.f1) == Array{Float64,1}
        @test typeof(dist64.f2) == Array{Float64,1} 
        @test size(dist64.f0, 1) == N
        @test size(dist64.f1, 1) == N
        @test size(dist64.f2, 1) == N
    end
    # For 32 bit floating numbers, one dimensional
    @testset "float32_1d" begin
        M = 10
        dist32 = JuSwalbe.dist32_1d(zeros(Float32,M),zeros(Float32,M),zeros(Float32,M))
        @test typeof(dist32) == JuSwalbe.dist32_1d
        @test typeof(dist32.f0) == Array{Float32,1}
        @test typeof(dist32.f1) == Array{Float32,1}
        @test typeof(dist32.f2) == Array{Float32,1} 
        @test size(dist32.f0, 1) == M
        @test size(dist32.f1, 1) == M
        @test size(dist32.f2, 1) == M
    end
    # For 64 bit floating numbers, two dimensional
    @testset "float64_2d" begin
        N = 10
        M = 5
        dist64 = JuSwalbe.dist64_2d(zeros(Float64, N, M),zeros(Float64, N, M),zeros(Float64, N, M),zeros(Float64, N, M),zeros(Float64, N, M),zeros(Float64, N, M),zeros(Float64, N, M),zeros(Float64, N, M),zeros(Float64, N, M))
        @test typeof(dist64) == JuSwalbe.dist64_2d
        test_dict = struct2dict(dist64)
        for key in keys(test_dict)
            @test typeof(test_dict[key]) == Array{Float64,2}
            @test size(test_dict[key]) == (N,M)
            @test size(test_dict[key], 1) == N
            @test size(test_dict[key], 2) == M    
        end
    end
    # For 32 bit floating numbers, two dimensional
    @testset "float32_2d" begin
        N = 10
        M = 5
        dist32 = JuSwalbe.dist32_2d(zeros(Float32, N, M),zeros(Float32, N, M),zeros(Float32, N, M),zeros(Float32, N, M),zeros(Float32, N, M),zeros(Float32, N, M),zeros(Float32, N, M),zeros(Float32, N, M),zeros(Float32, N, M))
        @test typeof(dist32) == JuSwalbe.dist32_2d
        test_dict = struct2dict(dist32)
        for key in keys(test_dict)
            @test typeof(test_dict[key]) == Array{Float32,2}
            @test size(test_dict[key]) == (N,M)
            @test size(test_dict[key], 1) == N
            @test size(test_dict[key], 2) == M    
        end
    end
end