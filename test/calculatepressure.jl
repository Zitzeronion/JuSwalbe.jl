@testset "Pressure" begin
    @testset "Laplacian_2d" begin
        N, M = (4,4)
        moment = JuSwalbe.macroquant64_2d(ones(N,M), JuSwalbe.velocity64_2d(fill(0.1, (N,M)),fill(0.2, (N,M))), zeros(N,M))
        # Fill the array a little bit more demanding, otherwise the laplacian would be zero
        for i in 1:size(moment.height,1) 
            for j in 1:size(moment.height,2)
                moment.height[i,j] = i+j
            end
        end
        # Make the laplace calculation
        test_lap = Δh(moment)
        # Test the dimensions
        @test size(test_lap) == (N,M)
        result = [8.0 4.0 4.0 0.0; 4.0 0.0 0.0 -4.0; 4.0 0.0 0.0 -4.0; 0.0 -4.0 -4.0 -8.0]
        # Compare the results
        for j in 1:M, i in 1:N
            @test result[i,j] == test_lap[i,j]
        end
    end

    @testset "Disjoining_2d_macroquant" begin
        N, M = (4,4)
        moment = JuSwalbe.macroquant64_2d(reshape(1:(N*M), N, M), JuSwalbe.velocity64_2d(fill(0.1, (N,M)),fill(0.2, (N,M))), zeros(N,M))
        
        # Make the laplace calculation
        test_disj = Π(moment)
        # Test the dimensions
        @test size(test_disj) == (N,M)
        result = [-1.6082e-5   -1.28656e-7  -2.20603e-8  -7.31997e-9;
                  -2.01025e-6  -7.44536e-8  -1.6082e-8   -5.86078e-9;
                  -5.95628e-7  -4.68862e-8  -1.20826e-8  -4.76503e-9;
                  -2.51281e-7  -3.14101e-8  -9.30669e-9  -3.92626e-9]
        # Compare the results
        for j in 1:M, i in 1:N
            @test result[i,j] ≈ test_disj[i,j] atol = 1e-10
        end
    end
    @testset "Disjoining_2d_array" begin
        N, M = (4,4)
        height = ones(4,4)
        # Make the disjoining pressure calculation
        test_disj = Π(height, exponents=[3,2],γ=1.0,θ=0.5)
        # Test the dimensions
        @test size(test_disj) == (N,M)
        result = 2.0/0.1 * (0.1^3 - 0.1^2)
        # Compare the results
        for j in 1:M, i in 1:N
            @test test_disj[i,j] == result
        end
    end
    
    @testset "Filmpressure_nothing_2d_moment" begin
        N, M = (4,4)
        moment = JuSwalbe.macroquant64_2d(reshape(collect(1.0:1.0:(N*M)),N,M), JuSwalbe.velocity64_2d(zeros(N,M),zeros(N,M)),zeros(N,M))
        # Make the laplace calculation
        # Make the laplace calculation
        test_pressure = pressure(moment, γ=1.0)
        # Test the dimensions
        @test size(test_pressure) == (N,M)
        result = 0.0
        # Compare the results
        for j in 1:M, i in 1:N
            @test test_pressure[i,j] == result
        end
    end
    @testset "Filmpressure_nodisjoining_2d_moment" begin
        N, M = (4,4)
        moment = JuSwalbe.macroquant64_2d(reshape(collect(1.0:1.0:(N*M)),N,M), JuSwalbe.velocity64_2d(zeros(N,M),zeros(N,M)),zeros(N,M))
        # Make the laplace calculation
        test_pressure = pressure(moment, θ=0.0, γ=0.1)
        # Test the dimensions
        @test size(test_pressure) == (N,M)
        # Compare the results
        @test test_pressure == -0.1*Δh(moment)
    end
    @testset "Filmpressure_2d_moment" begin
        N, M = (4,4)
        moment = JuSwalbe.macroquant64_2d(reshape(collect(1.0:1.0:(N*M)),N,M), JuSwalbe.velocity64_2d(zeros(N,M),zeros(N,M)),zeros(N,M))
        # Make the laplace calculation
        test_pressure = pressure(moment, θ=1.0/9.0, γ=0.1)
        # Test the dimensions
        @test size(test_pressure) == (N,M)
        # Compare the results
        @test test_pressure == -0.1*Δh(moment) .+ Π(moment, γ=0.1)
    end
    @testset "Filmpressure_nolaplacian_2d_moment" begin
        N, M = (4,4)
        moment = JuSwalbe.macroquant64_2d(ones(N,M), JuSwalbe.velocity64_2d(zeros(N,M),zeros(N,M)),zeros(N,M))
        # Make the laplace calculation
        test_pressure = pressure(moment, γ=1.0)
        # Test the dimensions
        @test size(test_pressure) == (N,M)
        
        # Compare the results
        @test test_pressure == Π(moment, γ=1.0)
        
    end
    @testset "Filmpressure_nothing_2d_array" begin
        N, M = (4,4)
        height = fill(0.1, (N,M))
        # Make the laplace calculation
        test_pressure = pressure(height, γ=1.0)
        # Test the dimensions
        @test size(test_pressure) == (N,M)
        result = 0.0
        # Compare the results
        for j in 1:M, i in 1:N
            @test test_pressure[i,j] == result
        end
    end
    @testset "Filmpressure_nolaplacian_2d_array" begin
        N, M = (4,4)
        height = ones(N,M)
        # Make the laplace calculation
        test_pressure = pressure(height, γ=1.0)
        # Test the dimensions
        @test size(test_pressure) == (N,M)
        
        # Compare the results
        @test test_pressure == Π(height, γ=1.0)
        
    end
    @testset "Filmpressure_nodisjoining_2d_array" begin
        N, M = (4,4)
        height = reshape(collect(1.0:1.0:(N*M)),N,M)
        # Make the laplace calculation
        test_pressure = pressure(height, θ=0.0, γ=0.1)
        # Test the dimensions
        @test size(test_pressure) == (N,M)
        # Compare the results
        @test test_pressure == -0.1*Δh(height)
    end
    @testset "Filmpressure_2d_array" begin
        N, M = (4,4)
        height = reshape(collect(1.0:1.0:(N*M)),N,M)
        # Make the laplace calculation
        test_pressure = pressure(height, θ=1.0/9.0, γ=0.1)
        # Test the dimensions
        @test size(test_pressure) == (N,M)
        # Compare the results
        @test test_pressure == -0.1*Δh(height) .+ Π(height, γ=0.1)
    end
end