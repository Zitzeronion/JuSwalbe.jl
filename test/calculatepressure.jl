@testset "Pressure" begin
    tolerances = Dict(Float16 => 1e-3, Float32 => 1e-5, Float64 => 1e-7, Real => 1e-7)
    @testset "Laplacian two dimensions" begin
        N, M = (4,4)
        typelist = [Float16 Float32 Float64 Real]
        for type in typelist
            testheight = zeros(type, (N,M))
            velocities = JuSwalbe.Twovector{Matrix{type}}(zeros(type, (N,M)),zeros(type, (N,M)))
            moment = JuSwalbe.Macroquant(ones(type, (N,M)), velocities, zeros(type, (N,M)), zeros(type, (N,M)))
            # Fill the array a little bit more demanding, otherwise the laplacian would be zero
            for i in 1:size(moment.height,1) 
                for j in 1:size(moment.height,2)
                    moment.height[i,j] = i+j
                    testheight[i,j] = i+j
                end
            end
            # Make the laplace calculation
            test_lap1 = Δh(moment)
            @test isa(test_lap1, Array{type, 2})
            test_lap2 = Δh(testheight)
            @test isa(test_lap2, Array{type, 2})
            # Test the dimensions
            @test size(test_lap1) == (N,M)
            @test size(test_lap2) == (N,M)
            result = [8.0 4.0 4.0 0.0; 4.0 0.0 0.0 -4.0; 4.0 0.0 0.0 -4.0; 0.0 -4.0 -4.0 -8.0]
            resultcorrected = zeros(type, (N,M))
            for (ind, ele) in enumerate(result)
                resultcorrected[ind] = result[ind]
            end
            # Compare the results
            for j in 1:M, i in 1:N
                @test result[i,j] == test_lap1[i,j]
                @test result[i,j] == test_lap2[i,j]
            end
        end
    end

    @testset "Laplacian one dimension" begin
        N = 8
        typelist = [Float16 Float32 Float64 Real]
        for type in typelist
            testheight = ones(type, N)
            moment = JuSwalbe.Macroquant(ones(type, N), zeros(type, N), zeros(type, N), zeros(type, N))
            # Fill the array a little bit more demanding, otherwise the laplacian would be zero
            for i in 1:length(moment.height) 
                moment.height[i] = i*i
                testheight[i] = i*i
            end
            # Make the laplace calculation
            test_lap1 = Δh(moment)
            @test isa(test_lap1, Array{type, 1})
            test_lap2 = Δh(testheight)
            @test isa(test_lap2, Array{type, 1})
            # Test the dimensions
            @test length(test_lap1) == N
            @test length(test_lap2) == N
            result = [66.0 2.0 2.0 2.0 2.0 2.0 2.0 -78.0]
            resultcorrected = zeros(type, N)
            for (ind, ele) in enumerate(result)
                resultcorrected[ind] = result[ind]
            end
            # Compare the results
            for i in 1:N
                @test result[i] == test_lap1[i]
                @test result[i] == test_lap2[i]
            end
        end
    end

    @testset "Disjoining pressure two dimensions" begin
        N, M = (4,4)

        typelist = [Float16 Float32 Float64]
        for type in typelist
            testheight = ones(type, (N,M)) .* reshape(collect(1:N*M),N,M)
            moment = JuSwalbe.Macroquant(testheight, JuSwalbe.Twovector(ones(type, (N,M)),ones(type, (N,M))), zeros(type, (N,M)), zeros(type, (N,M)))
            
            # Make the laplace calculation
            test_disj1 = Π(moment)
            @test isa(test_disj1, Array{type, 2})
            test_disj2 = Π(testheight)
            @test isa(test_disj1, Array{type, 2})
            # Test the dimensions
            @test size(test_disj1) == (N,M)
            @test size(test_disj2) == (N,M)
            result = [-1.6082e-5   -1.28656e-7  -2.20603e-8  -7.31997e-9;
                    -2.01025e-6  -7.44536e-8  -1.6082e-8   -5.86078e-9;
                    -5.95628e-7  -4.68862e-8  -1.20826e-8  -4.76503e-9;
                    -2.51281e-7  -3.14101e-8  -9.30669e-9  -3.92626e-9]
            # Compare the results
            resultcorrected = zeros(type, (N,M))
            for (ind, ele) in enumerate(result)
                resultcorrected[ind] = result[ind]
            end
            for i in 1:M, j in 1:N
                @test result[i,j] - test_disj1[i,j] ≈ 0 atol = tolerances[type]
                @test result[i,j] - test_disj2[i,j] ≈ 0 atol = tolerances[type]
            end
        end
    end
    @testset "Disjoining pressure one dimension" begin
        N = 4
        typelist = [Float16 Float32 Float64]
        for type in typelist
            testheight = ones(type, N) .* collect(1:N)
            moment = JuSwalbe.Macroquant(testheight, ones(type, N), zeros(type, N), zeros(type, N))
            
            # Make the laplace calculation
            test_disj1 = Π(moment)
            @test isa(test_disj1, Array{type, 1})
            test_disj2 = Π(testheight)
            @test isa(test_disj2, Array{type, 1})
            # Test the dimensions
            @test length(test_disj1) == N
            @test length(test_disj2) == N
            result = [-1.603e-5 -2.0e-6 -6.0e-7 -2.4e-7]

            # Compare the results
            resultcorrected = zeros(type, N)
            for (ind, ele) in enumerate(result)
                resultcorrected[ind] = result[ind]
            end
            for i in 1:N
                @test result[i] - test_disj1[i] ≈ 0 atol = tolerances[type]
                @test result[i] - test_disj2[i] ≈ 0 atol = tolerances[type]
            end
        end
    end
    
    @testset "Filmpressure nothing two dimensions" begin
        N, M = (4,4)
        typelist = [Float16 Float32 Float64]
        for type in typelist
            testheight = ones(type, (N,M))
            moment = JuSwalbe.Macroquant(testheight, JuSwalbe.Twovector(zeros(type ,(N,M)),zeros(type ,(N,M))), zeros(type ,(N,M)), zeros(type ,(N,M)))
            # Make the pressure calculation
            test_pressure1 = pressure(moment, γ=type(1.0), θ=zeros(type,(1,1)))
            test_pressure2 = pressure(moment, γ=type(1.0), θ=zeros(type,(1,1)))
            # Test the dimensions
            @test isa(test_pressure1, Array{type, 2})
            @test isa(test_pressure2, Array{type, 2})
            @test size(test_pressure1) == (N,M)
            @test size(test_pressure2) == (N,M)
            result = type(0.0)
            # Compare the results
            for j in 1:M, i in 1:N
                @test test_pressure1[i,j] == result
                @test test_pressure1[i,j] == result
            end
        end
    end

    @testset "Filmpressure nothing one dimension" begin
        N = 4
        typelist = [Float16 Float32 Float64]
        for type in typelist
            testheight = ones(type, N)
            moment = JuSwalbe.Macroquant(testheight, zeros(type ,N), zeros(type ,N), zeros(type ,N))
            # Make the pressure calculation
            test_pressure1 = pressure(moment, γ=type(1.0), θ=zeros(type, 1))
            test_pressure2 = pressure(moment, γ=type(1.0), θ=zeros(type, 1))
            # Test the dimensions
            @test isa(test_pressure1, Array{type, 1})
            @test isa(test_pressure2, Array{type, 1})
            @test length(test_pressure1) == N
            @test length(test_pressure2) == N
            result = type(0.0)
            # Compare the results
            for i in 1:N
                @test test_pressure1[i] == result
                @test test_pressure1[i] == result
            end
        end
    end

    @testset "Filmpressure nodisjoining two dimensional" begin
        N, M = (4,4)
        typelist = [Float16 Float32 Float64]
        for type in typelist
            testheight = ones(type, (N,M)) .* reshape(collect(1:N*M),N,M)
            moment = JuSwalbe.Macroquant(testheight, JuSwalbe.Twovector(zeros(type ,(N,M)),zeros(type ,(N,M))), zeros(type ,(N,M)), zeros(type ,(N,M)))
            
            test_pressure1 = pressure(moment, θ=zeros(type, (1,1)), γ=type(0.1))
            test_pressure2 = pressure(testheight, θ=zeros(type, (1,1)), γ=type(0.1))
            # Test the dimensions
            @test isa(test_pressure1, Array{type,2})
            @test isa(test_pressure2, Array{type,2})
            @test size(test_pressure1) == (N,M)
            @test size(test_pressure2) == (N,M)
            # Compare the results
            @test test_pressure1 - (type(-0.1)*Δh(moment)) ≈ 0 tolerances[type]
            @test test_pressure2 - (type(-0.1)*Δh(testheight)) ≈ 0 tolerances[type]
        end
    end

    @testset "Filmpressure nodisjoining one dimensional" begin
        N = 4
        typelist = [Float16 Float32 Float64]
        for type in typelist
            testheight = ones(type, N) .* collect(1:N)
            moment = JuSwalbe.Macroquant(testheight, zeros(type ,N), zeros(type ,N), zeros(type ,N))
            
            test_pressure1 = pressure(moment, θ=zeros(type, 1), γ=type(0.1))
            test_pressure2 = pressure(testheight, θ=zeros(type, 1), γ=type(0.1))
            # Test the dimensions
            @test isa(test_pressure1, Array{type,1})
            @test isa(test_pressure2, Array{type,1})
            @test length(test_pressure1) == N
            @test length(test_pressure2) == N
            # Compare the results
            @test test_pressure1 - (type(-0.1)*Δh(moment)) ≈ 0 tolerances[type]
            @test test_pressure2 - (type(-0.1)*Δh(testheight)) ≈ 0 tolerances[type]
        end
    end
    
    @testset "Filmpressure NoLaplacian two dimensional" begin
        N, M = (4,4)
        typelist = [Float16 Float32 Float64]
        for type in typelist
            testheight = ones(type, (N,M)) .* reshape(collect(1:N*M),N,M)
            moment = JuSwalbe.Macroquant(testheight, JuSwalbe.Twovector(zeros(type ,(N,M)),zeros(type ,(N,M))), zeros(type ,(N,M)), zeros(type ,(N,M)))
            # Make the laplace calculation
            test_pressure1 = pressure(moment, γ=1.0)
            test_pressure2 = pressure(testheight, γ=1.0)
            # Test the dimensions
            @test isa(test_pressure1, Array{type,2})
            @test isa(test_pressure2, Array{type,2})
            @test size(test_pressure1) == (N,M)
            @test size(test_pressure2) == (N,M)
            # Compare the results
            for i in 1:N, j in 1:M
                @test test_pressure1[i,j] - Π(moment, γ=1.0)[i,j] ≈ 0 tolerances[type]
                @test test_pressure2[i,j] - Π(moment, γ=1.0)[i,j] ≈ 0 tolerances[type]
        end
        
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