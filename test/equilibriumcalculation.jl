@testset "Equilibria" begin
    @testset "1d_64_nogravity_novelocity" begin
        N = 10
        # Generate test moments and a test distribution
        testmoments = JuSwalbe.macroquant64_1d(ones(Float64,N),zeros(Float64,N),zeros(Float64,N))
        equilibria = calc_equilibrium_distribution(testmoments)
        # Push them in a dictonary to easily loop through the keys (fields)
        test_dict = struct2dict(equilibria)
        # Test types and lengths
        @test isa(equilibria, JuSwalbe.dist64_1d)
        @test length(keys(test_dict)) == 3
        # Test the 
        for key in keys(test_dict)
            @test isa(test_dict[key], Array{Float64,1})
            @test length(test_dict[key]) == N
            # Check if the values are actually consistent
            if key == :f0
                for value in test_dict[key]
                    @test value == 1.0
                end
            else
                for value in test_dict[key]
                    @test value == 0.0
                end
            end
        end
    end
    @testset "1d_64_gravity_novelocity" begin
        N = 10
        # Generate test moments and a test distribution
        testmoments = JuSwalbe.macroquant64_1d(ones(Float64,N),zeros(Float64,N),zeros(Float64,N))
        equilibria = calc_equilibrium_distribution(testmoments, 0.1)
        # Push them in a dictonary to easily loop through the keys (fields)
        test_dict = struct2dict(equilibria)
        # Test types and lengths
        @test isa(equilibria, JuSwalbe.dist64_1d)
        @test length(keys(test_dict)) == 3
        # Test the 
        for key in keys(test_dict)
            @test isa(test_dict[key], Array{Float64,1})
            @test length(test_dict[key]) == N
            # Check if the values are actually consistent
            if key == :f0
                for value in test_dict[key]
                    @test value == 0.95
                end
            else
                for value in test_dict[key]
                    @test value == 1.0/40.0
                end
            end
        end
    end
    @testset "1d_64_gravity_velocity" begin
        N = 10
        # Generate test moments and a test distribution
        testmoments = JuSwalbe.macroquant64_1d(ones(Float64,N),fill(0.1,N),zeros(Float64,N))
        equilibria = calc_equilibrium_distribution(testmoments, 0.1)
        # Push them in a dictonary to easily loop through the keys (fields)
        test_dict = struct2dict(equilibria)
        # Test types and lengths
        @test isa(equilibria, JuSwalbe.dist64_1d)
        @test length(keys(test_dict)) == 3
        # Test the 
        for key in keys(test_dict)
            @test isa(test_dict[key], Array{Float64,1})
            @test length(test_dict[key]) == N
            # Check if the values are actually consistent
            if key == :f0
                for value in test_dict[key]
                    @test value == 1 - 1.0/20.0 - 0.01
                end
            elseif key == :f1
                for value in test_dict[key]
                    @test value == 1.0/40.0 + 1.0/20.0 + 1.0/200.0
                end
            elseif key == :f2
                for value in test_dict[key]
                    @test value == 1.0/40.0 - 1.0/20.0 + 1.0/200.0
                end
            end
        end
    end
    @testset "1d_64_nogravity_velocity" begin
        N = 10
        # Generate test moments and a test distribution
        testmoments = JuSwalbe.macroquant64_1d(ones(Float64,N),fill(0.1,N),zeros(Float64,N))
        equilibria = calc_equilibrium_distribution(testmoments)
        # Push them in a dictonary to easily loop through the keys (fields)
        test_dict = struct2dict(equilibria)
        # Test types and lengths
        @test isa(equilibria, JuSwalbe.dist64_1d)
        @test length(keys(test_dict)) == 3
        # Test the 
        for key in keys(test_dict)
            @test isa(test_dict[key], Array{Float64,1})
            @test length(test_dict[key]) == N
            # Check if the values are actually consistent
            if key == :f0
                for value in test_dict[key]
                    @test value == 1 - 0.01
                end
            elseif key == :f1
                for value in test_dict[key]
                    @test value ≈ 1.0/20.0 + 1.0/200.0 atol = 1e-5
                end
            elseif key == :f2
                for value in test_dict[key]
                    @test value ≈ - 1.0/20.0 + 1.0/200.0 atol = 1e-5
                end
            end
        end
    end
    # With lower precision
    @testset "1d_32_nogravity_novelocity" begin
        N = 10
        # Generate test moments and a test distribution
        testmoments = JuSwalbe.macroquant32_1d(ones(Float32,N),zeros(Float32,N),zeros(Float32,N))
        equilibria = calc_equilibrium_distribution(testmoments)
        # Push them in a dictonary to easily loop through the keys (fields)
        test_dict = struct2dict(equilibria)
        # Test types and lengths
        @test isa(equilibria, JuSwalbe.dist32_1d)
        @test length(keys(test_dict)) == 3
        # Test the 
        for key in keys(test_dict)
            @test isa(test_dict[key], Array{Float32,1})
            @test length(test_dict[key]) == N
            # Check if the values are actually consistent
            if key == :f0
                for value in test_dict[key]
                    @test value == 1.0
                end
            else
                for value in test_dict[key]
                    @test value == 0.0
                end
            end
        end
    end
    @testset "1d_32_gravity_novelocity" begin
        N = 10
        # Generate test moments and a test distribution
        testmoments = JuSwalbe.macroquant32_1d(ones(Float32,N),zeros(Float32,N),zeros(Float32,N))
        equilibria = calc_equilibrium_distribution(testmoments, 0.1)
        # Push them in a dictonary to easily loop through the keys (fields)
        test_dict = struct2dict(equilibria)
        # Test types and lengths
        @test isa(equilibria, JuSwalbe.dist32_1d)
        @test length(keys(test_dict)) == 3
        # Test the 
        for key in keys(test_dict)
            @test isa(test_dict[key], Array{Float32,1})
            @test length(test_dict[key]) == N
            # Check if the values are actually consistent
            if key == :f0
                for value in test_dict[key]
                    @test value == 0.95f0
                end
            else
                for value in test_dict[key]
                    @test value == 0.025f0
                end
            end
        end
    end
    @testset "1d_32_gravity_velocity" begin
        N = 10
        # Generate test moments and a test distribution
        testvelocity = Vector{Float32}(undef, N)
        testmoments = JuSwalbe.macroquant32_1d(ones(Float32,N),fill!(testvelocity, 0.1),zeros(Float32,N))
        equilibria = calc_equilibrium_distribution(testmoments, 0.1)
        # Push them in a dictonary to easily loop through the keys (fields)
        test_dict = struct2dict(equilibria)
        # Test types and lengths
        @test isa(equilibria, JuSwalbe.dist32_1d)
        @test length(keys(test_dict)) == 3
        # Test the 
        for key in keys(test_dict)
            @test isa(test_dict[key], Array{Float32,1})
            @test length(test_dict[key]) == N
            # Check if the values are actually consistent
            if key == :f0
                for value in test_dict[key]
                    @test value == 1 - 0.05f0 - 0.01f0
                end
            elseif key == :f1
                for value in test_dict[key]
                    @test value == 0.08f0
                end
            elseif key == :f2
                for value in test_dict[key]
                    @test value == - 0.02f0
                end
            end
        end
    end
    @testset "1d_32_nogravity_velocity" begin
        N = 10
        # Generate test moments and a test distribution
        testvelocity = Vector{Float32}(undef, N)
        testmoments = JuSwalbe.macroquant32_1d(ones(Float32,N),fill!(testvelocity, 0.1),zeros(Float32,N))
        equilibria = calc_equilibrium_distribution(testmoments)
        # Push them in a dictonary to easily loop through the keys (fields)
        test_dict = struct2dict(equilibria)
        # Test types and lengths
        @test isa(equilibria, JuSwalbe.dist32_1d)
        @test length(keys(test_dict)) == 3
        # Test the 
        for key in keys(test_dict)
            @test isa(test_dict[key], Array{Float32,1})
            @test length(test_dict[key]) == N
            # Check if the values are actually consistent
            if key == :f0
                for value in test_dict[key]
                    @test value == 1 - 0.01f0
                end
            elseif key == :f1
                for value in test_dict[key]
                    @test value ≈ 1.0/20.0 + 1.0/200.0 atol = 1e-5
                end
            elseif key == :f2
                for value in test_dict[key]
                    @test value ≈ - 1.0/20.0 + 1.0/200.0 atol = 1e-5
                end
            end
        end
    end

    @testset "2d_64_nogravity_novelocity" begin
        N = 3
        M = 2
        # Generate test moments and a test distribution
        testvelocity_x = Array{Float64,2}(undef, (N,M))
        testvelocity_y = Array{Float64,2}(undef, (N,M))
        testmoments = JuSwalbe.macroquant64_2d(ones(Float64,(N,M)),JuSwalbe.velocity64_2d(fill!(testvelocity_x, 0.0),fill!(testvelocity_y, 0.0)),zeros(Float64,(N,M)))
        equilibria = calc_equilibrium_distribution(testmoments)
        # Push them in a dictonary to easily loop through the keys (fields)
        test_dict = struct2dict(equilibria)
        # Test types and lengths
        @test isa(equilibria, JuSwalbe.dist64_2d)
        @test length(keys(test_dict)) == 9
        # Test the 
        for key in keys(test_dict)
            @test isa(test_dict[key], Array{Float64,2})
            @test size(test_dict[key]) == (N,M)
            @test size(test_dict[key],1) == N
            @test size(test_dict[key],2) == M
            # Check if the values are actually consistent
            if key == :f0
                for value in test_dict[key]
                    @test value == 1
                end
            else
                for value in test_dict[key]
                    @test value == 0
                end
            end
        end
    end

    @testset "2d_64_gravity_novelocity" begin
        N = 3
        M = 2
        weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
        # I use dicts as they are easier to iterate over
        weightsdict = Dict()
        for (dist, weight) in enumerate(weights)
            weightsdict[Symbol("f$(dist-1)")] = weight
        end
        # Generate test moments and a test distribution
        testvelocity_x = Array{Float64,2}(undef, (N,M))
        testvelocity_y = Array{Float64,2}(undef, (N,M))
        testmoments = JuSwalbe.macroquant64_2d(ones(Float64,(N,M)),JuSwalbe.velocity64_2d(fill!(testvelocity_x, 0.0),fill!(testvelocity_y, 0.0)),zeros(Float64,(N,M)))
        equilibria = calc_equilibrium_distribution(testmoments, 0.1)
        # Push them in a dictonary to easily loop through the keys (fields)
        test_dict = struct2dict(equilibria)
        # Test types and lengths
        @test isa(equilibria, JuSwalbe.dist64_2d)
        @test length(keys(test_dict)) == 9
        # Test the 
        for key in keys(test_dict)
            @test isa(test_dict[key], Array{Float64,2})
            @test size(test_dict[key]) == (N,M)
            @test size(test_dict[key],1) == N
            @test size(test_dict[key],2) == M
            # Check if the values are actually consistent
            if key == :f0
                for value in test_dict[key]
                    @test value == 1 - 1/3
                end
            else
                for value in test_dict[key]
                    @test value == 3/2*0.1 * weightsdict[key]
                end
            end
        end
    end

    @testset "2d_64_gravity_velocity" begin
        N = 3
        M = 2
        weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
        lattice_vel = [0 0; 1 0; 0 1; -1 0; 0 -1; 1 1; -1 1; -1 -1; 1 -1]
        # I use dicts as they are easier to iterate over
        weightsdict = Dict()
        latticeveldict = Dict()
        for (dist, weight) in enumerate(weights)
            weightsdict[Symbol("f$(dist-1)")] = weight
            latticeveldict[Symbol("f$(dist-1)")]=lattice_vel[dist,:]
        end
        # Generate test moments and a test distribution
        testvelocity_x = Array{Float64,2}(undef, (N,M))
        testvelocity_y = Array{Float64,2}(undef, (N,M))
        testmoments = JuSwalbe.macroquant64_2d(ones(Float64,(N,M)),JuSwalbe.velocity64_2d(fill!(testvelocity_x, 0.1),fill!(testvelocity_y, 0.2)),zeros(Float64,(N,M)))
        equilibria = calc_equilibrium_distribution(testmoments, 0.1)
        # Push them in a dictonary to easily loop through the keys (fields)
        test_dict = struct2dict(equilibria)
        # Test types and lengths
        @test isa(equilibria, JuSwalbe.dist64_2d)
        @test length(keys(test_dict)) == 9
        # Test the 
        for key in keys(test_dict)
            @test isa(test_dict[key], Array{Float64,2})
            @test size(test_dict[key]) == (N,M)
            @test size(test_dict[key],1) == N
            @test size(test_dict[key],2) == M
            # Check if the values are actually consistent
            if key == :f0
                for value in test_dict[key]
                    @test value == 1 - 3/10
                end
            else
                for value in test_dict[key]
                    @test value ≈ weightsdict[key] * (3.0/2.0 * 0.1 + 3 * (latticeveldict[key][1]*0.1 + latticeveldict[key][2]*0.2) 
                                                      + 9.0/2.0 * (latticeveldict[key][1]*0.1 + latticeveldict[key][2]*0.2)^2 
                                                      - 3.0/2.0 * 1.0/20.0) atol = 1e-6
                end
            end
        end
    end

    @testset "2d_64_nogravity_velocity" begin
        N = 3
        M = 2
        weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
        lattice_vel = [0 0; 1 0; 0 1; -1 0; 0 -1; 1 1; -1 1; -1 -1; 1 -1]
        # I use dicts as they are easier to iterate over
        weightsdict = Dict()
        latticeveldict = Dict()
        for (dist, weight) in enumerate(weights)
            weightsdict[Symbol("f$(dist-1)")] = weight
            latticeveldict[Symbol("f$(dist-1)")]=lattice_vel[dist,:]
        end
        # Generate test moments and a test distribution
        testvelocity_x = Array{Float64,2}(undef, (N,M))
        testvelocity_y = Array{Float64,2}(undef, (N,M))
        testmoments = JuSwalbe.macroquant64_2d(ones(Float64,(N,M)),JuSwalbe.velocity64_2d(fill!(testvelocity_x, 0.1),fill!(testvelocity_y, 0.2)),zeros(Float64,(N,M)))
        equilibria = calc_equilibrium_distribution(testmoments)
        # Push them in a dictonary to easily loop through the keys (fields)
        test_dict = struct2dict(equilibria)
        # Test types and lengths
        @test isa(equilibria, JuSwalbe.dist64_2d)
        @test length(keys(test_dict)) == 9
        # Test the 
        for key in keys(test_dict)
            @test isa(test_dict[key], Array{Float64,2})
            @test size(test_dict[key]) == (N,M)
            @test size(test_dict[key],1) == N
            @test size(test_dict[key],2) == M
            # Check if the values are actually consistent
            if key == :f0
                for value in test_dict[key]
                    @test value == 1 + 1.0/30.0
                end
            else
                for value in test_dict[key]
                    @test value ≈ weightsdict[key] * (3 * (latticeveldict[key][1]*0.1 + latticeveldict[key][2]*0.2) 
                                                      + 9.0/2.0 * (latticeveldict[key][1]*0.1 + latticeveldict[key][2]*0.2)^2 
                                                      - 3.0/2.0 * 1.0/20.0) atol = 1e-6
                end
            end
        end
    end
   
end