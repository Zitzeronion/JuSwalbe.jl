@testset "Equilibria" begin
tolerances = Dict(Float16 => 1e-3, Float32 => 1e-5, Float64 => 1e-7, Real => 1e-7)
    @testset "One dimension -velocity -gravity" begin
        N = 10
        # Generate test moments and a test distribution
        typelist = [Float16 Float32 Float64 Real]
        for type in typelist
            testmoments = JuSwalbe.Macroquant{Vector{type}, Vector{type}}(ones(type,N),zeros(type,N),zeros(type,N),zeros(type,N))
            equilibria = calc_equilibrium_distribution(testmoments)
            # Push them in a dictonary to easily loop through the keys (fields)
            test_dict = struct2dict(equilibria)
            # Test types and lengths
            @test isa(equilibria, JuSwalbe.DistributionD1Q3)
            @test length(keys(test_dict)) == 3
            # Test the 
            for key in keys(test_dict)
                @test isa(test_dict[key], Array{type,1})
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
    end
    @testset "One dimension -velocity +gravity" begin
        N = 10
        # Generate test moments and a test distribution
        typelist = [Float16 Float32 Float64 Real]
        for type in typelist
            testmoments = JuSwalbe.Macroquant{Vector{type}, Vector{type}}(ones(type,N),zeros(type,N),zeros(type,N),zeros(type,N))
            equilibria = calc_equilibrium_distribution(testmoments; gravity=type(0.1))
            # Push them in a dictonary to easily loop through the keys (fields)
            test_dict = struct2dict(equilibria)
            # Test types and lengths
            @test isa(equilibria, JuSwalbe.DistributionD1Q3)
            @test length(keys(test_dict)) == 3
            # Test the 
            for key in keys(test_dict)
                @test isa(test_dict[key], Array{type,1})
                @test length(test_dict[key]) == N
                # Check if the values are actually consistent
                if key == :f0
                    for value in test_dict[key]
                        @test value == type(0.95)
                    end
                else
                    for value in test_dict[key]
                        @test value == type(1.0/40.0)
                    end
                end
            end
        end
    end
    @testset "One dimension +velocity +gravity" begin
        N = 10
        # Generate test moments and a test distribution
        typelist = [Float16 Float32 Float64 Real]
        for type in typelist
            test_velocity = ones(type, N) * type(0.1)
            testmoments = JuSwalbe.Macroquant{Vector{type}, Vector{type}}(ones(type,N),test_velocity,zeros(type,N),zeros(type,N))
            equilibria = calc_equilibrium_distribution(testmoments; gravity=type(0.1))
            # Push them in a dictonary to easily loop through the keys (fields)
            test_dict = struct2dict(equilibria)
            # Test types and lengths
            @test isa(equilibria, JuSwalbe.DistributionD1Q3)
            @test length(keys(test_dict)) == 3
            # Test the 
            for key in keys(test_dict)
                @test isa(test_dict[key], Array{type,1})
                @test length(test_dict[key]) == N
                # Check if the values are actually consistent
                if key == :f0
                    for value in test_dict[key]
                        @test value - type(1 - 1.0/20.0 - 0.01) ≈ type(0.0) atol=tolerances[type]
                    end
                elseif key == :f1
                    for value in test_dict[key]
                        @test value - type(1.0/40.0 + 1.0/20.0 + 1.0/200.0) ≈ type(0.0) atol=tolerances[type]
                    end
                elseif key == :f2
                    for value in test_dict[key]
                        @test value - type(1.0/40.0 - 1.0/20.0 + 1.0/200.0) ≈ type(0.0) atol=tolerances[type]
                    end
                end
            end
        end
    end
    @testset "One dimension +velocity -gravity" begin
        N = 10
        # Generate test moments and a test distribution
        typelist = [Float16 Float32 Float64 Real]
        for type in typelist
            test_velocity = ones(type, N) * type(0.1)
            testmoments = JuSwalbe.Macroquant{Vector{type}, Vector{type}}(ones(type,N),test_velocity,zeros(type,N),zeros(type,N))
            equilibria = calc_equilibrium_distribution(testmoments)
            # Push them in a dictonary to easily loop through the keys (fields)
            test_dict = struct2dict(equilibria)
            # Test types and lengths
            @test isa(equilibria, JuSwalbe.DistributionD1Q3)
            @test length(keys(test_dict)) == 3
            # Test the 
            for key in keys(test_dict)
                @test isa(test_dict[key], Array{type,1})
                @test length(test_dict[key]) == N
                # Check if the values are actually consistent
                if key == :f0
                    for value in test_dict[key]
                        @test value ≈ type(1 - 0.01) atol = 1e-3
                    end
                elseif key == :f1
                    for value in test_dict[key]
                        @test value ≈ type(1.0/20.0 + 1.0/200.0) atol = 1e-3
                    end
                elseif key == :f2
                    for value in test_dict[key]
                        @test value ≈ type(- 1.0/20.0 + 1.0/200.0) atol = 1e-3
                    end
                end
            end
        end
    end
    

    @testset "Two dimension -velocity -gravity" begin
        N = 3
        M = 2
        # Generate test moments and a test distribution
        typelist = [Float16 Float32 Float64]
        for type in typelist
            testvelocity_x = zeros(type, (N,M))
            testvelocity_y = zeros(type, (N,M))
            test_vel = JuSwalbe.Twovector(x=testvelocity_x, y=testvelocity_y) 
            testmoments = JuSwalbe.Macroquant(height=ones(type,(N,M)), velocity=test_vel, pressure=zeros(type,(N,M)), energy=zeros(type,(N,M)))
            equilibria = calc_equilibrium_distribution(testmoments, gravity=type(0.0))
            # Push them in a dictonary to easily loop through the keys (fields)
            test_dict = struct2dict(equilibria)
            # Test types and lengths
            @test isa(equilibria, JuSwalbe.DistributionD2Q9)
            @test length(keys(test_dict)) == 9
            # Test the 
            for key in keys(test_dict)
                @test isa(test_dict[key], Array{type,2})
                @test size(test_dict[key]) == (N,M)
                @test size(test_dict[key],1) == N
                @test size(test_dict[key],2) == M
                # Check if the values are actually consistent
                if key == :f0
                    for value in test_dict[key]
                        @test value == type(1)
                    end
                else
                    for value in test_dict[key]
                        @test value == type(0)
                    end
                end
            end
        end
    end

    @testset "Two dimension -velocity +gravity" begin
        N = 3
        M = 2
        weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
        typelist = [Float16 Float32 Float64]
        for type in typelist
            # I use dicts as they are easier to iterate over
            weightsdict = Dict()
            # For type correction
            weightcorrected = Array{type,1}(undef, 9)
            for (dist, weight) in enumerate(weights)
                weightcorrected[dist] = weights[dist]
                weightsdict[Symbol("f$(dist-1)")] = weightcorrected[dist]
            end
            # Generate test moments and a test distribution
            testvelocity_x = zeros(type, (N,M))
            testvelocity_y = zeros(type, (N,M))
            test_vel = JuSwalbe.Twovector(testvelocity_x, testvelocity_y)
            testmoments = JuSwalbe.Macroquant{Matrix{type}, JuSwalbe.Twovector{Matrix{type}}}(ones(type,(N,M)), test_vel, zeros(type,(N,M)), zeros(type,(N,M)))
            equilibria = calc_equilibrium_distribution(testmoments; gravity=type(0.1))
            # Push them in a dictonary to easily loop through the keys (fields)
            test_dict = struct2dict(equilibria)
            # Test types and lengths
            @test isa(equilibria, JuSwalbe.DistributionD2Q9)
            @test length(keys(test_dict)) == 9
            # Test the 
            for key in keys(test_dict)
                @test isa(test_dict[key], Array{type,2})
                @test size(test_dict[key]) == (N,M)
                @test size(test_dict[key],1) == N
                @test size(test_dict[key],2) == M
                # Check if the values are actually consistent
                if key == :f0
                    for value in test_dict[key]
                        @test value - type(1 - 1/12) ≈ type(0.0) atol = tolerances[type]
                    end
                else
                    for value in test_dict[key]
                        @test value ≈ type(3/2*0.1 * weightsdict[key]) atol = 1e-3
                    end
                end
            end
        end
    end

    @testset "Two dimension +velocity +gravity" begin
        N = 3
        M = 2
        weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
        lattice_vel = [0 0; 1 0; 0 1; -1 0; 0 -1; 1 1; -1 1; -1 -1; 1 -1]
        # I use dicts as they are easier to iterate over
        weightsdict = Dict()
        latticeveldict = Dict()
        typelist = [Float16 Float32 Float64]
        for type in typelist
            # For type correction
            weightcorrected = Array{type,1}(undef, 9)
            for (dist, weight) in enumerate(weights)
                weightcorrected[dist] = weights[dist]
                weightsdict[Symbol("f$(dist-1)")] = weightcorrected[dist]
                latticeveldict[Symbol("f$(dist-1)")]=lattice_vel[dist,:]
            end
            # Generate test moments and a test distribution
            testvelocity_x = ones(type, (N,M)) * type(0.1)
            testvelocity_y = ones(type, (N,M)) * type(0.2)
            test_vel = JuSwalbe.Twovector(testvelocity_x, testvelocity_y)
            testmoments = JuSwalbe.Macroquant{Matrix{type}, JuSwalbe.Twovector{Matrix{type}}}(ones(type,(N,M)), test_vel, zeros(type,(N,M)), zeros(type,(N,M)))
            equilibria = calc_equilibrium_distribution(testmoments; gravity=type(0.1))
            # Push them in a dictonary to easily loop through the keys (fields)
            test_dict = struct2dict(equilibria)
            # Test types and lengths
            @test isa(equilibria, JuSwalbe.DistributionD2Q9)
            @test length(keys(test_dict)) == 9
            # Test the 
            for key in keys(test_dict)
                @test isa(test_dict[key], Array{type,2})
                @test size(test_dict[key]) == (N,M)
                @test size(test_dict[key],1) == N
                @test size(test_dict[key],2) == M
                # Check if the values are actually consistent
                if key == :f0
                    for value in test_dict[key]
                        @test value - type(1 - 7/60) ≈ type(0) atol=tolerances[type]
                    end
                else
                    for value in test_dict[key]
                        @test value ≈ type(weightsdict[key] * (3.0/2.0 * 0.1 + 3 * (latticeveldict[key][1]*0.1 + latticeveldict[key][2]*0.2) 
                                                        + 9.0/2.0 * (latticeveldict[key][1]*0.1 + latticeveldict[key][2]*0.2)^2 
                                                        - 3.0/2.0 * 1.0/20.0)) atol = 1e-3
                    end
                end
            end
        end
    end

    @testset "Two dimension +velocity -gravity" begin
        N = 3
        M = 2
        weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
        lattice_vel = [0 0; 1 0; 0 1; -1 0; 0 -1; 1 1; -1 1; -1 -1; 1 -1]
        # I use dicts as they are easier to iterate over
        weightsdict = Dict()
        latticeveldict = Dict()
        typelist = [Float16 Float32 Float64]
        for type in typelist
            # For type correction
            weightcorrected = Array{type,1}(undef, 9)
            for (dist, weight) in enumerate(weights)
                weightcorrected[dist] = weights[dist]
                weightsdict[Symbol("f$(dist-1)")] = weightcorrected[dist]
                latticeveldict[Symbol("f$(dist-1)")]=lattice_vel[dist,:]
            end
            # Generate test moments and a test distribution
            testvelocity_x = ones(type, (N,M)) * type(0.1)
            testvelocity_y = ones(type, (N,M)) * type(0.2)
            test_vel = JuSwalbe.Twovector(testvelocity_x, testvelocity_y)
            testmoments = JuSwalbe.Macroquant{Matrix{type}, JuSwalbe.Twovector{Matrix{type}}}(ones(type,(N,M)), test_vel, zeros(type,(N,M)), zeros(type,(N,M)))
            equilibria = calc_equilibrium_distribution(testmoments, gravity=type(0.0))
            # Push them in a dictonary to easily loop through the keys (fields)
            test_dict = struct2dict(equilibria)
            # Test types and lengths
            @test isa(equilibria, JuSwalbe.DistributionD2Q9)
            @test length(keys(test_dict)) == 9
            # Test the 
            for key in keys(test_dict)
                @test isa(test_dict[key], Array{type,2})
                @test size(test_dict[key]) == (N,M)
                @test size(test_dict[key],1) == N
                @test size(test_dict[key],2) == M
                # Check if the values are actually consistent
                if key == :f0
                    for value in test_dict[key]
                        @test value - type(1 - 1.0/30.0) ≈ type(0.0) atol=tolerances[type]
                    end
                else
                    for value in test_dict[key]
                        @test value ≈ type(weightsdict[key] * (3 * (latticeveldict[key][1]*0.1 + latticeveldict[key][2]*0.2) 
                                                        + 9.0/2.0 * (latticeveldict[key][1]*0.1 + latticeveldict[key][2]*0.2)^2 
                                                        - 3.0/2.0 * 1.0/20.0)) atol = 1e-3
                    end
                end
            end
        end
    end
end