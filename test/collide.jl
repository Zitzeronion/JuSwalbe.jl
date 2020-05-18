@testset "Collision operator" begin
    tolerances = Dict(Float16 => 1e-4, Float32 => 1e-5, Float64 => 1e-7, Real => 1e-7)
    @testset "One spatial dimension" begin
        @testset "Collision -forcing -temp" begin
            N = 8
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testmom = JuSwalbe.Macroquant(height = ones(type, N), velocity=zeros(type, N), pressure=zeros(type, N), energy=zeros(type, N))
                testforce = JuSwalbe.Forces(slip = zeros(type, N), h∇p=zeros(type, N), bathymetry=zeros(type, N), thermal=zeros(type, N))
                testtempdist = JuSwalbe.DistributionD1Q3(f0 = zeros(type, N), f1 = zeros(type, N), f2 = zeros(type, N))
                input = JuSwalbe.Inputconstants(τ = type(1.0), gravity=type(0.0))
                @test isa(testmom, JuSwalbe.Macroquant{Vector{type}, Vector{type}})
                @test isa(testtempdist, JuSwalbe.DistributionD1Q3{Vector{type}})
                @test isa(input, JuSwalbe.Inputconstants)
                testequi = calc_equilibrium_distribution(testmom, gravity=type(input.gravity))
                testdist = collisionBGK(testmom, testforce, testtempdist, input)
                @test isa(testdist, JuSwalbe.DistributionD1Q3{Vector{type}})
                for i in 1:N
                    @test testdist.f0[i] -  testequi.f0[i] ≈ type(0.0) atol=tolerances[type]
                    @test testdist.f1[i] -  testequi.f1[i] ≈ type(0.0) atol=tolerances[type]
                    @test testdist.f2[i] -  testequi.f2[i] ≈ type(0.0) atol=tolerances[type]
                end
            end
        end
        @testset "Collision +forcing -temp" begin
            N = 8
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testmom = JuSwalbe.Macroquant(height = ones(type, N), velocity=zeros(type, N), pressure=zeros(type, N), energy=zeros(type, N))
                testforce = JuSwalbe.Forces(slip = ones(type, N), h∇p=ones(type, N), bathymetry=fill(type(0.1), N), thermal=zeros(type, N))
                testtempdist = JuSwalbe.DistributionD1Q3(f0 = zeros(type, N), f1 = zeros(type, N), f2 = zeros(type, N))
                input = JuSwalbe.Inputconstants(τ = type(1.0), gravity=type(0.1))
                @test isa(testmom, JuSwalbe.Macroquant{Vector{type}, Vector{type}})
                @test isa(testtempdist, JuSwalbe.DistributionD1Q3{Vector{type}})
                @test isa(input, JuSwalbe.Inputconstants)
                testequi = calc_equilibrium_distribution(testmom, gravity=type(input.gravity))
                testdist = collisionBGK(testmom, testforce, testtempdist, input)
                @test isa(testdist, JuSwalbe.DistributionD1Q3{Vector{type}})
                for i in 1:N
                    @test testdist.f0[i] -  testequi.f0[i] ≈ type(0.0) atol=tolerances[type]
                    @test (testdist.f1[i] -  testequi.f1[i] - type(3/6) * type(2.1)) ≈ type(0.0) atol=tolerances[type]
                    @test (testdist.f2[i] -  testequi.f2[i] + type(3/6) * type(2.1)) ≈ type(0.0) atol=tolerances[type]
                end
            end
        end
        @testset "Collision +forcing +temp" begin
            N = 8
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testmom = JuSwalbe.Macroquant(height = ones(type, N), velocity=zeros(type, N), pressure=zeros(type, N), energy=zeros(type, N))
                testforce = JuSwalbe.Forces(slip = ones(type, N), h∇p=ones(type, N), bathymetry=fill(type(0.1), N), thermal=zeros(type, N))
                testtempdist = JuSwalbe.DistributionD1Q3(f0 = fill(type(0.1), N), f1 = fill(type(0.2), N), f2 = fill(type(-0.2), N))
                input = JuSwalbe.Inputconstants(τ = type(0.75), gravity=type(0.1))
                @test isa(testmom, JuSwalbe.Macroquant{Vector{type}, Vector{type}})
                @test isa(testtempdist, JuSwalbe.DistributionD1Q3{Vector{type}})
                @test isa(input, JuSwalbe.Inputconstants)
                testequi = calc_equilibrium_distribution(testmom, gravity=type(input.gravity))
                testdist = collisionBGK(testmom, testforce, testtempdist, input)
                @test isa(testdist, JuSwalbe.DistributionD1Q3{Vector{type}})
                for i in 1:N
                    @test testdist.f0[i] - (type(1-1/input.τ)*testtempdist.f0[i] + 1/input.τ * testequi.f0[i]) ≈ type(0.0) atol=tolerances[type]
                    @test (testdist.f1[i] -  (type(1-1/input.τ)*testtempdist.f1[i] + 1/input.τ * testequi.f1[i] + type(3/6) * type(2.1))) ≈ type(0.0) atol=tolerances[type]
                    @test (testdist.f2[i] -  (type(1-1/input.τ)*testtempdist.f2[i] + 1/input.τ * testequi.f2[i] - type(3/6) * type(2.1))) ≈ type(0.0) atol=tolerances[type]
                end
            end
        end
        @testset "Collision -forcing +temp" begin
            N = 8
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testmom = JuSwalbe.Macroquant(height = ones(type, N), velocity=zeros(type, N), pressure=zeros(type, N), energy=zeros(type, N))
                testforce = JuSwalbe.Forces(slip=zeros(type, N), h∇p=zeros(type, N), bathymetry=zeros(type, N), thermal=zeros(type, N))
                testtempdist = JuSwalbe.DistributionD1Q3(f0 = fill(type(0.1), N), f1 = fill(type(0.2), N), f2 = fill(type(-0.2), N))
                input = JuSwalbe.Inputconstants(τ = type(0.75), gravity=type(0.0))
                @test isa(testmom, JuSwalbe.Macroquant{Vector{type}, Vector{type}})
                @test isa(testtempdist, JuSwalbe.DistributionD1Q3{Vector{type}})
                @test isa(input, JuSwalbe.Inputconstants)
                testequi = calc_equilibrium_distribution(testmom, gravity=type(input.gravity))
                testdist = collisionBGK(testmom, testforce, testtempdist, input)
                @test isa(testdist, JuSwalbe.DistributionD1Q3{Vector{type}})
                for i in 1:N
                    @test testdist.f0[i] - (type(1-1/input.τ)*testtempdist.f0[i] + 1/input.τ * testequi.f0[i]) ≈ type(0.0) atol=tolerances[type]
                    @test (testdist.f1[i] -  (type(1-1/input.τ)*testtempdist.f1[i] + 1/input.τ * testequi.f1[i])) ≈ type(0.0) atol=tolerances[type]
                    @test (testdist.f2[i] -  (type(1-1/input.τ)*testtempdist.f2[i] + 1/input.τ * testequi.f2[i])) ≈ type(0.0) atol=tolerances[type]
                end
            end
        end
    end
    @testset "Two spatial dimensions" begin
        @testset "Collision -forcing -temp" begin
            N, M = 8, 6
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testvel = JuSwalbe.Twovector(x=fill(type(0.1),(N,M)), y=fill(type(0.1),(N,M)))
                zerovec = JuSwalbe.Twovector(x=zeros(type, (N,M)), y=zeros(type, (N,M)))
                testmom = JuSwalbe.Macroquant(height = ones(type, (N,M)), velocity=zerovec, pressure=zeros(type, (N,M)), energy=zeros(type, (N,M)))
                testforce = JuSwalbe.Forces(slip = zerovec, h∇p=zerovec, bathymetry=zerovec, thermal=zerovec)
                testtempdist = JuSwalbe.DistributionD2Q9(f0 = zeros(type, (N,M)), 
                                                         f1 = zeros(type, (N,M)), 
                                                         f2 = zeros(type, (N,M)),
                                                         f3 = zeros(type, (N,M)),
                                                         f4 = zeros(type, (N,M)),
                                                         f5 = zeros(type, (N,M)),
                                                         f6 = zeros(type, (N,M)),
                                                         f7 = zeros(type, (N,M)),
                                                         f8 = zeros(type, (N,M)))

                input = JuSwalbe.Inputconstants(τ = type(1.0), gravity=type(0.0))
                @test isa(testmom, JuSwalbe.Macroquant{Matrix{type}, JuSwalbe.Twovector{Matrix{type}}})
                @test isa(testtempdist, JuSwalbe.DistributionD2Q9{Matrix{type}})
                @test isa(input, JuSwalbe.Inputconstants)
                testequi = calc_equilibrium_distribution(testmom, gravity=type(input.gravity))
                testdist = collisionBGK(testmom, testforce, testtempdist, input)
                @test isa(testdist, JuSwalbe.DistributionD2Q9{Matrix{type}})
                testdict = struct2dict(testdist)
                equidict = struct2dict(testequi)
                for key in keys(testdict)
                    for i in 1:N, j in 1:M
                        @test testdict[key][i,j] -  equidict[key][i,j] ≈ type(0.0) atol=tolerances[type]
                    end
                end
            end
        end
        @testset "Collision +forcing -temp" begin
            N, M = 5, 5
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testvel = JuSwalbe.Twovector(x=fill(type(0.1),(N,M)), y=fill(type(0.1),(N,M)))
                zerovec = JuSwalbe.Twovector(x=zeros(type, (N,M)), y=zeros(type, (N,M)))
                testmom = JuSwalbe.Macroquant(height = ones(type, (N,M)), velocity=zerovec, pressure=zeros(type, (N,M)), energy=zeros(type, (N,M)))
                testforce = JuSwalbe.Forces(slip = testvel, h∇p=zerovec, bathymetry=zerovec, thermal=zerovec)
                testtempdist = JuSwalbe.DistributionD2Q9(f0 = zeros(type, (N,M)), f1 = zeros(type, (N,M)), f2 = zeros(type, (N,M)),
                                                         f3 = zeros(type, (N,M)), f4 = zeros(type, (N,M)), f5 = zeros(type, (N,M)),
                                                         f6 = zeros(type, (N,M)), f7 = zeros(type, (N,M)), f8 = zeros(type, (N,M)))

                input = JuSwalbe.Inputconstants(τ = type(1.0), gravity=type(0.0))
                @test isa(testmom, JuSwalbe.Macroquant{Matrix{type}, JuSwalbe.Twovector{Matrix{type}}})
                @test isa(testtempdist, JuSwalbe.DistributionD2Q9{Matrix{type}})
                @test isa(testforce, JuSwalbe.Forces{JuSwalbe.Twovector{Matrix{type}}})
                @test isa(input, JuSwalbe.Inputconstants)
                testequi = calc_equilibrium_distribution(testmom, gravity=type(input.gravity))
                testdist = collisionBGK(testmom, testforce, testtempdist, input)
                @test isa(testdist, JuSwalbe.DistributionD2Q9{Matrix{type}})
                testdict = struct2dict(testdist)
                equidict = struct2dict(testequi)
                # Solution is just the force: 3w_i/c_i^2 * c_i * F
                solutiondict = Dict(:f0 => type(0), :f1 => type(1/30), :f2 => type(1/30), :f3 => type(-1/30), :f4 => type(-1/30),
                                    :f5 => type(1/120), :f6 => type(0), :f7 => type(-1/120), :f8 => type(0))
                for key in keys(testdict)
                    for i in 1:N, j in 1:M
                        @test testdict[key][i,j] -  (equidict[key][i,j] + solutiondict[key])  ≈ type(0.0) atol=tolerances[type]
                    end
                end
            end
        end
        @testset "Collision +forcing +temp" begin
            N, M = 5, 5
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testvel = JuSwalbe.Twovector(x=fill(type(0.1),(N,M)), y=fill(type(0.1),(N,M)))
                zerovec = JuSwalbe.Twovector(x=zeros(type, (N,M)), y=zeros(type, (N,M)))
                testmom = JuSwalbe.Macroquant(height = ones(type, (N,M)), velocity=zerovec, pressure=zeros(type, (N,M)), energy=zeros(type, (N,M)))
                testforce = JuSwalbe.Forces(slip = testvel, h∇p=zerovec, bathymetry=zerovec, thermal=zerovec)
                testtempdist = JuSwalbe.DistributionD2Q9(f0 = fill(type(0.1), (N,M)), f1 = fill(type(0.2), (N,M)), f2 = fill(type(-0.2), (N,M)),
                                                         f3 = zeros(type, (N,M)), f4 = zeros(type, (N,M)), f5 = zeros(type, (N,M)),
                                                         f6 = zeros(type, (N,M)), f7 = zeros(type, (N,M)), f8 = zeros(type, (N,M)))

                input = JuSwalbe.Inputconstants(τ = type(0.75), gravity=type(0.0))
                @test isa(testmom, JuSwalbe.Macroquant{Matrix{type}, JuSwalbe.Twovector{Matrix{type}}})
                @test isa(testtempdist, JuSwalbe.DistributionD2Q9{Matrix{type}})
                @test isa(testforce, JuSwalbe.Forces{JuSwalbe.Twovector{Matrix{type}}})
                @test isa(input, JuSwalbe.Inputconstants)
                testequi = calc_equilibrium_distribution(testmom, gravity=type(input.gravity))
                testdist = collisionBGK(testmom, testforce, testtempdist, input)
                @test isa(testdist, JuSwalbe.DistributionD2Q9{Matrix{type}})
                testdict = struct2dict(testdist)
                equidict = struct2dict(testequi)
                tempdict = struct2dict(testtempdist)
                soldict = Dict(:f0 => type(0), :f1 => type(1/30), :f2 => type(1/30), :f3 => type(-1/30), :f4 => type(-1/30),
                               :f5 => type(1/120), :f6 => type(0), :f7 => type(-1/120), :f8 => type(0))
                for key in keys(testdict)
                    for i in 1:N, j in 1:M
                        @test testdict[key][i,j] - (type(1-1/input.τ)*tempdict[key][i,j] + 1/input.τ * equidict[key][i,j] + soldict[key]) ≈ type(0.0) atol=tolerances[type]
                    end
                end
            end
        end
        @testset "Collision -forcing +temp" begin
            N, M = 5, 5
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testvel = JuSwalbe.Twovector(x=fill(type(0.1),(N,M)), y=fill(type(0.1),(N,M)))
                zerovec = JuSwalbe.Twovector(x=zeros(type, (N,M)), y=zeros(type, (N,M)))
                testmom = JuSwalbe.Macroquant(height = ones(type, (N,M)), velocity=zerovec, pressure=zeros(type, (N,M)), energy=zeros(type, (N,M)))
                testforce = JuSwalbe.Forces(slip = zerovec, h∇p=zerovec, bathymetry=zerovec, thermal=zerovec)
                testtempdist = JuSwalbe.DistributionD2Q9(f0 = fill(type(0.1), (N,M)), f1 = fill(type(0.2), (N,M)), f2 = fill(type(-0.2), (N,M)),
                                                        f3 = zeros(type, (N,M)), f4 = zeros(type, (N,M)), f5 = zeros(type, (N,M)),
                                                        f6 = zeros(type, (N,M)), f7 = zeros(type, (N,M)), f8 = zeros(type, (N,M)))

                input = JuSwalbe.Inputconstants(τ = type(0.75), gravity=type(0.0))
                @test isa(testmom, JuSwalbe.Macroquant{Matrix{type}, JuSwalbe.Twovector{Matrix{type}}})
                @test isa(testtempdist, JuSwalbe.DistributionD2Q9{Matrix{type}})
                @test isa(testforce, JuSwalbe.Forces{JuSwalbe.Twovector{Matrix{type}}})
                @test isa(input, JuSwalbe.Inputconstants)
                testequi = calc_equilibrium_distribution(testmom, gravity=type(input.gravity))
                testdist = collisionBGK(testmom, testforce, testtempdist, input)
                @test isa(testdist, JuSwalbe.DistributionD2Q9{Matrix{type}})
                testdict = struct2dict(testdist)
                equidict = struct2dict(testequi)
                tempdict = struct2dict(testtempdist)
                
                for key in keys(testdict)
                    for i in 1:N, j in 1:M
                        @test testdict[key][i,j] - (type(1-1/input.τ)*tempdict[key][i,j] + 1/input.τ * equidict[key][i,j]) ≈ type(0.0) atol=tolerances[type]
                    end
                end
            end
        end
    end
end

@testset "Streaming function" begin
    @testset "One dimension periodic" begin
        N = 8
        typelist = [Float16 Float32 Float64]
        for type in typelist
            testtempdist = JuSwalbe.DistributionD1Q3(f0 = fill(type(0.1), N), f1 = fill(type(0.2), N), f2 = fill(type(-0.2), N))
            testtempdist.f1[3] = type(-0.1)
            testtempdist.f2[3] = type(0.1)

            @test isa(testtempdist, JuSwalbe.DistributionD1Q3{Vector{type}})
            streamdistperiodic!(testtempdist)
            @test isa(testtempdist, JuSwalbe.DistributionD1Q3{Vector{type}})
            @test length(testtempdist.f0) == 8
            @test length(testtempdist.f1) == 8
            @test length(testtempdist.f2) == 8
            for i in testtempdist.f0
                @test i == type(0.1)
            end
            @test testtempdist.f1[4] == type(-0.1)
            @test testtempdist.f1[3] == type(0.2)
            @test testtempdist.f1[5] == type(0.2)
            @test testtempdist.f2[2] == type(0.1)
            @test testtempdist.f2[3] == type(-0.2)
            @test testtempdist.f2[4] == type(-0.2)
            
        end
    end
    @testset "Two dimensions periodic" begin
        N, M = 8, 6
        typelist = [Float16 Float32 Float64]
        for type in typelist
            testtempdist = JuSwalbe.DistributionD2Q9(f0 = fill(type(0.1), (N,M)), f1 = fill(type(0.2), (N,M)), f2 = fill(type(-0.2), (N,M)),
                                                     f3 = fill(type(0.3), (N,M)), f4 = fill(type(0.1), (N,M)), f5 = fill(type(-0.1), (N,M)),
                                                     f6 = fill(type(0.6), (N,M)), f7 = fill(type(0.5), (N,M)), f8 = fill(type(-0.5), (N,M)))
            testtempdist.f1[3,4] = type(1)
            testtempdist.f2[3,3] = type(1)
            testtempdist.f3[3,4] = type(1)
            testtempdist.f4[3,3] = type(1)
            testtempdist.f5[3,4] = type(1)
            testtempdist.f6[3,3] = type(1)
            testtempdist.f7[3,4] = type(1)
            testtempdist.f8[3,3] = type(1)

            @test isa(testtempdist, JuSwalbe.DistributionD2Q9{Matrix{type}})
            testdist = streamdistperiodic!(testtempdist)
            @test isa(testtempdist, JuSwalbe.DistributionD2Q9{Matrix{type}})
            @test size(testdist.f0) == (N,M)
            @test size(testdist.f1) == (N,M)
            @test size(testdist.f2) == (N,M)
            @test size(testdist.f3) == (N,M)
            @test size(testdist.f4) == (N,M)
            @test size(testdist.f5) == (N,M)
            @test size(testdist.f6) == (N,M)
            @test size(testdist.f7) == (N,M)
            @test size(testdist.f8) == (N,M)
            for i in testdist.f0
                @test i == type(0.1)
            end
            for i in 1:N, j in 1:M
                if i == 4 && j == 4
                    @test testdist.f1[i,j] == type(1)
                else
                    @test testdist.f1[i,j] == type(0.2) 
                end
            end
            for i in 1:N, j in 1:M
                if i == 3 && j == 4
                    @test testdist.f2[i,j] == type(1)
                else
                    @test testdist.f2[i,j] == type(-0.2) 
                end
            end
            for i in 1:N, j in 1:M
                if i == 2 && j == 4
                    @test testdist.f3[i,j] == type(1)
                else
                    @test testdist.f3[i,j] == type(0.3) 
                end
            end
            for i in 1:N, j in 1:M
                if i == 3 && j == 2
                    @test testdist.f4[i,j] == type(1)
                else
                    @test testdist.f4[i,j] == type(0.1) 
                end
            end
            for i in 1:N, j in 1:M
                if i == 4 && j == 5
                    @test testdist.f5[i,j] == type(1)
                else
                    @test testdist.f5[i,j] == type(-0.1) 
                end
            end
            for i in 1:N, j in 1:M
                if i == 2 && j == 4
                    @test testdist.f6[i,j] == type(1)
                else
                    @test testdist.f6[i,j] == type(0.6) 
                end
            end
            for i in 1:N, j in 1:M
                if i == 2 && j == 3
                    @test testdist.f7[i,j] == type(1)
                else
                    @test testdist.f7[i,j] == type(0.5) 
                end
            end
            for i in 1:N, j in 1:M
                if i == 4 && j == 2
                    @test testdist.f8[i,j] == type(1)
                else
                    @test testdist.f8[i,j] == type(-0.5) 
                end
            end
        end
    end
end