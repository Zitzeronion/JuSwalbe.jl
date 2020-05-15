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
end