@testset "Forcing" begin
    tolerances = Dict(Float16 => 1e-4, Float32 => 1e-5, Float64 => 1e-7, Real => 1e-7)
    @testset "Slippage one dimension with velocity" begin
        N = 8
        typelist = [Float16 Float32 Float64]
        for type in typelist
            delta_set = [1.0 0.5 0.0]
            sol_set = [1/110 2/115 1/20]
            # Test different delta values
            for delta in 1:3
                label = delta_set[delta]
                @testset "δ = $label" begin
                    moment = JuSwalbe.Macroquant(height=ones(type, N), velocity=fill(type(0.1), N), pressure=zeros(type, N), energy=zeros(type, N))
                    forces = JuSwalbe.Forces(slip=zeros(type, N), thermal=zeros(type, N), h∇p=zeros(type, N), bathymetry=zeros(type, N))
                    input = JuSwalbe.Inputconstants(δ=type(delta_set[delta]), τ=type(1.0))
                    slippage = computeslip(moment, forces, input)
                    @test isa(forces.slip, Array{type, 1})
                    @test length(slippage) == N
                    
                    for i in 1:N
                        @test slippage[i] - type(sol_set[delta]) ≈ type(0.0) atol=tolerances[type]
                    end
                end
            end
        end
    end
    @testset "Slippage one dimension without velocity" begin
        N, M = 8, 9
        typelist = [Float16 Float32 Float64]
        for type in typelist
            # Test different delta values
            delta_set = [1.0 0.5 0.0]
            for delta in 1:3
                label = delta_set[delta]
                @testset "δ = $label" begin
                    moment = JuSwalbe.Macroquant(height=ones(type, N), velocity=zeros(type, N), pressure=zeros(type, N), energy=zeros(type, N))
                    forces = JuSwalbe.Forces(slip=zeros(type, N), thermal=zeros(type, N), h∇p=zeros(type, N), bathymetry=zeros(type, N))
                    input = JuSwalbe.Inputconstants(δ=type(delta_set[delta]), τ=type(1.0))
                    slippage = computeslip(moment, forces, input)
                    @test isa(forces.slip, Array{type, 1})
                    @test length(slippage) == N
                    for i in 1:N
                        @test slippage[i] ≈ type(0.0) atol=tolerances[type]
                    end
                end
            end
        end
    end    
    @testset "Slippage two dimensions with velocity" begin
        N, M = 8, 9
        typelist = [Float16 Float32 Float64]
        for type in typelist
            # Test different delta values
            delta_set = [1.0 0.5 0.0]
            sol_set = [1/110 2/115 1/20]
            for delta in 1:3
                label = delta_set[delta]
                @testset "δ = $label" begin
                    vel = JuSwalbe.Twovector(x=fill(type(0.1), (N,M)), y=fill(type(-0.1), (N,M)))
                    zerovec = JuSwalbe.Twovector(x=zeros(type ,(N,M)), y=zeros(type ,(N,M)))
                    moment = JuSwalbe.Macroquant(height=ones(type, (N,M)), velocity=vel, pressure=zeros(type, (N,M)), energy=zeros(type, (N,M)))
                    forces = JuSwalbe.Forces(slip=zerovec, thermal=zerovec, h∇p=zerovec, bathymetry=zerovec)
                    input = JuSwalbe.Inputconstants(δ=type(delta_set[delta]), τ=type(1.0))
                    slippage = computeslip(moment, forces, input)
                    @test isa(zerovec, JuSwalbe.Twovector{Array{type, 2}})
                    @test isa(moment, JuSwalbe.Macroquant{Array{type, 2}, JuSwalbe.Twovector{Array{type, 2}}})
                    @test isa(forces, JuSwalbe.Forces{JuSwalbe.Twovector{Array{type, 2}}})
                    @test isa(input, JuSwalbe.Inputconstants)
                    
                    for i in 1:N, j in 1:M
                        @test slippage.x[i,j] - type(sol_set[delta]) ≈ type(0.0) atol=tolerances[type]
                        @test slippage.y[i,j] + type(sol_set[delta]) ≈ type(0.0) atol=tolerances[type]
                    end
                end
            end
        end
    end
    @testset "Slippage two dimensions without velocity" begin
        N, M = 8, 9
        typelist = [Float16 Float32 Float64]
        for type in typelist
            # Test different delta values
            delta_set = [1.0 0.5 0.0]
            sol_set = [1/110 2/115 1/20]
            for delta in 1:3
                label = delta_set[delta]
                @testset "δ = $label" begin
                    vel = JuSwalbe.Twovector(x=fill(type(0.1), (N,M)), y=fill(type(-0.1), (N,M)))
                    zerovec = JuSwalbe.Twovector(x=zeros(type ,(N,M)), y=zeros(type ,(N,M)))
                    moment = JuSwalbe.Macroquant(height=ones(type, (N,M)), velocity=vel, pressure=zeros(type, (N,M)), energy=zeros(type, (N,M)))
                    forces = JuSwalbe.Forces(slip=zerovec, thermal=zerovec, h∇p=zerovec, bathymetry=zerovec)
                    input = JuSwalbe.Inputconstants(δ=type(delta_set[delta]), τ=type(1.0))
                    slippage = computeslip(moment, forces, input)
                    
                    for i in 1:N, j in 1:M
                        @test slippage.x[i,j] - type(sol_set[delta]) ≈ type(0.0) atol=tolerances[type]
                        @test slippage.y[i,j] + type(sol_set[delta]) ≈ type(0.0) atol=tolerances[type]
                    end
                end
            end
        end
    end
    @testset "Thermo capillary forcing two dimensions" begin
        N, M = 8, 9
        typelist = [Float16 Float32 Float64]
        for type in typelist
            # Test different delta values
            mom = simplemoment2d(N, M; T=type)
            force = JuSwalbe.Forces(slip=JuSwalbe.Twovector(x=fill(type(0.1), (N,M)), y=fill(type(0.1), (N,M))), h∇p = simpleTwovector(N, M; T=type),
                                    thermal=JuSwalbe.Twovector(x=zeros(type, (N,M)), y=zeros(type, (N,M))), bathymetry=JuSwalbe.Twovector(x=zeros(type, (N,M)), y=zeros(type, (N,M))))
            @test isa(force, JuSwalbe.Forces{JuSwalbe.Twovector{Matrix{type}}})
            @test isa(force.thermal.x, Array{type, 2})
            @test isa(force.thermal.y, Array{type, 2})
            @test size(force.thermal.x) == (N,M)
            @test size(force.thermal.y) == (N,M)
        end
    end
    @testset "Thermo capillary forcing one dimension" begin
        N = 8
        typelist = [Float16 Float32 Float64]
        for type in typelist
            # Test different delta values
            mom = simplemoment1d(N; T=type)
            force = JuSwalbe.Forces(slip=fill(type(0.1), N) , h∇p = zeros(type, N),
                                    thermal=zeros(type, N), bathymetry=zeros(type, N))
            @test isa(force, JuSwalbe.Forces{Vector{type}})
            @test isa(force.thermal, Array{type, 1})
            @test length(force.thermal) == N
        end
    end
end
