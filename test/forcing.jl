@testset "Forcing" begin
    tolerances = Dict(Float16 => 1e-4, Float32 => 1e-5, Float64 => 1e-7, Real => 1e-7)
    @testset "Slippage one dimension" begin
        N = 8
        typelist = [Float16 Float32 Float64]
        for type in typelist
            @testset "δ = 1" begin
                moment = JuSwalbe.Macroquant(height=ones(type, N), velocity=fill(type(0.1), N), pressure=zeros(type, N), energy=zeros(type, N))
                forces = JuSwalbe.Forces(slip=zeros(type, N), thermal=zeros(type, N), h∇p=zeros(type, N), bathymetry=zeros(type, N))
                input = JuSwalbe.Inputconstants(δ=type(1.0), τ=type(1.0))
                slippage = computeslip(moment, forces, input)
                
                for i in 1:N
                    @test slippage[i] - type(1/110) ≈ type(0.0) atol=tolerances[type]
                end
            end
            # Different δ
            @testset "δ = 0.5" begin
                moment = JuSwalbe.Macroquant(height=ones(type, N), velocity=fill(type(0.1), N), pressure=zeros(type, N), energy=zeros(type, N))
                forces = JuSwalbe.Forces(slip=zeros(type, N), thermal=zeros(type, N), h∇p=zeros(type, N), bathymetry=zeros(type, N))
                input = JuSwalbe.Inputconstants(δ=type(0.5), τ=type(1.0))
                slippage = computeslip(moment, forces, input)
                
                for i in 1:N
                    @test slippage[i] - type(2/115) ≈ type(0.0) atol=tolerances[type]
                end
            end
            # No Slip case
            @testset "δ = 0" begin
                moment = JuSwalbe.Macroquant(height=ones(type, N), velocity=fill(type(0.1), N), pressure=zeros(type, N), energy=zeros(type, N))
                forces = JuSwalbe.Forces(slip=zeros(type, N), thermal=zeros(type, N), h∇p=zeros(type, N), bathymetry=zeros(type, N))
                input = JuSwalbe.Inputconstants(δ=type(0.0), τ=type(1.0))
                slippage = computeslip(moment, forces, input)
                
                for i in 1:N
                    @test slippage[i] - type(1/20) ≈ type(0.0) atol=tolerances[type]
                end
            end
            
        end
    end
    @testset "Slippage two dimensions" begin
        N, M = 8, 9
        typelist = [Float16 Float32 Float64]
        for type in typelist
            @testset "δ = 1" begin
                moment = JuSwalbe.Macroquant(height=ones(type, N), velocity=fill(type(0.1), N), pressure=zeros(type, N), energy=zeros(type, N))
                forces = JuSwalbe.Forces(slip=zeros(type, N), thermal=zeros(type, N), h∇p=zeros(type, N), bathymetry=zeros(type, N))
                input = JuSwalbe.Inputconstants(δ=type(1.0), τ=type(1.0))
                slippage = computeslip(moment, forces, input)
                
                for i in 1:N
                    @test slippage[i] - type(1/110) ≈ type(0.0) atol=tolerances[type]
                end
            end
            # Different δ
            @testset "δ = 0.5" begin
                moment = JuSwalbe.Macroquant(height=ones(type, N), velocity=fill(type(0.1), N), pressure=zeros(type, N), energy=zeros(type, N))
                forces = JuSwalbe.Forces(slip=zeros(type, N), thermal=zeros(type, N), h∇p=zeros(type, N), bathymetry=zeros(type, N))
                input = JuSwalbe.Inputconstants(δ=type(0.5), τ=type(1.0))
                slippage = computeslip(moment, forces, input)
                
                for i in 1:N
                    @test slippage[i] - type(2/115) ≈ type(0.0) atol=tolerances[type]
                end
            end
            # No Slip case
            @testset "δ = 0" begin
                moment = JuSwalbe.Macroquant(height=ones(type, N), velocity=fill(type(0.1), N), pressure=zeros(type, N), energy=zeros(type, N))
                forces = JuSwalbe.Forces(slip=zeros(type, N), thermal=zeros(type, N), h∇p=zeros(type, N), bathymetry=zeros(type, N))
                input = JuSwalbe.Inputconstants(δ=type(0.0), τ=type(1.0))
                slippage = computeslip(moment, forces, input)
                
                for i in 1:N
                    @test slippage[i] - type(1/20) ≈ type(0.0) atol=tolerances[type]
                end
            end
            
        end
    end    
end
