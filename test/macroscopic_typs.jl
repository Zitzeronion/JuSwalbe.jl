@testset "Datastructures" begin
    # Test the types and the sizes of the moments
    # For 64 bit floating numbers, one dimensional
    @testset "One dimensional moments" begin
        n = 10
        typelist = [Float16 Float32 Float64 Real]
        for type in typelist
            moments = JuSwalbe.Macroquant{Vector{type},Vector{type}}(zeros(type,n),zeros(type,n),zeros(type,n),zeros(type,n))
            @test isa(moments, JuSwalbe.Macroquant)
            test_dict = struct2dict(moments)
            for key in keys(test_dict)
                @test isa(test_dict[key], Array{type,1})
                @test size(test_dict[key], 1) == n
            end
        end
    end
    @testset "One dimensional forces" begin
        n = 10
        typelist = [Float16 Float32 Float64 Real]
        for type in typelist
            forces = JuSwalbe.Forces{Vector{type}}(zeros(type,n),zeros(type,n),zeros(type,n),zeros(type,n))
            @test isa(forces, JuSwalbe.Forces)
            test_dict = struct2dict(forces)
            for key in keys(test_dict)
                @test isa(test_dict[key], Array{type,1})
                @test size(test_dict[key], 1) == n
            end
        end
    end
    # For 32 bit floating numbers, one dimensional
    @testset "Two dimensional moments" begin
        m,n = 10,8
        typelist = [Float16 Float32 Float64 Real]
        for type in typelist
            twovec = JuSwalbe.Twovector{Matrix{type}}(zeros(type, (m,n)), zeros(type, (m,n)))
            moments = JuSwalbe.Macroquant{Matrix{type},JuSwalbe.Twovector{Matrix{type}}}(zeros(type, (m,n)),twovec ,zeros(type, (m,n)),zeros(type, (m,n)))
            @test isa(moments, JuSwalbe.Macroquant)
            test_dict = struct2dict(moments)
            for key in keys(test_dict)
                if key == :velocity
                    @test isa(test_dict[key], JuSwalbe.Twovector)
                    @test isa(test_dict[key].x, Matrix{type})
                    @test isa(test_dict[key].y, Matrix{type})
                    @test size(test_dict[key].x) == (m,n)
                    @test size(test_dict[key].y) == (m,n)
                else
                    @test isa(test_dict[key], Array{type,2})
                    @test size(test_dict[key], 1) == m
                    @test size(test_dict[key], 2) == n
                end
            end
        end
    end
    @testset "Two dimensional forces" begin
        n, m = 10, 6
        typelist = [Float16 Float32 Float64 Real]
        for type in typelist
            twovec = JuSwalbe.Twovector{Matrix{type}}(zeros(type, (n,m)),zeros(type, (n,m))) 
            forces = JuSwalbe.Forces{JuSwalbe.Twovector{Matrix{type}}}(twovec, twovec, twovec, twovec)
            @test isa(forces, JuSwalbe.Forces)
            test_dict = struct2dict(forces)
            for key in keys(test_dict)
                @test isa(test_dict[key], JuSwalbe.Twovector)
                @test isa(test_dict[key].x, Matrix{type})
                @test isa(test_dict[key].y, Matrix{type})
                @test size(test_dict[key].x) == (n,m)
                @test size(test_dict[key].y) == (n,m)
            end
        end
    end
end