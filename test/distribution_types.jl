@testset "Distribution types" begin
    # Test the types and the sizes of the distributions
    @testset "Distributions D1Q3" begin
        N = 10 
        typelist = [Float16 Float32 Float64 Real]
        for type in typelist
            test_distribution = JuSwalbe.DistributionD1Q3(zeros(type,N),zeros(type,N),zeros(type,N))
            test_dict = struct2dict(test_distribution)
            @test isa(test_distribution, JuSwalbe.DistributionD1Q3)
            for key in keys(test_dict)
                @test isa(test_dict[key], Array{type,1}) 
                @test length(test_dict[key]) == N
            end
        end
    end
    @testset "Distributions D2Q9" begin
        N, M = 10, 14 
        typelist = [Float16 Float32 Float64 Real]
        for type in typelist
            test_distribution = JuSwalbe.DistributionD2Q9(zeros(type,(N,M)),
                                                          zeros(type,(N,M)),
                                                          zeros(type,(N,M)),
                                                          zeros(type,(N,M)),
                                                          zeros(type,(N,M)),
                                                          zeros(type,(N,M)),
                                                          zeros(type,(N,M)),
                                                          zeros(type,(N,M)),
                                                          zeros(type,(N,M)))
            test_dict = struct2dict(test_distribution)
            @test isa(test_distribution, JuSwalbe.DistributionD2Q9)
            for key in keys(test_dict)
                @test isa(test_dict[key], Array{type,2}) 
                @test size(test_dict[key]) == (N,M)
            end
        end
    end
end