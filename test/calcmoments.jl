@testset "Moment calculation" begin
    tolerances = Dict(Float16 => 1e-3, Float32 => 1e-5, Float64 => 1e-7, Real => 1e-7)
    @testset "D2Q9" begin
        @testset "With Guo Forcing but no forces" begin
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testarray1 = fill(type(0.1), (5,5))
                testarray2 = fill(type(0.2), (5,5))
                testdist = JuSwalbe.DistributionD2Q9(f0=testarray1, f1=testarray2, f2=testarray1, f3=testarray1, f4=testarray1, f5=testarray1, f6=testarray1, f7=testarray1, f8=testarray1)
                @test isa(testdist, JuSwalbe.DistributionD2Q9)
                zerovec = JuSwalbe.Twovector(x=zeros(type, (5,5)), y=zeros(type, (5,5)))
                testforce = JuSwalbe.Forces(slip=zerovec, thermal=zerovec, h∇p=zerovec, bathymetry=zerovec)
                mom = calculatemoments(testdist, testforce)
                @test isa(mom, JuSwalbe.Macroquant)
                @test size(mom.height) == (5,5)
                @test size(mom.velocity.x) == (5,5)
                @test size(mom.velocity.y) == (5,5)
                for i in 1:5, j in 1:5
                    @test (mom.height[i,j] - type(1.0)) ≈ type(0.0) atol = tolerances[type] 
                    @test (mom.velocity.x[i,j] - type(0.1/1.0)) ≈ type(0.0) atol = tolerances[type]
                    @test (mom.velocity.y[i,j] - type(0.0/1.0)) ≈ type(0.0) atol = tolerances[type]  
                end
            end
        end
        @testset "With velocity in ux" begin
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testarray1 = fill(type(0.1), (5,5))
                testarray2 = fill(type(0.2), (5,5))
                testdist = JuSwalbe.DistributionD2Q9(testarray1, testarray2, testarray1, testarray1, testarray1, testarray1, testarray1, testarray1, testarray1)
                @test isa(testdist, JuSwalbe.DistributionD2Q9)
                mom = calculatemoments(testdist)
                @test isa(mom, JuSwalbe.Macroquant)
                @test size(mom.height) == (5,5)
                @test size(mom.velocity.x) == (5,5)
                @test size(mom.velocity.y) == (5,5)
                for i in 1:5, j in 1:5
                    @test (mom.height[i,j] - type(1.0)) ≈ type(0.0) atol = tolerances[type] 
                    @test (mom.velocity.x[i,j] - type(0.1/1.0)) ≈ type(0.0) atol = tolerances[type]
                    @test (mom.velocity.y[i,j] - type(0.0/1.0)) ≈ type(0.0) atol = tolerances[type]  
                end
            end
        end
        @testset "With velocity in uy" begin
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testarray1 = fill(type(0.1), (5,5))
                testarray2 = fill(type(0.2), (5,5))
                testdist = JuSwalbe.DistributionD2Q9(testarray1, testarray1, testarray2, testarray1, testarray1, testarray1, testarray1, testarray1, testarray1)
                @test isa(testdist, JuSwalbe.DistributionD2Q9)
                mom = calculatemoments(testdist)
                @test isa(mom, JuSwalbe.Macroquant)
                @test size(mom.height) == (5,5)
                @test size(mom.velocity.x) == (5,5)
                @test size(mom.velocity.y) == (5,5)
                for i in 1:5, j in 1:5
                    @test (mom.height[i,j] - type(1.0)) ≈ type(0.0) atol = tolerances[type] 
                    @test (mom.velocity.x[i,j] - type(0.0/1.0)) ≈ type(0.0) atol = tolerances[type]
                    @test (mom.velocity.y[i,j] - type(0.1/1.0)) ≈ type(0.0) atol = tolerances[type]  
                end
            end
        end
        @testset "With velocity in ux and uy" begin
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testarray1 = fill(type(0.1), (5,5))
                testarray2 = fill(type(0.2), (5,5))
                testdist = JuSwalbe.DistributionD2Q9(testarray1, testarray2, testarray2, testarray1, testarray1, testarray1, testarray1, testarray1, testarray1)
                @test isa(testdist, JuSwalbe.DistributionD2Q9)
                mom = calculatemoments(testdist)
                @test isa(mom, JuSwalbe.Macroquant)
                @test size(mom.height) == (5,5)
                @test size(mom.velocity.x) == (5,5)
                @test size(mom.velocity.y) == (5,5)
                for i in 1:5, j in 1:5
                    @test (mom.height[i,j] - type(1.1)) ≈ type(0.0) atol = tolerances[type] 
                    @test (mom.velocity.x[i,j] - type(0.1/1.1)) ≈ type(0.0) atol = tolerances[type]
                    @test (mom.velocity.y[i,j] - type(0.1/1.1)) ≈ type(0.0) atol = tolerances[type]  
                end
            end
        end
        @testset "Without velocity" begin
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testarray1 = fill(type(0.1), (5,5))
                testarray2 = fill(type(0.2), (5,5))
                testdist = JuSwalbe.DistributionD2Q9(testarray1, testarray1, testarray1, testarray1, testarray1, testarray1, testarray1, testarray1, testarray1)
                @test isa(testdist, JuSwalbe.DistributionD2Q9)
                mom = calculatemoments(testdist)
                @test isa(mom, JuSwalbe.Macroquant)
                @test size(mom.height) == (5,5)
                @test size(mom.velocity.x) == (5,5)
                @test size(mom.velocity.y) == (5,5)
                for i in 1:5, j in 1:5
                    @test (mom.height[i,j] - type(0.9)) ≈ type(0.0) atol = tolerances[type] 
                    @test (mom.velocity.x[i,j] - type(0.0)) ≈ type(0.0) atol = tolerances[type]
                    @test (mom.velocity.y[i,j] - type(0.0)) ≈ type(0.0) atol = tolerances[type]  
                end
            end
        end
    end

    @testset "D1Q3" begin
        @testset "With velocity" begin
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testarray1 = fill(type(0.1), 10)
                testarray2 = fill(type(0.2), 10)
                testdist = JuSwalbe.DistributionD1Q3(testarray1, testarray2, testarray1)
                @test isa(testdist, JuSwalbe.DistributionD1Q3)
                mom = calculatemoments(testdist)
                @test isa(mom, JuSwalbe.Macroquant)
                @test length(mom.height) == 10
                @test length(mom.velocity) == 10
                for i in 1:10
                    @test (mom.height[i] - type(0.4)) ≈ 0.0 atol = tolerances[type] 
                    @test (mom.velocity[i] - type(0.1/0.4)) ≈ 0.0 atol = tolerances[type]  
                end
            end
        end
        @testset "Without velocity" begin
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testarray1 = fill(type(0.1), 10)
                testarray2 = fill(type(0.2), 10)
                testdist = JuSwalbe.DistributionD1Q3(testarray1, testarray1, testarray1)
                @test isa(testdist, JuSwalbe.DistributionD1Q3)
                mom = calculatemoments(testdist)
                @test isa(mom, JuSwalbe.Macroquant)
                @test length(mom.height) == 10
                @test length(mom.velocity) == 10
                for i in 1:10
                    @test (mom.height[i] - type(0.3)) ≈ 0.0 atol = tolerances[type] 
                    @test (mom.velocity[i] - type(0.0)) ≈ 0.0 atol = tolerances[type]
                end
            end
        end
    end

    @testset "Distribution to Array" begin
        typelist = [Float16 Float32 Float64]
        for type in typelist
            testarray1 = fill(type(0.1), (5,5))
            testarray2 = fill(type(0.2), (5,5))
            testdist = JuSwalbe.DistributionD2Q9(testarray1, testarray1, testarray2, testarray1, testarray1, testarray1, testarray1, testarray1, testarray1)
            @test isa(testdist, JuSwalbe.DistributionD2Q9)
            arr = dist2array(testdist)
            @test isa(arr, Array{type, 3})
            @test size(arr) == (9,5,5)
            for i in 1:5, j in 1:5
                @test testarray1[i,j] == arr[1,i,j]
                @test testarray2[i,j] == arr[3,i,j]
            end
        end
    end
end

