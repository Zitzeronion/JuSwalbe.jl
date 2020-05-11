@testset "Moment calculation" begin
    @testset "D2Q9" begin
        @testset "With velocity in ux" begin
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testarray1 = fill(type(0.1), (5,5))
                testarray2 = fill(type(0.2), (5,5))
                testdist = JuSwalbe.DistributionD2Q9(testarray1, testarray2, testarray1, testarray1, testarray1, testarray1, testarray1, testarray1, testarray1)
                @test isa(testdist, JuSwalbe.DistributionD2Q9)
                mom = calculatemoments(testdist)
                @test isa(mom, JuSwalbe.Macroquant)
                @test size(mom.height) = (5,5)
                @test size(mom.velocity.x) = (5,5)
                @test size(mom.velocity.y) = (5,5)
                for i in 1:5, j in 1:5
                    @test (mom.height[i,j] - type(1.0)) ≈ 0.0 atol = 1e-5 
                    @test (mom.velocity.x[i,j] - type(0.1)) ≈ 0.0 atol = 1e-5
                    @test (mom.velocity.y[i,j] - type(0.0)) ≈ 0.0 atol = 1e-5  
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
                @test size(mom.height) = (5,5)
                @test size(mom.velocity.x) = (5,5)
                @test size(mom.velocity.y) = (5,5)
                for i in 1:5, j in 1:5
                    @test (mom.height[i,j] - type(1.0)) ≈ 0.0 atol = 1e-5 
                    @test (mom.velocity.x[i,j] - type(0.0)) ≈ 0.0 atol = 1e-5
                    @test (mom.velocity.y[i,j] - type(0.1)) ≈ 0.0 atol = 1e-5  
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
                @test size(mom.height) = (5,5)
                @test size(mom.velocity.x) = (5,5)
                @test size(mom.velocity.y) = (5,5)
                for i in 1:5, j in 1:5
                    @test (mom.height[i,j] - type(1.1)) ≈ 0.0 atol = 1e-5 
                    @test (mom.velocity.x[i,j] - type(0.1)) ≈ 0.0 atol = 1e-5
                    @test (mom.velocity.y[i,j] - type(0.1)) ≈ 0.0 atol = 1e-5  
                end
            end
        end

        @testset "Without velocity" begin
            typelist = [Float16 Float32 Float64]
            for type in typelist
                testarray1 = fill(type(0.1), (5,5))
                testarray2 = fill(type(0.2), (5,5))
                testdist = JuSwalbe.DistributionD2Q9(testarray1, testarray2, testarray2, testarray1, testarray1, testarray1, testarray1, testarray1, testarray1)
                @test isa(testdist, JuSwalbe.DistributionD2Q9)
                mom = calculatemoments(testdist)
                @test isa(mom, JuSwalbe.Macroquant)
                @test size(mom.height) = (5,5)
                @test size(mom.velocity.x) = (5,5)
                @test size(mom.velocity.y) = (5,5)
                for i in 1:5, j in 1:5
                    @test (mom.height[i,j] - type(1.1)) ≈ 0.0 atol = 1e-5 
                    @test (mom.velocity.x[i,j] - type(0.1)) ≈ 0.0 atol = 1e-5
                    @test (mom.velocity.y[i,j] - type(0.1)) ≈ 0.0 atol = 1e-5  
                end
            end
        end
    end
end