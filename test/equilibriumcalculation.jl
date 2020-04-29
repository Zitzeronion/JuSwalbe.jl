@testset "Moments" begin
    N = 10
    moments = JuSwalbe.macroquant64_1d(zeros(Float64,N),zeros(Float64,N),zeros(Float64,N))
    @test typeof(moments.height) == Array{Float64,1}
    @test typeof(moments.velocity) == Array{Float64,1}
    @test typeof(moments.energy) == Array{Float64,1} 
    @test size(moments.height, 1) == N
end