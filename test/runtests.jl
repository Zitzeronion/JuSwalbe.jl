using JuSwalbe
using Test

# Pkg.test("JuSwalbe", test_args=["call"])
@testset "JuSwalbe.jl" begin
    include("call.jl")
    include("equilibriumcalculation.jl")
end
