using JuSwalbe
using Test, DrWatson, DelimitedFiles

# Pkg.test("JuSwalbe", test_args=["call"])
@testset "JuSwalbe.jl" begin
    include("call.jl")
    include("macroscopic_typs.jl")
    include("distribution_types.jl")
    include("readinput.jl")
    include("equilibriumcalculation.jl")
end
