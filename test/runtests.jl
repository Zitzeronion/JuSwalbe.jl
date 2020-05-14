using JuSwalbe
using Test, DrWatson, DelimitedFiles, Images

@testset "JuSwalbe.jl" begin
    include("call.jl")
    include("macroscopic_typs.jl")
    include("distribution_types.jl")
    include("readinput.jl")
    include("equilibriumcalculation.jl")
    include("calculatepressure.jl")
    include("calcmoments.jl")
    include("forcing.jl")
end
