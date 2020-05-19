using JuSwalbe
using Test, DrWatson, DelimitedFiles, Images

@testset "JuSwalbe.jl" begin
    @testset "LBM Parts" begin
        include("macroscopic_typs.jl")
        include("distribution_types.jl")
        include("readinput.jl")
        include("equilibriumcalculation.jl")
        include("calculatepressure.jl")
        include("calcmoments.jl")
        include("forcing.jl")
        include("collide.jl")
    end
    @testset "IO and Others" begin
        include("minimalsetup.jl")
        include("readwritedata.jl")
    end
end


