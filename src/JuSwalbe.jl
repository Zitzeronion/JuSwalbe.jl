module JuSwalbe

using Pkg
using DelimitedFiles, JSON, DrWatson, Images, Parameters

include("call.jl")
include("macroscopic_typs.jl")
include("distribution_types.jl")
include("readinput.jl")
include("equilibriumcalculation.jl")
include("calculatepressure.jl")
include("calcmoments.jl")
include("forcing.jl")
include("collide.jl")

function minimalsetup1d(N::Int)
    input = Inputconstants(lx=N)
    mom = Macroquant(height = ones(N), velocity=fill(0.1, N), pressure=zeros(N), energy=zeros(N))
    forces = Forces(slip=fill(0.1, N), thermal=zeros(N), h∇p=fill(-0.1, N), bathymetry=zeros(N))
    dist = DistributionD1Q3(f0 = ones(N), f1 = zeros(N), f2 = zeros(N))
    return input, mom, forces, dist
end

function minimalsetup2d(N::Int, M::Int)
    input = Inputconstants(lx=N, ly=M)
    vel = Twovector(x=fill(0.02, (N,M)), y=fill(-0.02, (N,M)))
    zerovec = Twovector(x=zeros(N,M), y=zeros(N,M))
    mom = Macroquant(height = ones(N,M), velocity=vel, pressure=zeros(N,M), energy=zeros(N,M))
    forces = Forces(slip=vel, thermal=zerovec, h∇p=vel, bathymetry=zerovec)
    dist = DistributionD2Q9(f0 = ones(N,M), f1 = zeros(N,M), f2 = zeros(N,M), 
                            f3 = zeros(N,M), f4 = zeros(N,M), f5 = zeros(N,M), 
                            f6 = zeros(N,M), f7 = zeros(N,M), f8 = zeros(N,M))
    return input, mom, forces, dist
end

# IO related
export call, readinput, findargument, minimalsetup1d, minimalsetup2d
# Pressure related 
export Δh, Π, pressure, ∇p
# Distributions functions and moments
export calc_equilibrium_distribution, calculatemoments, dist2array, collisionBGK
# Forces
export computeslip

end # module
