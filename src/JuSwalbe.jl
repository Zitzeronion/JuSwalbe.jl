module JuSwalbe

using Pkg
using DelimitedFiles, JSON, DrWatson, Images, Parameters, BSON, Glob

include("minimalsetup.jl")
include("macroscopic_typs.jl")
include("distribution_types.jl")
include("readinput.jl")
include("equilibriumcalculation.jl")
include("calculatepressure.jl")
include("calcmoments.jl")
include("forcing.jl")
include("collide.jl")
include("readwritedata.jl")

# IO related
export readinput, findargument, minimalsetup1d, minimalsetup2d, savecheckpoint, loadcheckpoint
# Pressure related 
export Δh, Π, pressure, ∇p
# Distributions functions and moments
export calc_equilibrium_distribution, calculatemoments, dist2array, collisionBGK, streamdistperiodic!
# Forces
export computeslip

end # module
