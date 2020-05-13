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
include("slippage.jl")

# IO related
export call, readinput
# Pressure related 
export Δh, Π, pressure, ∇p
# Distributions functions and moments
export calc_equilibrium_distribution, calculatemoments, dist2array
# Forces
export computeslip

end # module
