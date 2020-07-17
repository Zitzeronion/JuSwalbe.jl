module JuSwalbe
    import Base.Threads.@spawn
    using Pkg
    using DelimitedFiles, JSON, DrWatson, Images, Parameters, BSON, Glob, Distributions, Makie
    using Revise, CUDA

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
    include("startrun.jl")

    # IO related
    export readinput, findargument, minimalsetup1d, minimalsetup2d, simplemoment2d, simpleTwovector, zeroTwovector, simpledistD2Q9, simplemoment1d
    # File IO
    export savecheckpoint, savecheckpoints, loadcheckpoint, height2file, velocity2file, velocityandheight2file
    # Pressure related 
    export Δh, Π, pressure, ∇p
    # Distributions functions and moments
    export calc_equilibrium_distribution, updateequilibrium!, calculatemoments, dist2array, collisionBGK, streamdistperiodic!, streamdistperiodic
    # Forces
    export computeslip, computethermalcapillarywaves

    # Full simulations
    # export runsimulation_1D, runsimulation_2D

end # module
