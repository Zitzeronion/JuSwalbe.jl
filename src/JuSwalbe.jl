module JuSwalbe

    using Pkg
    using DelimitedFiles, JSON, DrWatson, Images, Parameters, BSON, Glob, Distributions
    using Revise

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
    export readinput, findargument, minimalsetup1d, minimalsetup2d, simplemoment2d, simpleTwovector, simpledistD2Q9, simplemoment1d
    # File IO
    export savecheckpoint, savecheckpoints, loadcheckpoint, height2file, velocity2file, velocityandheight2file
    # Pressure related 
    export Δh, Π, pressure, ∇p
    # Distributions functions and moments
    export calc_equilibrium_distribution, calculatemoments, dist2array, collisionBGK, streamdistperiodic!
    # Forces
    export computeslip, computethermalcapillarywaves

end # module
