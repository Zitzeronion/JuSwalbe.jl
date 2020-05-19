"""
    savecheckpoint(dist, t)

Saves the distribution function to a binary file. 

This file can be seen as a checkpoint or snapshot of the simulations internal state and be used to restart the simulation.

# Example
```jldoctest
julia> using JuSwalbe

julia> dist = JuSwalbe.DistributionD1Q3(f0 = ones(10), f1 = fill(0.1, 10), f2 = fill(-0.1, 10))
JuSwalbe.DistributionD1Q3{Array{Float64,1}}
  f0: Array{Float64}((10,)) [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  f1: Array{Float64}((10,)) [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
  f2: Array{Float64}((10,)) [-0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1]

julia> savecheckpoint(dist)

julia> isfile("data/checkpoint_tmp_0.bson")
true

julia> load_dist = loadcheckpoint("data/checkpoint_tmp_0.bson")
JuSwalbe.DistributionD1Q3{Array{Float64,1}}
  f0: Array{Float64}((10,)) [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  f1: Array{Float64}((10,)) [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
  f2: Array{Float64}((10,)) [-0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1]

julia> load_dist.f0 == dist.f0; load_dist.f1 == dist.f1; load_dist.f2 == dist.f2
true

julia> rm("data/checkpoint_tmp_0.bson")
```

See also: [`loadcheckpoint`](@ref)
"""
function savecheckpoint(dist::JuSwalbe.Distributionfunction; name::String="tmp", t::Int=0)
    outpath = mkpath("../data/$name")
    distdict = struct2dict(dist)
    location = outpath * "/checkpoint_" * "$t" * ".bson"
    # println("Saved checkpoint @ $location")
    bson(location, distdict)
end

"""
    loadcheckpoint(file)

Loads a checkpoint file and generates a distribution function from it.

See also: [`savecheckpoint`](@ref)
"""
function loadcheckpoint(file)
    # Save the binary file to a dictonary
    distdict = BSON.load(file)
    # Check if it is one or two dimensional
    # and create a distribution 
    if length(keys(distdict)) == 3
        dist = JuSwalbe.DistributionD1Q3(f0 = distdict[:f0], f1 = distdict[:f1], f2 = distdict[:f2])
    elseif length(keys(distdict)) == 9
        dist = JuSwalbe.DistributionD2Q9(f0 = distdict[:f0], f1 = distdict[:f1], f2 = distdict[:f2],
                                         f3 = distdict[:f3], f4 = distdict[:f4], f5 = distdict[:f5],
                                         f6 = distdict[:f6], f7 = distdict[:f7], f8 = distdict[:f8])
    end
    # Lastly return the distribution
    return dist
end