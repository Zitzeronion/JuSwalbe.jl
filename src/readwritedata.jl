"""
    savecheckpoint(dist, t)

Saves the distribution function to a binary file. 

This file can be seen as a checkpoint or snapshot of the simulations internal state.
"""
function savecheckpoint(dist::JuSwalbe.Distributionfunction; name::String="tmp", t::Int=0)
    distdict = struct2dict(dist)
    location = "data/checkpoint_" * name * "_$t" * ".bson"
    bson(location, distdict)
end

"""
    loadcheckpoint(file)

Loads a checkpoint file and generates a distribution function from it.
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