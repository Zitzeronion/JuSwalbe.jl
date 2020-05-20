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

julia> isfile("data/tmp/checkpoint_0.bson")
true

julia> load_dist = loadcheckpoint("data/tmp/checkpoint_0.bson")
JuSwalbe.DistributionD1Q3{Array{Float64,1}}
  f0: Array{Float64}((10,)) [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  f1: Array{Float64}((10,)) [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
  f2: Array{Float64}((10,)) [-0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1]

julia> load_dist.f0 == dist.f0; load_dist.f1 == dist.f1; load_dist.f2 == dist.f2
true

julia> rm("data/tmp/checkpoint_0.bson")
```

See also: [`loadcheckpoint`](@ref)

# Additional information
The size of the checkpoint heavily depends on the system size.
In case of two spatial dimensions and a 4096^2 lattice a single checkpoint file for a Float64 simulation is about 1.2GB data.

"""
function savecheckpoint(dist::JuSwalbe.Distributionfunction; name::String="tmp", t::Int=0)
    # Generate a path for the checkpoint and a file name
    outpath = mkpath("data/$name")
    location = outpath * "/checkpoint_" * "$t" * ".bson"
    # BSON saves only dictonaries, so transform the distribution to a dictonary
    distdict = struct2dict(dist)
    # Remove all previous checkpoints
    rm.(glob("data/$name/checkpoint_*"))
    # Save the distribution and generate a checkpoint file
    bson(location, distdict)
end

"""
    savecheckpoints(dist, name, t)

Same as savecheckpoint but does not delete the previous checkpoint files.

See also [`savecheckpoint`](@ref), [`loadcheckpoint`](@ref)
"""
function savecheckpoints(dist::JuSwalbe.Distributionfunction; name::String="tmp", t::Int=0)
    # Generate a path for the checkpoint
    outpath = mkpath("data/$name")
    location = outpath * "/checkpoint_" * "$t" * ".bson"
    # BSON saves only dictonaries, so transform the distribution to a dictonary
    distdict = struct2dict(dist)

    # Save the distribution and generate a checkpoint file
    bson(location, distdict)
end

"""
    loadcheckpoint(file)

Loads a checkpoint file and generates a distribution function from it.

With the distribution function one can rerun a simulation from the point where the distribution function was saved.
This is helpful for the long runs, when it takes several hours to generate convergence.
Such it is adviced to generate a checkpoint roughly every three hours.

# Example
For an example see `savecheckpoint`.

See also: [`savecheckpoint`](@ref), [`savecheckpoints`](@ref)
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

"""
    height2file(mom, name, time)

Saves the heightfield to a file in a binary format.

For data purposes and reproducibility one can save the heightfield to a file.
Such the file can be accessed with another software for e.g. post processing.
    
# Example
```jldoctest
julia> using JuSwalbe

julia> mom = JuSwalbe.Macroquant(height = ones(5,5), velocity=JuSwalbe.Twovector(x=zeros(5,5), y=zeros(5,5)), pressure=zeros(5,5), energy=zeros(5,5))
JuSwalbe.Macroquant{Array{Float64,2},JuSwalbe.Twovector{Array{Float64,2}}}
  height: Array{Float64}((5, 5)) [1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0; … ; 1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0]
  velocity: JuSwalbe.Twovector{Array{Float64,2}}
  pressure: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  energy: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]

julia> height2file(mom)

julia> h = BSON.load("data/tmp/height_0.bson")
Dict{Symbol,Any} with 1 entry:
  :height => [1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0; … ; 1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0]

julia> h[:height]
5×5 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0

julia> mom.height == h[:height]
true
```

See also: [`savecheckpoint`](@ref), [`velocity2file`](@ref), [`velocityandheight2file`](@ref)
"""
function height2file(mom::Macroscopic_quantity; name::String="tmp", time::Int=0)
    outpath = mkpath("data/$name")
    location = outpath * "/height_" * "$time" * ".bson"
    # println("Saved checkpoint @ $location")
    bson(location, height = mom.height)
end

"""
    velocity2file(mom, name, time)

Saves the velocity field to a binary file called velocity_t.bson.

Independent of the dimensionality of the system, `D1Q3`, `D2Q9` saves the velocity to a file.
In case of `D2Q9` saves both velocity components to the file such `velocity.x` and `velocity.y`.

# Example
```jldoctest
julia> using JuSwalbe

julia> mom = simplemoment2d(5,5)
JuSwalbe.Macroquant{Array{Float64,2},JuSwalbe.Twovector{Array{Float64,2}}}
  height: Array{Float64}((5, 5)) [1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0; … ; 1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0]
  velocity: JuSwalbe.Twovector{Array{Float64,2}}
  pressure: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  energy: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]


```
"""
function velocity2file(mom::Macroquant{Vector{T}, Vector{T}}; name::String="tmp", time::Int=0) where {T<:Number}
    outpath = mkpath("../data/$name")
    location = outpath * "/velocity_" * "$time" * ".bson"
    
    bson(location, velocity = mom.velocity)
    
end
function velocity2file(mom::Macroquant{Matrix{T}, JuSwalbe.Twovector{Matrix{T}}}; name::String="tmp", time::Int=0) where {T<:Number}
    outpath = mkpath("../data/$name")
    location = outpath * "/velocity_" * "$time" * ".bson"
        
    bson(location, velocity_x = mom.velocity.x, velocity_y = mom.velocity.y)
    
end

function velocityandheight2file(mom::Macroquant{Vector{T}, Vector{T}}; name::String="tmp", time::Int=0) where {T<:Number}
    outpath = mkpath("../data/$name")
    location = outpath * "/height_velocity_" * "$time" * ".bson"
    
    bson(location, height = mom.height, velocity = mom.velocity)

end
function velocityandheight2file(mom::Macroquant{Matrix{T}, JuSwalbe.Twovector{Matrix{T}}}; name::String="tmp", time::Int=0) where {T<:Number}
    outpath = mkpath("../data/$name")
    location = outpath * "/height_velocity_" * "$time" * ".bson"
    
    bson(location, height = mom.height, velocity_x = mom.velocity.x, velocity_y = mom.velocity.y)

end
