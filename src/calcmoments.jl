"""
    calculatemoments(dist)

Computes the lattice Boltzmann moments (h,v) from a distribution function `dist`.

Based on the input either a `D1Q3` or a `D2Q9` distribution this function performs a summation over the particle speeds.
In one spatial dimension the velocity set [-1, 0, 1] allows for three speeds in two for nine speeds.
See also: [`dist2array`](@ref)

# Math
The moments are calculated according to 

`` \\rho = \\sum_{i=0}^N f_i ``

`` \\rho\\mathbf{u} = \\sum_{i=0}^N \\mathbf{c}_i f_i ``

`` e = \\sum_{i=0}^N \\mathbf{c}_i \\mathbf{c}_i f_i ``

Where N is either 3 or 9. In one spatial dimension the velocity is of course not a vector of ``(x,y)``.


# Example
```jldoctest
julia> using JuSwalbe

julia> a = fill(0.1, (10, 10))
10×10 Array{Float64,2}:
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1

julia> b = fill(0.2, (10,10))
10×10 Array{Float64,2}:
 0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2
 0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2
 0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2
 0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2
 0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2
 0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2
 0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2
 0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2
 0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2
 0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2  0.2

julia> dist = JuSwalbe.DistributionD2Q9(f0=a,f1=b,f2=a,f3=a,f4=a,f5=a,f6=a,f7=a,f8=a)
JuSwalbe.DistributionD2Q9{Array{Float64,2}}
  f0: Array{Float64}((10, 10)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f1: Array{Float64}((10, 10)) [0.2 0.2 … 0.2 0.2; 0.2 0.2 … 0.2 0.2; … ; 0.2 0.2 … 0.2 0.2; 0.2 0.2 … 0.2 0.2]
  f2: Array{Float64}((10, 10)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f3: Array{Float64}((10, 10)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f4: Array{Float64}((10, 10)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f5: Array{Float64}((10, 10)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f6: Array{Float64}((10, 10)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f7: Array{Float64}((10, 10)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f8: Array{Float64}((10, 10)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]

julia> mom = calculatemoments(dist)
JuSwalbe.Macroquant{Array{Float64,2},JuSwalbe.Twovector{Array{Float64,2}}}
  height: Array{Float64}((10, 10)) [0.9999999999999999 0.9999999999999999 … 0.9999999999999999 0.9999999999999999; 0.9999999999999999 0.9999999999999999 … 0.9999999999999999 0.9999999999999999; … ; 0.9999999999999999 0.9999999999999999 … 0.9999999999999999 0.9999999999999999; 0.9999999999999999 0.9999999999999999 … 0.9999999999999999 0.9999999999999999]
  velocity: JuSwalbe.Twovector{Array{Float64,2}}
  pressure: Array{Float64}((10, 10)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  energy: Array{Float64}((10, 10)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]

julia> mom.height
10×10 Array{Float64,2}:
 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0

julia> mom.velocity.x
10×10 Array{Float64,2}:
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1  0.1

julia> mom.velocity.y
10×10 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> dist1d = JuSwalbe.DistributionD1Q3(f0=a[1:10],f1=b[1:10],f2=a[1:10]) # Vectorize a and b
JuSwalbe.DistributionD1Q3{Array{Float64,1}}
  f0: Array{Float64}((10,)) [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
  f1: Array{Float64}((10,)) [0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2]
  f2: Array{Float64}((10,)) [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

julia> mom1d = calculatemoments(dist1d)
JuSwalbe.Macroquant{Array{Float64,1},Array{Float64,1}}
  height: Array{Float64}((10,)) [0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]
  velocity: Array{Float64}((10,)) [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
  pressure: Array{Float64}((10,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  energy: Array{Float64}((10,)) [0.30000000000000004, 0.30000000000000004, 0.30000000000000004, 0.30000000000000004, 0.30000000000000004, 0.30000000000000004, 0.30000000000000004, 0.30000000000000004, 0.30000000000000004, 0.30000000000000004]

julia> mom1d.height
10-element Array{Float64,1}:
 0.4
 0.4
 0.4
 0.4
 0.4
 0.4
 0.4
 0.4
 0.4
 0.4


julia> mom1d.velocity
10-element Array{Float64,1}:
 0.25
 0.25
 0.25
 0.25
 0.25
 0.25
 0.25
 0.25
 0.25
 0.25

```

# References
## Code
See also: [`dist2array`](@ref)

## Two spatial dimensions
Here I go with Dellar
- [Nonhydrodynamic modes and *a priori* construction of shallow water lattice Boltzmann equations](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.036309)

## One spatial dimension
For the one dimensional simulation the paper by Chopard is quite nice
- [Study of the 1D lattice Boltzmann shallow water equation and its coupling to build a canal network.](https://www.sciencedirect.com/science/article/pii/S0021999110003372)

"""
function calculatemoments(dist::JuSwalbe.DistributionD1Q3{Vector{T}}) where {T<:Number}
    # Get a size for the output arrays
    len = length(dist.f0)
    # Allocate the output
    height = zeros(T, len)
    velocity = zeros(T, len)
    pressure = zeros(T, len)
    energy = zeros(T, len)
    # Lattice velocities
    c = [0 1 -1]
    csquare = [0 1 1]
    # Store all distributions in an array
    distarray = hcat(dist.f0, dist.f1, dist.f2) 
    height = sum(distarray, dims=2)[:, 1]
    velocity = sum(distarray .* c, dims=2)[:, 1]
    energy = sum(distarray .* csquare, dims=2)[:, 1]
    
    velocity ./= height

    result = JuSwalbe.Macroquant(height, velocity, pressure, energy)
    return result
end

function calculatemoments(dist::JuSwalbe.DistributionD1Q3{Vector{T}}, force::JuSwalbe.Forces{Vector{T}}) where {T<:Number}
    # Get a size for the output arrays
    len = length(dist.f0)
    # Allocate the output
    height = zeros(T, len)
    velocity = zeros(T, len)
    pressure = zeros(T, len)
    energy = zeros(T, len)
    # Lattice velocities
    c = [0 1 -1]
    csquare = [0 1 1]
    # Store all distributions in an array
    distarray = hcat(dist.f0, dist.f1, dist.f2)
    # Store forces in array 
    forcearray = hcat(force.slip, force.thermal, force.h∇p)
    # Moment summation is acutally performed here
    height = sum(distarray, dims=2)[:, 1] 
    velocity = sum(distarray .* c, dims=2)[:, 1] + T(0.5) * sum(forcearray, dims=2)[:, 1]
    energy = sum(distarray .* csquare, dims=2)[:, 1]
    # Last step is to divide the velocity by h
    velocity ./= height
    # Return the output
    result = JuSwalbe.Macroquant(height, velocity, pressure, energy)
    return result
end
# Below is the two dimensional thing
function calculatemoments(dist::JuSwalbe.DistributionD2Q9{Matrix{T}}) where {T<:Number}
    # Get a size for the output arrays
    width, thick = size(dist.f0)
    # Allocate the output
    height = zeros(T, (width,thick))
    velocityx = zeros(T, (width,thick))
    velocityy = zeros(T, (width,thick))
    pressure = zeros(T, (width,thick))
    energy = zeros(T, (width,thick))
    # Lattice velocities
    lattice_vel = [0.0 0.0; 1.0 0.0; 0.0 1.0; -1.0 0.0; 0.0 -1.0; 1.0 1.0; -1.0 1.0; -1.0 -1.0; 1.0 -1.0]
    indices = [Symbol("f$(i-1)") for i in 1:9]
    # Store all distributions in an array
    distarray = dist2array(dist)
    # The moments are just a summation of distribution functions
    height = sum(distarray, dims=1)[1, :, :]
    velocityx = sum(distarray .* lattice_vel[:, 1], dims=1)[1, :, :]
    velocityy = sum(distarray .* lattice_vel[:, 2], dims=1)[1, :, :]
    
    velocityx ./= height
    velocityy ./= height

    velocity = JuSwalbe.Twovector(velocityx, velocityy)

    result = JuSwalbe.Macroquant(height, velocity, pressure, energy)
    return result
end

function calculatemoments(dist::JuSwalbe.DistributionD2Q9{Matrix{T}}, force::JuSwalbe.Forces{Twovector{Matrix{T}}}) where {T<:Number}
    # Get a size for the output arrays
    width, thick = size(dist.f0)
    # Allocate the output
    height = zeros(T, (width,thick))
    velocityx = zeros(T, (width,thick))
    velocityy = zeros(T, (width,thick))
    pressure = zeros(T, (width,thick))
    energy = zeros(T, (width,thick))
    # Lattice velocities
    lattice_vel = [0.0 0.0; 1.0 0.0; 0.0 1.0; -1.0 0.0; 0.0 -1.0; 1.0 1.0; -1.0 1.0; -1.0 -1.0; 1.0 -1.0]
    indices = [Symbol("f$(i-1)") for i in 1:9]
    # Store all distributions in an array
    distarray = dist2array(dist)
    forcearrayx, forcearrayy = dist2array(force)
    # The moments are just a summation of distribution functions
    height .= sum(distarray, dims=1)[1, :, :]
    velocityx .= sum(distarray .* lattice_vel[:, 1], dims=1)[1, :, :] .+ T(0.5) * sum(forcearrayx, dims=3)[:, :, 1]
    velocityy .= sum(distarray .* lattice_vel[:, 2], dims=1)[1, :, :] .+ T(0.5) * sum(forcearrayy, dims=3)[:, :, 1]
    
    velocityx ./= height
    velocityy ./= height

    velocity = JuSwalbe.Twovector(velocityx, velocityy)

    result = JuSwalbe.Macroquant(height, velocity, pressure, energy)
    return result
end

"""
    dist2array(dist)

Transforms the `D2Q9` struct into a three dimensional array.

Axis are chosen such that it can be easily multiplied with the nine speed lattice speeds.
The function also excepts a `Forces` struct as input and returns a three dimensional array for both `x` and `y`.

# Example
```jldoctest
julia> using JuSwalbe

julia> a = fill(0.1, (5,5))
5×5 Array{Float64,2}:
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1

julia> b = fill(0.2, (5,5))
5×5 Array{Float64,2}:
 0.2  0.2  0.2  0.2  0.2
 0.2  0.2  0.2  0.2  0.2
 0.2  0.2  0.2  0.2  0.2
 0.2  0.2  0.2  0.2  0.2
 0.2  0.2  0.2  0.2  0.2

julia> dist = JuSwalbe.DistributionD2Q9(f0=a, f1=b, f2=a, f3=a, f4=a, f5=a, f6=a, f7=a, f8=a)
JuSwalbe.DistributionD2Q9{Array{Float64,2}}
  f0: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f1: Array{Float64}((5, 5)) [0.2 0.2 … 0.2 0.2; 0.2 0.2 … 0.2 0.2; … ; 0.2 0.2 … 0.2 0.2; 0.2 0.2 … 0.2 0.2]
  f2: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f3: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f4: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f5: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f6: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f7: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f8: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]

julia> arr = dist2array(dist)
9×5×5 Array{Float64,3}:
[:, :, 1] =
 0.1  0.1  0.1  0.1  0.1
 0.2  0.2  0.2  0.2  0.2
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1

[:, :, 2] =
 0.1  0.1  0.1  0.1  0.1
 0.2  0.2  0.2  0.2  0.2
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1

[:, :, 3] =
 0.1  0.1  0.1  0.1  0.1
 0.2  0.2  0.2  0.2  0.2
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1

[:, :, 4] =
 0.1  0.1  0.1  0.1  0.1
 0.2  0.2  0.2  0.2  0.2
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1

[:, :, 5] =
 0.1  0.1  0.1  0.1  0.1
 0.2  0.2  0.2  0.2  0.2
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1

```

"""
function dist2array(dist::JuSwalbe.DistributionD2Q9{Matrix{T}}) where {T<:Number}
    # Get the size to allocate memory
    width, thick = size(dist.f0)
    # Symbols to access the fields: :f0, :f1, :f2 ...
    indices = [Symbol("f$(i-1)") for i in 1:9]
    # Allocate the output array
    distarray = Array{T, 3}(undef, (9,width,thick))
    for (i, j) in enumerate(indices)
        distarray[i, :, :] = getfield(dist, j)
    end
    return distarray
end

function dist2array(force::JuSwalbe.Forces{JuSwalbe.Twovector{Matrix{T}}}) where {T<:Number}
    # Get the size to allocate memory
    width, thick = size(force.slip.x)
    # Symbols to access the fields:
    indices = [:slip :h∇p :bathymetry :thermal]
    # Allocate the output array
    forcearrayx = Array{T, 3}(undef, (width,thick,4))
    forcearrayy = Array{T, 3}(undef, (width,thick,4))
    for (i, j) in enumerate(indices)
        forcearrayx[:, :, i] = getfield(getfield(force, j),:x)
        forcearrayy[:, :, i] = getfield(getfield(force, j),:y)
    end
    return forcearrayx,forcearrayy
end