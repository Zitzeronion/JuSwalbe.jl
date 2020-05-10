"""
    calc_equilibrium_distribution(mom::JuSwalbe.Macroquant; gravity=0)

Calculates the equilibrium distributions based on the macroscopic quantities `mom` and gravitational acceleration `gravity`.

The equilibrium distribtutions are at the heart of the lattice Boltzmann method.
As the expansion is made around the equilibrium, lattice Boltzmann always assumes that the flow field is close to equilibrium.
Therefore the equilibrium is calculated from macroscopic quantities, i.e. height `.height` and velocity `.velocity`.

# Math:
## One spatial dimension
Defining equations in one spatial dimensions are expressed in the following way

``f_0^{eq} =  h - \\frac{1}{2v^2}gh^2 - \\frac{1}{v^2}hu^2`` 

``f_1^{eq} =  \\frac{1}{4v^2}gh^2 + \\frac{1}{2v}hu + \\frac{1}{2v^2}hu^2``

``f_2^{eq} =  \\frac{1}{4v^2}gh^2 - \\frac{1}{2v}hu + \\frac{1}{2v^2}hu^2 ``

## Two spatial dimension
In two spatial dimensions the velocity because a vector with a `x` and `y` component.
Such the equations look a little more complex

``f_0^{eq} =  h - \\frac{4}{9}h(\\frac{15}{2}gh - \\frac{3}{2}u^2)``

``f_i^{eq} = w_i h(\\frac{3}{2}gh + 3 \\mathbf{c}_i\\cdot\\mathbf{u} + \\frac{9}{2}(\\mathbf{c}_i\\cdot\\mathbf{u})^2 - \\frac{3}{2}u^2))``

with i are the number of lattice speeds and w_i and c_i are the weights and sets of lattice velocities. 

# Example
TBD

# References
## One spatial dimension
There are plenty of ways to calculate them I mainly the ones derived in Eq.(13) of [Study of the 1D lattice Boltzmann shallow water equation and its coupling to build a canal network.](https://www.sciencedirect.com/science/article/pii/S0021999110003372)

## Two spatial dimensions
Here I stick to the paper written by Paul Dellar, Eq.(26) of [Nonhydrodynamic modes and *a priori* construction of shallow water lattice Boltzmann equations](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.65.036309)
Originally these equilibria have been worked out by Paul Salmon.
"""
function calc_equilibrium_distribution(mom::JuSwalbe.Macroquant{Vector{T}, Vector{T}}; gravity::T=T(0)) where {T<:Number}
    len = length(mom.height)
    v = 1.0
    vsquare = v^2
    equi_dist = JuSwalbe.DistributionD1Q3(zeros(T, len), zeros(T, len), zeros(T, len))
    equi_dist.f0 = (mom.height .- 1/(2*vsquare) * gravity .* mom.height.^2 
                    .- 1/vsquare .* mom.height .* mom.velocity.^2)
    equi_dist.f1 = (1/(4*vsquare) * gravity .* mom.height.^2 
                    .+ 1/(2*v) .* mom.height .* mom.velocity 
                    .+ 1/(2*vsquare) .* mom.height .* mom.velocity.^2)
    equi_dist.f2 = (1/(4*vsquare) * gravity .* mom.height.^2 
                    .- 1/(2*v) .* mom.height .* mom.velocity 
                    .+ 1/(2*vsquare) .* mom.height .* mom.velocity.^2)
        
    return equi_dist
end

function calc_equilibrium_distribution(mom::JuSwalbe.Macroquant{Matrix{T}, JuSwalbe.Twovector{Matrix{T}}}; gravity::T=T(0)) where {T<:Number}
    # In two dimensions there is a little more work needed
    width = size(mom.height,1)
    thick = size(mom.height,2)

    v = 1.0
    vsquare = v^2
    # Standard D2Q9 weights and lattice velocities
    # For type security, otherwise weights will always be Float64
    if T != Float64
        weights = Array{T, 1}(undef, 9)
        w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
        for (ind, ele) in enumerate(weights)
            weights[ind] = w[ind]
        end
    else
        weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
    end
    lattice_vel = [0 0; 1 0; 0 1; -1 0; 0 -1; 1 1; -1 1; -1 -1; 1 -1]
    # I use dicts as they are easier to iterate over
    weightsdict = Dict()
    latticeveldict = Dict()
    for (dist, weight) in enumerate(weights)
        weightsdict[Symbol("f$(dist-1)")] = weight
        latticeveldict[Symbol("f$(dist-1)")]=lattice_vel[dist,:]
    end
    # f0 has a special formular to be calculated
    relevantkeys = [:f1, :f2, :f3, :f4, :f5, :f6, :f7, :f8]
    # Initialize a dummy distribution
    equi_dist = JuSwalbe.DistributionD2Q9(zeros(T,(width, thick)), 
                                          zeros(T,(width, thick)), 
                                          zeros(T,(width, thick)),
                                          zeros(T,(width, thick)),
                                          zeros(T,(width, thick)),
                                          zeros(T,(width, thick)),
                                          zeros(T,(width, thick)),
                                          zeros(T,(width, thick)),
                                          zeros(T,(width, thick)))
    # Calculate the velocity squared
    udotu = velocitysquared(mom)
    # Again make use dicts for iteration purpose
    equidict = struct2dict(equi_dist)
    # Calculate f0
    equidict[:f0] = mom.height .- weights[1] * mom.height .*(15.0/2.0 * gravity * mom.height .- 3.0/2.0 * udotu)
    # And the other 8 distribtuion functions
    for key in relevantkeys
        latticeveldotu = mom.velocity.x * latticeveldict[key][1] .+ mom.velocity.y * latticeveldict[key][2]
        equidict[key] = weightsdict[key] * mom.height .* (3.0/2.0 * gravity * mom.height
                                                          .+ 3 * latticeveldotu 
                                                          .+ 9.0/2.0 * latticeveldotu.^2 
                                                          .- 3.0/2.0 * udotu)
    end
    # Save them to the distrubtion type which was created earlier
    equi_dist = JuSwalbe.DistributionD2Q9(equidict[:f0],
                                          equidict[:f1],
                                          equidict[:f2],
                                          equidict[:f3],
                                          equidict[:f4],
                                          equidict[:f5],
                                          equidict[:f6],
                                          equidict[:f7],
                                          equidict[:f8])

    # And return it :)
    return equi_dist
end

"""
    velocitysquared(mom::Macroquant{Matrix{T}, JuSwalbe.Twovector{Matrix{T}}})

Computes the square of the velocity vector (ux, uy) at every lattice point.

The magnitude of the velocity is needed to calculate the equilibrium distribution in the dimensional case.
In the one dimensional case the velocity is just a vector and therefore has no `.x` and `.y` component.
See also [`calc_equilibrium_distribution`](@ref)

# Math
The velocity squared ``u^2(x,y)`` is computed according to 

`` u^2(x,y) = (u_x, u_y)^2(x,y) = u_x^2(x,y) + u_y^2(x,y).``

With lower case `x` and `y` the respective component of the velocity vector is addressed.

# Example
```jldoctest
julia> using JuSwalbe

julia> velocities = JuSwalbe.Twovector{Matrix{Float64}}(fill(0.1, (4,4)),fill(0.2, (4,4)))
JuSwalbe.Twovector{Array{Float64,2}}([0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1], [0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2])

julia> moment = JuSwalbe.Macroquant{Matrix{Float64},JuSwalbe.Twovector{Matrix{Float64}}}(ones(Float64, (4,4)), velocities, ones(Float64, (4,4)), ones(Float64, (4,4)))
JuSwalbe.Macroquant{Array{Float64,2},JuSwalbe.Twovector{Array{Float64,2}}}([1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0], JuSwalbe.Twovector{Array{Float64,2}}([0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1], [0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2]), [1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0], [1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0])

julia> moment.velocity.x
4×4 Array{Float64,2}:
 0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1

julia> JuSwalbe.velocitysquared(moment)
4×4 Array{Float64,2}:
 0.05  0.05  0.05  0.05
 0.05  0.05  0.05  0.05
 0.05  0.05  0.05  0.05
 0.05  0.05  0.05  0.05
```
"""
function velocitysquared(mom::Macroquant{Matrix{T}, JuSwalbe.Twovector{Matrix{T}}}) where {T<:Number}
    velsquared = mom.velocity.x.^2 .+ mom.velocity.y.^2
    return velsquared
end

