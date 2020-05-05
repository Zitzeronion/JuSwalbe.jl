"""
    calc_equilibrium_distribution(mom:macroquant64_1d, gravity=0)

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
``f_i^{eq} = w_i h(\\frac{3}{2}gh + 3 \\mathbf{c}_i\\cdot\\mathbf{u} + \\frac{9}{2}(\\mathbf{c}_i\\cdot\\mathbf{u})^2 - \\frac{3}{2}u^2))

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
function calc_equilibrium_distribution(mom::macroquant64_1d, gravity=0)
    len = length(mom.height)
    v = 1.0
    vsquare = v^2
    equi_dist = dist64_1d(zeros(Float64, len), zeros(Float64, len), zeros(Float64, len))
    for i in 1:len
        equi_dist.f0[i] = (mom.height[i] 
                        - 1/(2*vsquare) * gravity * mom.height[i]^2 
                        - 1/vsquare * mom.height[i] * mom.velocity[i]^2)
        equi_dist.f1[i] = (1/(4*vsquare) * gravity * mom.height[i]^2 
                        + 1/(2*v) * mom.height[i] * mom.velocity[i] 
                        + 1/(2*vsquare) * mom.height[i] * mom.velocity[i]^2)
        equi_dist.f2[i] = (1/(4*vsquare) * gravity * mom.height[i]^2 
                        - 1/(2*v) * mom.height[i] * mom.velocity[i] 
                        + 1/(2*vsquare) * mom.height[i] * mom.velocity[i]^2)
    end
    return equi_dist
end

function calc_equilibrium_distribution(mom::macroquant32_1d, gravity=0)
    len = length(mom.height)
    v = 1.0
    vsquare = v^2
    equi_dist = dist32_1d(zeros(Float32, len), zeros(Float32, len), zeros(Float32, len))

    for i in 1:len
        equi_dist.f0[i] = (mom.height[i] - 1/(2*vsquare) * gravity * mom.height[i]^2 
                        - 1/vsquare * mom.height[i] * mom.velocity[i]^2)
        equi_dist.f1[i] = (1/(4*vsquare) * gravity * mom.height[i]^2 
                        + 1/(2*v) * mom.height[i] * mom.velocity[i] 
                        + 1/(2*vsquare) * mom.height[i] * mom.velocity[i]^2)
        equi_dist.f2[i] = (1/(4*vsquare) * gravity * mom.height[i]^2 
                        - 1/(2*v) * mom.height[i] * mom.velocity[i] 
                        + 1/(2*vsquare) * mom.height[i] * mom.velocity[i]^2)
    end
    return equi_dist
end

function calc_equilibrium_distribution(mom::macroquant64_2d, gravity=0)
    # In two dimensions there is a little more work needed
    width = size(mom.height,1)
    thick = size(mom.height,2)

    v = 1.0
    vsquare = v^2
    # Standard D2Q9 weights and lattice velocities
    weights = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]
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
    equi_dist = dist64_2d(zeros(Float64,(width, thick)), 
                          zeros(Float64,(width, thick)), 
                          zeros(Float64,(width, thick)),
                          zeros(Float64,(width, thick)),
                          zeros(Float64,(width, thick)),
                          zeros(Float64,(width, thick)),
                          zeros(Float64,(width, thick)),
                          zeros(Float64,(width, thick)),
                          zeros(Float64,(width, thick)))
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
    equi_dist = dist64_2d(equidict[:f0],
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
    velocitysquared(mom::macroquant64_2d)

Computes the square of the velocity vector (ux, uy) at every lattice point.

The magnitude of the velocity is needed to calculate the equilibrium distribution.

# Math
The velocity squared ``u^2(x,y)`` is computed according to 

`` u^2(x,y) = (u_x, u_y)^2(x,y) = u_x^2(x,y) + u_y^2(x,y).``

With lower case `x` and `y` the respective component of the velocity vector is addressed.

# Example
```jldoctest
julia> using JuSwalbe

julia> moment = JuSwalbe.macroquant64_2d(ones(4,4), JuSwalbe.velocity64_2d(fill(0.1, (4,4)),fill(0.2, (4,4))), zeros(4,4))
JuSwalbe.macroquant64_2d([1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0], JuSwalbe.velocity64_2d([0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1; 0.1 0.1 0.1 0.1], [0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2; 0.2 0.2 0.2 0.2]), [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0])

julia> moment.velocity.x
4×4 Array{Float64,2}:
 0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1

julia> velocitysquared(moment)
4×4 Array{Float64,2}:
 0.05  0.05  0.05  0.05
 0.05  0.05  0.05  0.05
 0.05  0.05  0.05  0.05
 0.05  0.05  0.05  0.05
"""
function velocitysquared(mom::macroquant64_2d)
    velsquared = mom.velocity.x.^2 + mom.velocity.y.^2
    return velsquared
end

