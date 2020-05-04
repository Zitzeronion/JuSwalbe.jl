"""
    calc_equilibrium_distribution(mom:macroquant64_1d, gravity=0)

Calculates the equilibrium distributions based on the macroscopic quantities `mom` and gravitational acceleration `gravity`.

The equilibrium distribtutions are at the heart of the lattice Boltzmann method.
As the expansion is made around the equilibrium, lattice Boltzmann always assumes that the flow field is close to equilibrium.
Therefore the equilibrium is calculated from macroscopic quantities, i.e. height `.height` and velocity `.velocity`.
There are plenty of ways to calculate them I mainly the ones derived in Eq.(13) of [Study of the 1D lattice Boltzmann shallow water equation and its coupling to build a canal network.](https://www.sciencedirect.com/science/article/pii/S0021999110003372)

# Math:
Defining equations in one spatial dimensions are expressed in the following way
``\\begin{align*} 
f_0^{eq} &=  h - \\frac{1}{2v^2}gh^2 - \\frac{1}{v^2}hu^2 \\ 
f_1^{eq} &=  \\frac{1}{4v^2}gh^2 + \\frac{1}{2v}hu + \\frac{1}{2v^2}hu^2 \\
f_2^{eq} &=  \\frac{1}{4v^2}gh^2 - \\frac{1}{2v}hu + \\frac{1}{2v^2}hu^2 
\\end{align*}``

# Example
TBD
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
    width = size(mom.height,1)
    thick = size(mom.height,2)

    v = 1.0
    vsquare = v^2

    weights = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36]
    lattice_vel = [0 0; 1 0; 0 1; -1 0; 0 -1; 1 1; -1 1; -1 -1; 1 -1]

    equi_dist = dist64_2d(zeros(Float64,(width, thick)), 
                          zeros(Float64,(width, thick)), 
                          zeros(Float64,(width, thick)),
                          zeros(Float64,(width, thick)),
                          zeros(Float64,(width, thick)),
                          zeros(Float64,(width, thick)),
                          zeros(Float64,(width, thick)),
                          zeros(Float64,(width, thick)),
                          zeros(Float64,(width, thick)))
    
    udotu = velocitysquared(mom)
    equidict = struct2dict(equi_dist)

    equi_dist.f0 = mom.height .- weights[1] .* mom.height .*(15.0/2.0 * gravity .* mom.height - 3.0/2.0 .* udotu)

end

function velocitysquared(mom::macroquant64_2d)
    velsquared = mom.velocity.x.^2 + mom.velocity.y.^2
    return velsquared
end

