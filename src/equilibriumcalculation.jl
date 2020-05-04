"""
    calc_equilibrium_distribution(mom:macroquant64_1d, gravity=0)

Calculates the equilibrium distributions based on the macroscopic quantities `mom` and gravitational acceleration `gravity`.

The equilibrium distribtutions are at the heart of the lattice Boltzmann method.
As the expansion is made around the equilibrium, lattice Boltzmann always assumes that the flow field is close to equilibrium.
Therefore the equilibrium is calculated from macroscopic quantities, i.e. height `.height` and velocity `.velocity`.
There are plenty of ways to calculate them I mainly used: [Pham van Thang, Bastien Chopard, Laurent Lefevre, Diemer Anda Ondo, Eduardo Mendes. 
Study of the 1D lattice Boltzmann shallow water equation and its coupling to build a canal network. 
Journal of Computational Physics, Elsevier, 2010, 229 (19), pp.7373-7400. 10.1016/j.jcp.2010.06.022](https://www.sciencedirect.com/science/article/pii/S0021999110003372)
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