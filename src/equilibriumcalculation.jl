function calc_equilibrium_distribution(mom::macroquant64_1d, gravity=0)
    len = length(mom.height)
    v = 1.0
    vsquare = v^2
    equi_dist = dist64_1d(zeros(Float64, len), zeros(Float64, len), zeros(Float64, len))
    for i in 1:len
        equi_dist.f0[i] = mom.height[i] 
                        - 1/(2*vsquare) * gravity * mom.height[i]^2 
                        - 1/vsquare * mom.height[i] * mom.velocity[i]^2
        equi_dist.f1[i] = 1/(4*vsquare) * gravity * mom.height[i]^2 
                        + 1/(2*v) * mom.height[i] * mom.velocity[i] 
                        + 1/(2*vsquare) * mom.height[i] * mom.velocity[i]^2
        equi_dist.f2[i] = 1/(4*vsquare) * gravity * mom.height[i]^2 
                        - 1/(2*v) * mom.height[i] * mom.velocity[i] 
                        + 1/(2*vsquare) * mom.height[i] * mom.velocity[i]^2
    end
    return equi_dist
end

function calc_equilibrium_distribution(mom::macroquant32_1d, gravity=0)
    len = length(mom.height)
    v = 1.0
    vsquare = v^2
    equi_dist = dist32_1d(zeros(Float32, len), zeros(Float32, len), zeros(Float32, len))

    for i in 1:len
        equi_dist.f0[i] = mom.height[i] - 1/(2*vsquare) * gravity * mom.height[i]^2 
                        - 1/vsquare * mom.height[i] * mom.velocity[i]^2
        equi_dist.f1[i] = 1/(4*vsquare) * gravity * mom.height[i]^2 
                        + 1/(2*v) * mom.height[i] * mom.velocity[i] 
                        + 1/(2*vsquare) * mom.height[i] * mom.velocity[i]^2
        equi_dist.f2[i] = 1/(4*vsquare) * gravity * mom.height[i]^2 
                        - 1/(2*v) * mom.height[i] * mom.velocity[i] 
                        + 1/(2*vsquare) * mom.height[i] * mom.velocity[i]^2
    end
    return equi_dist
end