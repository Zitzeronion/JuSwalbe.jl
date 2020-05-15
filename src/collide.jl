"""
    collisionBGK(mom, tempdist, input)

Computes the BGK collision operation for a `D1Q3` and a `D2Q9` lattice Boltzmann algorithm.

The collision operator is the central object in kinetic theory. 
However there is no general solution for this operator. 
One work around is to make assumptions and try to model them. 
The BGK operator is one of the fairly successful operators. 
It stems from the idea that the distributions are close to their equilibria and that all of them *relax* with a single relaxation time `τ`.

# Math
The collision operator is usually the right hand side of the kinetic equation 

`` \\partial_t f + \\frac{\\mathbf{p}}{m}\\cdot\\nabla f + \\mathbf{F}\\cdot \\partial_{\\mathbf{p}} f = \\Omega . ``

After tedious mathematics and expansions one can obtain the simplest collision operator the BGK as

`` \\Omega_{BGK} = -\\frac{f - f^{eq}}{\\tau}, ``

with `τ` being the relaxation parameter. 

# Example

## One spatial dimension
```jldoctest
julia> using JuSwalbe

julia> 


```

# References

"""
function collisionBGK(mom::JuSwalbe.Macroquant{Vector{T},Vector{T}}, forces::JuSwalbe.Forces{Vector{T}}, tempdist::JuSwalbe.DistributionD1Q3{Vector{T}}, input::JuSwalbe.Inputconstants) where {T<:Number}
    # Get the size, type and allocate result array
    len = length(mom.height)
    weights = [T(2/3) T(1/6) T(1/6)]
    clat = [T(0) T(1) T(-1)]
    τ = input.τ
    gravity = input.gravity  
    csquared = T(3)
    allforces = zeros(T, len)
    allforces .= sum(hcat(forces.slip, forces.h∇p, forces.thermal, forces.bathymetry), dims=2)[:,1]

    newdist = JuSwalbe.DistributionD1Q3(f0=zeros(T, len),
                                        f1=zeros(T, len),
                                        f2=zeros(T, len))

    eqdist = calc_equilibrium_distribution(mom; gravity=T(gravity))

    newdist.f0 .= T(1-1/τ) * tempdist.f0 .+ T(1/τ) * eqdist.f0 
    newdist.f1 .= T(1-1/τ) * tempdist.f1 .+ T(1/τ) * eqdist.f1 .+ weights[2] * csquared * clat[2] * allforces
    newdist.f2 .= T(1-1/τ) * tempdist.f2 .+ T(1/τ) * eqdist.f2 .+ weights[3] * csquared * clat[3] * allforces

    return newdist
end

"""
    streamdistperiodic!(dist)

Streams the lattice Boltzmann distribution accoriding to their velocity vector.
"""
function streamdistperiodic!(dist::JuSwalbe.DistributionD1Q3{Vector{T}}) where {T<:Number}
    len = length(dist.f0)
    f1dummy = zeros(T, len)
    f2dummy = zeros(T, len)
    f1dummy .= circshift(dist.f1, 1)
    f2dummy .= circshift(dist.f2, -1)
end
