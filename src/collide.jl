"""
    collisionBGK(mom, tempdist, input)

Computes the collision operation for a `D1Q3` and a `D2Q9` lattice Boltzmann algorithm.

# Math

# Example

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
