"""
    collisionBGK(mom, tempdist, input)

Computes the collision operation for a `D1Q3` and a `D2Q9` lattice Boltzmann algorithm.

# Math

# Example

# References

"""
function collisionBGK(mom::JuSwalbe.Macroquant{Vector{T},Vector{T}}, forces::JuSwalbe.Forces{Vector{T}}, tempdist::JuSwalbe.DistributionD1Q3{Vector{T}}, input::JuSwalbe.Inputconstants) where {T<:Number}
    # Get the size, type and allocate result array
    width, thick = size(mom.height)
    newdist = zeros(T, (9,width, thick))
    


end

function streamdist(dist::JuSwalbe.DistributionD1Q3{Vector{T}}) where {T<:Number}
    
end