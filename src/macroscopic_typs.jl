abstract type Macroscopic_quantity end

mutable struct Macroquant{T, K} <: Macroscopic_quantity
    # In one spatial dimension T::Vector{Float32,64} <: Array{Float32, 1}
    # In two spatial dimension T::Matrix{Float32,64} <: Array{Float32, 2}
    height::T
    velocity::K
    pressure::T
    energy::T
end

mutable struct Forces{T} <: Macroscopic_quantity
    slip::T
    ∇p::T
    bathymetry::T
    thermal::T
end

mutable struct Twovector{T} <: Macroscopic_quantity
    x::T
    y::T
end

