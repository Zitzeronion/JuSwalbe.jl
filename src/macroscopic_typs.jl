abstract type Macroscopic_quantity end

@with_kw mutable struct Macroquant{T, K} <: Macroscopic_quantity
    # In one spatial dimension T::Vector{Float32,64} <: Array{Float32,64, 1}
    # In two spatial dimension T::Matrix{Float32,64} <: Array{Float32.64, 2}
    # The macro with_kw adds a specific constructor operation to the struct.
    height::T
    velocity::K
    pressure::T
    energy::T
end

@with_kw mutable struct Forces{T} <: Macroscopic_quantity
    slip::T
    h∇p::T
    bathymetry::T
    thermal::T
end

@with_kw mutable struct Twovector{T} <: Macroscopic_quantity
    x::T
    y::T
end

@with_kw mutable struct GPUvec <: Macroscopic_quantity
    x::CuArray{Float32,2,Nothing}
    y::CuArray{Float32,2,Nothing}
end

@with_kw mutable struct MacroquantGPU <: Macroscopic_quantity
    height::CuArray{Float32,2,Nothing}
    velocity::GPUvec
    pressure::CuArray{Float32,2,Nothing}
    energy::CuArray{Float32,2,Nothing}
end

@with_kw mutable struct ForcesGPU <: Macroscopic_quantity
    slip::GPUvec
    h∇p::GPUvec
    bathymetry::GPUvec
    thermal::GPUvec
end
