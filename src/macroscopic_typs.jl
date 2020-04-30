abstract type macroscopic_quantity end

mutable struct macroquant64_1d <: macroscopic_quantity
    height::Vector{Float64}
    velocity::Vector{Float64}
    energy::Vector{Float64}
end

mutable struct macroquant32_1d <: macroscopic_quantity
    height::Vector{Float32}
    velocity::Vector{Float32}
    energy::Vector{Float32}
end
# TODO: Fix energy type, as it should be a tensor...
mutable struct velocity64_2d <: macroscopic_quantity
    x::Matrix{Float64}
    y::Matrix{Float64}
end

mutable struct velocity32_2d <: macroscopic_quantity
    x::Matrix{Float32}
    y::Matrix{Float32}
end

mutable struct macroquant64_2d <: macroscopic_quantity
    height::Matrix{Float64}
    velocity::velocity64_2d
    energy::Matrix{Float64}
end

mutable struct macroquant32_2d <: macroscopic_quantity
    height::Matrix{Float32}
    velocity::velocity32_2d
    energy::Matrix{Float32}
end