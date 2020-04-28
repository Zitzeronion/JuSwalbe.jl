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

mutable struct macroquant64_2d <: macroscopic_quantity
    height::Matrix{Float64}
    velocity::Matrix{Float64}
    energy::Matrix{Float64}
end

mutable struct macroquant32_2d <: macroscopic_quantity
    height::Matrix{Float32}
    velocity::Matrix{Float32}
    energy::Matrix{Float32}
end