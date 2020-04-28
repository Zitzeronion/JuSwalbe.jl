abstract type dist end

mutable struct dist64_1d <: dist
    f0::Vector{Float64}
    f1::Vector{Float64}
    f2::Vector{Float64}
end

mutable struct dist32_1d <: dist
    f0::Vector{Float32}
    f1::Vector{Float32}
    f2::Vector{Float32}
end

mutable struct dist64_2d <: dist
    f0::Matrix{Float64}
    f1::Matrix{Float64}
    f2::Matrix{Float64}
    f3::Matrix{Float64}
    f4::Matrix{Float64}
    f5::Matrix{Float64}
    f6::Matrix{Float64}
    f7::Matrix{Float64}
    f8::Matrix{Float64}
end

mutable struct dist32_2d <: dist
    f0::Matrix{Float32}
    f1::Matrix{Float32}
    f2::Matrix{Float32}
    f3::Matrix{Float32}
    f4::Matrix{Float32}
    f5::Matrix{Float32}
    f6::Matrix{Float32}
    f7::Matrix{Float32}
    f8::Matrix{Float32}
end