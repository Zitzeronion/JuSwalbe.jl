#============================================================================#
#                               Structures                                   #
#============================================================================#
@with_kw struct params{T}
    lx::T = 512
    ly::T = 512
    mt::T = 20000
end

@with_kw mutable struct velocity{T}
    para::params 
    x::Array{T,2} = zeros(T, (para.lx, para.ly)) 
    y::Array{T,2} = zeros(T, (para.lx, para.ly))
end

@with_kw mutable struct force{T}
    para::params 
    x::Array{T,2} = zeros(T, (para.lx, para.ly)) 
    y::Array{T,2} = zeros(T, (para.lx, para.ly))
end

@with_kw mutable struct cuforce{T}
    para::params 
    x = CUDA.zeros(T, (para.lx, para.ly)) 
    y = CUDA.zeros(T, (para.lx, para.ly))
end

@with_kw mutable struct moments{T}
    para::params
    height::Array{T,2} = ones(T, (para.lx, para.ly))
    pressure::Array{T,2} = zeros(T, (para.lx, para.ly))
end