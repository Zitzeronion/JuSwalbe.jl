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

julia> input, mom, force, dist = minimalsetup1d(10)
(JuSwalbe.Inputconstants
  lx: Int64 10
  ly: Int64 512
  maxruntime: Int64 100000
  dumping: Int64 1000
  τ: Float64 1.0
  gravity: Float64 0.0
  γ: Float64 0.01
  δ: Float64 1.0
  μ: Float64 0.16666666666666666
, JuSwalbe.Macroquant{Array{Float64,1},Array{Float64,1}}
  height: Array{Float64}((10,)) [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  velocity: Array{Float64}((10,)) [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
  pressure: Array{Float64}((10,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  energy: Array{Float64}((10,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
, JuSwalbe.Forces{Array{Float64,1}}
  slip: Array{Float64}((10,)) [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
  h∇p: Array{Float64}((10,)) [-0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1]
  bathymetry: Array{Float64}((10,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  thermal: Array{Float64}((10,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
, JuSwalbe.DistributionD1Q3{Array{Float64,1}}
  f0: Array{Float64}((10,)) [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  f1: Array{Float64}((10,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  f2: Array{Float64}((10,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
)

julia> force.slip = zeros(10)
10-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0

julia> force.h∇p = zeros(10)
10-element Array{Float64,1}:
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0
 0.0

julia> newdist = collisionBGK(mom, force, dist, input)
JuSwalbe.DistributionD1Q3{Array{Float64,1}}
  f0: Array{Float64}((10,)) [0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99]
  f1: Array{Float64}((10,)) [0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001]
  f2: Array{Float64}((10,)) [-0.045, -0.045, -0.045, -0.045, -0.045, -0.045, -0.045, -0.045, -0.045, -0.045]

julia> equi = calc_equilibrium_distribution(mom, gravity=input.gravity)
JuSwalbe.DistributionD1Q3{Array{Float64,1}}
  f0: Array{Float64}((10,)) [0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99, 0.99]
  f1: Array{Float64}((10,)) [0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001, 0.05500000000000001]
  f2: Array{Float64}((10,)) [-0.045, -0.045, -0.045, -0.045, -0.045, -0.045, -0.045, -0.045, -0.045, -0.045]

julia> using Test; @test equi.f0 == newdist.f0 # without forces and τ = 1 they have to be equal
Test Passed
```

# References
## One spatial dimension
- [Asymmetric lattice Boltzmann model for shallow water flows](https://www.sciencedirect.com/science/article/abs/pii/S0045793013003599)

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

function collisionBGK(mom::JuSwalbe.Macroquant{Matrix{T},JuSwalbe.Twovector{Matrix{T}}}, forces::JuSwalbe.Forces{JuSwalbe.Twovector{Matrix{T}}}, tempdist::JuSwalbe.DistributionD2Q9{Matrix{T}}, input::JuSwalbe.Inputconstants) where {T<:Number}
  # Get the size, type and allocate result array
  width, thick = size(mom.height)
  weights = [T(4/9) T(1/9) T(1/9) T(1/9) T(1/9) T(1/36) T(1/36) T(1/36) T(1/36)]
  clat = [T(0) T(0); T(1) T(0); T(0) T(1); T(-1) T(0); T(0) T(-1); T(1) T(1); T(-1) T(1); T(-1) T(-1); T(1) T(-1)]
  clatsquare = [T(0) T(1) T(1) T(1) T(1) T(2) T(2) T(2) T(2)]
  τ = input.τ
  gravity = input.gravity  
  csquared = T(3)
  allforces = dist2array(forces)

  newdist = JuSwalbe.DistributionD2Q9(f0=zeros(T, (width, thick)),
                                      f1=zeros(T, (width, thick)),
                                      f2=zeros(T, (width, thick)),
                                      f3=zeros(T, (width, thick)),
                                      f4=zeros(T, (width, thick)),
                                      f5=zeros(T, (width, thick)),
                                      f6=zeros(T, (width, thick)),
                                      f7=zeros(T, (width, thick)),
                                      f8=zeros(T, (width, thick)))

  eqdist = calc_equilibrium_distribution(mom; gravity=T(gravity))
  
  newdist.f0 .= T(1-1/τ) * tempdist.f0 .+ T(1/τ) * eqdist.f0 

  newdist.f1 .= (T(1-1/τ) * tempdist.f1 .+ T(1/τ) * eqdist.f1 
              .+ (weights[2] * csquared)/clatsquare[2] * (clat[2,1] * sum(allforces.x, dims=3)[:, :, 1] .+ clat[2,2] * sum(allforces.y, dims=3)[:, :, 1]))
  newdist.f2 .= (T(1-1/τ) * tempdist.f2 .+ T(1/τ) * eqdist.f2 
              .+ (weights[3] * csquared)/clatsquare[3] * (clat[3,1] * sum(allforces.x, dims=3)[:, :, 1] .+ clat[3,2] * sum(allforces.y, dims=3)[:, :, 1]))
  newdist.f3 .= (T(1-1/τ) * tempdist.f3 .+ T(1/τ) * eqdist.f3 
              .+ (weights[4] * csquared)/clatsquare[4] * (clat[4,1] * sum(allforces.x, dims=3)[:, :, 1] .+ clat[4,2] * sum(allforces.y, dims=3)[:, :, 1]))
  newdist.f4 .= (T(1-1/τ) * tempdist.f4 .+ T(1/τ) * eqdist.f4 
              .+ (weights[5] * csquared)/clatsquare[5] * (clat[5,1] * sum(allforces.x, dims=3)[:, :, 1] .+ clat[5,2] * sum(allforces.y, dims=3)[:, :, 1]))
  newdist.f5 .= (T(1-1/τ) * tempdist.f5 .+ T(1/τ) * eqdist.f5 
              .+ (weights[6] * csquared)/clatsquare[6] * (clat[6,1] * sum(allforces.x, dims=3)[:, :, 1] .+ clat[6,2] * sum(allforces.y, dims=3)[:, :, 1]))
  newdist.f6 .= (T(1-1/τ) * tempdist.f6 .+ T(1/τ) * eqdist.f6 
              .+ (weights[7] * csquared)/clatsquare[7] * (clat[7,1] * sum(allforces.x, dims=3)[:, :, 1] .+ clat[7,2] * sum(allforces.y, dims=3)[:, :, 1]))
  newdist.f7 .= (T(1-1/τ) * tempdist.f7 .+ T(1/τ) * eqdist.f7 
              .+ (weights[8] * csquared)/clatsquare[8] * (clat[8,1] * sum(allforces.x, dims=3)[:, :, 1] .+ clat[8,2] * sum(allforces.y, dims=3)[:, :, 1]))
  newdist.f8 .= (T(1-1/τ) * tempdist.f8 .+ T(1/τ) * eqdist.f8 
              .+ (weights[9] * csquared)/clatsquare[9] * (clat[9,1] * sum(allforces.x, dims=3)[:, :, 1] .+ clat[9,2] * sum(allforces.y, dims=3)[:, :, 1]))
  
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

    dist.f1 .= f1dummy
    dist.f2 .= f2dummy
    
    return dist
end
