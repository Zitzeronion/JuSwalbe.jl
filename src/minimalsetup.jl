"""
    minimalsetup1d(N)

Creates a minimal setup of all important variables.

Such it generates the moments `height`, `velocity` as well as a dummy distribution with its speeds `f0`, `f1` and `f2`.
Further a set of forcing is supplied, they can be used with `foces.slip`, `.h∇p`, `.bathymetry` and `.thermal`.

# Example
```jldoctest
julia> using JuSwalbe

julia> constants, mom, f, distribution = minimalsetup1d(10)
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

julia> mom.height
10-element Array{Float64,1}:
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
 1.0
```
See also: [`minimalsetup2d`](@ref)
"""
function minimalsetup1d(N::Int)
    input = Inputconstants(lx=N)
    mom = Macroquant(height = ones(N), velocity=fill(0.1, N), pressure=zeros(N), energy=zeros(N))
    forces = Forces(slip=fill(0.1, N), thermal=zeros(N), h∇p=fill(-0.1, N), bathymetry=zeros(N))
    dist = DistributionD1Q3(f0 = ones(N), f1 = zeros(N), f2 = zeros(N))
    return input, mom, forces, dist
end

"""
    minimalsetup2d(N,M)

Creates a minimal working example set of macroscopic moments, forces and a distribution.

# Example
```jldoctest
julia> using JuSwalbe

julia> input, mom, f, dist = minimalsetup2d(5,5)
(JuSwalbe.Inputconstants
  lx: Int64 5
  ly: Int64 5
  maxruntime: Int64 100000
  dumping: Int64 1000
  τ: Float64 1.0
  gravity: Float64 0.0
  γ: Float64 0.01
  δ: Float64 1.0
  μ: Float64 0.16666666666666666
, JuSwalbe.Macroquant{Array{Float64,2},JuSwalbe.Twovector{Array{Float64,2}}}
  height: Array{Float64}((5, 5)) [1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0; … ; 1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0]
  velocity: JuSwalbe.Twovector{Array{Float64,2}}
  pressure: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  energy: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
, JuSwalbe.Forces{JuSwalbe.Twovector{Array{Float64,2}}}
  slip: JuSwalbe.Twovector{Array{Float64,2}}
  h∇p: JuSwalbe.Twovector{Array{Float64,2}}
  bathymetry: JuSwalbe.Twovector{Array{Float64,2}}
  thermal: JuSwalbe.Twovector{Array{Float64,2}}
, JuSwalbe.DistributionD2Q9{Array{Float64,2}}
  f0: Array{Float64}((5, 5)) [1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0; … ; 1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0]
  f1: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f2: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f3: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f4: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f5: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f6: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f7: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f8: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
)

```

See also: [`minimalsetup1d`](@ref)
"""
function minimalsetup2d(N::Int, M::Int)
    input = Inputconstants(lx=N, ly=M)
    vel = Twovector(x=fill(0.02, (N,M)), y=fill(-0.02, (N,M)))
    zerovec = Twovector(x=zeros(N,M), y=zeros(N,M))
    mom = Macroquant(height = ones(N,M), velocity=vel, pressure=zeros(N,M), energy=zeros(N,M))
    forces = Forces(slip=vel, thermal=zerovec, h∇p=vel, bathymetry=zerovec)
    dist = DistributionD2Q9(f0 = ones(N,M), f1 = zeros(N,M), f2 = zeros(N,M), 
                            f3 = zeros(N,M), f4 = zeros(N,M), f5 = zeros(N,M), 
                            f6 = zeros(N,M), f7 = zeros(N,M), f8 = zeros(N,M))
    return input, mom, forces, dist
end

"""
  simplemoment2d(n,m,type)

Creates moments of type JuSwalbe.Macroquant with dimensions (n,m).

# Example
```jldoctest
julia> using JuSwalbe

julia> mom = simplemoment2d(5,5)
```
"""
function simplemoment2d(n::Int, m::Int; T=Float64)
  mom = JuSwalbe.Macroquant(height=ones(T, (n,m)), velocity=JuSwalbe.Twovector(x=fill(T(0.1),(n,m)), y=fill(T(0.1),(n,m))), pressure=zeros(T, (n,m)), energy=zeros(T,(n,m)))
  return mom
end