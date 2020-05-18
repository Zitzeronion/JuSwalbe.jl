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
Such the collision process in the lattice Boltzmann method tries to equilibrate the system.

If there is an additional force acting on the fluid, i.e. `gravity` than this part is usually included in the numerical implementation of the collision operator.
The forcing vector of source term is often called `` \\mathbf{S}_{\\alpha} `` and given by

`` \\mathbf{S}_{\\alpha} = \\Delta t \\frac{3 w_{\\alpha}}{\\mathvbf{c}_{\\alpha}^2}\\mathbf{c}_{\\alpha,i} \\mathbf{F}_i. ``


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

## Two spatial dimensions
```jldoctest
julia> using JuSwalbe, Test

julia> input, mom, force, dist = minimalsetup2d(10,5)
(JuSwalbe.Inputconstants
  lx: Int64 10
  ly: Int64 5
  maxruntime: Int64 100000
  dumping: Int64 1000
  τ: Float64 1.0
  gravity: Float64 0.0
  γ: Float64 0.01
  δ: Float64 1.0
  μ: Float64 0.16666666666666666
, JuSwalbe.Macroquant{Array{Float64,2},JuSwalbe.Twovector{Array{Float64,2}}}
  height: Array{Float64}((10, 5)) [1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0; … ; 1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0]
  velocity: JuSwalbe.Twovector{Array{Float64,2}}
  pressure: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  energy: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
, JuSwalbe.Forces{JuSwalbe.Twovector{Array{Float64,2}}}
  slip: JuSwalbe.Twovector{Array{Float64,2}}
  h∇p: JuSwalbe.Twovector{Array{Float64,2}}
  bathymetry: JuSwalbe.Twovector{Array{Float64,2}}
  thermal: JuSwalbe.Twovector{Array{Float64,2}}
, JuSwalbe.DistributionD2Q9{Array{Float64,2}}
  f0: Array{Float64}((10, 5)) [1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0; … ; 1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0]
  f1: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f2: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f3: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f4: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f5: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f6: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f7: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f8: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
)

julia> dist
JuSwalbe.DistributionD2Q9{Array{Float64,2}}
  f0: Array{Float64}((10, 5)) [1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0; … ; 1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0]
  f1: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f2: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f3: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f4: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f5: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f6: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f7: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f8: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]

julia> equi = calc_equilibrium_distribution(mom, gravity=input.gravity)
JuSwalbe.DistributionD2Q9{Array{Float64,2}}
  f0: Array{Float64}((10, 5)) [0.9994666666666666 0.9994666666666666 … 0.9994666666666666 0.9994666666666666; 0.9994666666666666 0.9994666666666666 … 0.9994666666666666 0.9994666666666666; … ; 0.9994666666666666 0.9994666666666666 … 0.9994666666666666 0.9994666666666666; 0.9994666666666666 0.9994666666666666 … 0.9994666666666666 0.9994666666666666]
  f1: Array{Float64}((10, 5)) [0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333; 0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333; … ; 0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333; 0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333]
  f2: Array{Float64}((10, 5)) [-0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999; -0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999; … ; -0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999; -0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999]
  f3: Array{Float64}((10, 5)) [-0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999; -0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999; … ; -0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999; -0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999]
  f4: Array{Float64}((10, 5)) [0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333; 0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333; … ; 0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333; 0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333]
  f5: Array{Float64}((10, 5)) [-3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5; -3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5; … ; -3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5; -3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5]
  f6: Array{Float64}((10, 5)) [-0.0031666666666666666 -0.0031666666666666666 … -0.0031666666666666666 -0.0031666666666666666; -0.0031666666666666666 -0.0031666666666666666 … -0.0031666666666666666 -0.0031666666666666666; … ; -0.0031666666666666666 -0.0031666666666666666 … -0.0031666666666666666 -0.0031666666666666666; -0.0031666666666666666 -0.0031666666666666666 … -0.0031666666666666666 -0.0031666666666666666]
  f7: Array{Float64}((10, 5)) [-3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5; -3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5; … ; -3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5; -3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5]
  f8: Array{Float64}((10, 5)) [0.0034999999999999996 0.0034999999999999996 … 0.0034999999999999996 0.0034999999999999996; 0.0034999999999999996 0.0034999999999999996 … 0.0034999999999999996 0.0034999999999999996; … ; 0.0034999999999999996 0.0034999999999999996 … 0.0034999999999999996 0.0034999999999999996; 0.0034999999999999996 0.0034999999999999996 … 0.0034999999999999996 0.0034999999999999996]

julia> force.slip = JuSwalbe.Twovector(x=zeros(10,5), y=zeros(10,5))
JuSwalbe.Twovector{Array{Float64,2}}
  x: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  y: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]

julia> force.h∇p = JuSwalbe.Twovector(x=zeros(10,5), y=zeros(10,5))
JuSwalbe.Twovector{Array{Float64,2}}
  x: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  y: Array{Float64}((10, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]

julia> newdist = collisionBGK(mom, force, dist, input)
JuSwalbe.DistributionD2Q9{Array{Float64,2}}
  f0: Array{Float64}((10, 5)) [0.9994666666666666 0.9994666666666666 … 0.9994666666666666 0.9994666666666666; 0.9994666666666666 0.9994666666666666 … 0.9994666666666666 0.9994666666666666; … ; 0.9994666666666666 0.9994666666666666 … 0.9994666666666666 0.9994666666666666; 0.9994666666666666 0.9994666666666666 … 0.9994666666666666 0.9994666666666666]
  f1: Array{Float64}((10, 5)) [0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333; 0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333; … ; 0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333; 0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333]
  f2: Array{Float64}((10, 5)) [-0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999; -0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999; … ; -0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999; -0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999]
  f3: Array{Float64}((10, 5)) [-0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999; -0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999; … ; -0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999; -0.006599999999999999 -0.006599999999999999 … -0.006599999999999999 -0.006599999999999999]
  f4: Array{Float64}((10, 5)) [0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333; 0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333; … ; 0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333; 0.006733333333333333 0.006733333333333333 … 0.006733333333333333 0.006733333333333333]
  f5: Array{Float64}((10, 5)) [-3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5; -3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5; … ; -3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5; -3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5]
  f6: Array{Float64}((10, 5)) [-0.0031666666666666666 -0.0031666666666666666 … -0.0031666666666666666 -0.0031666666666666666; -0.0031666666666666666 -0.0031666666666666666 … -0.0031666666666666666 -0.0031666666666666666; … ; -0.0031666666666666666 -0.0031666666666666666 … -0.0031666666666666666 -0.0031666666666666666; -0.0031666666666666666 -0.0031666666666666666 … -0.0031666666666666666 -0.0031666666666666666]
  f7: Array{Float64}((10, 5)) [-3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5; -3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5; … ; -3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5; -3.3333333333333335e-5 -3.3333333333333335e-5 … -3.3333333333333335e-5 -3.3333333333333335e-5]
  f8: Array{Float64}((10, 5)) [0.0034999999999999996 0.0034999999999999996 … 0.0034999999999999996 0.0034999999999999996; 0.0034999999999999996 0.0034999999999999996 … 0.0034999999999999996 0.0034999999999999996; … ; 0.0034999999999999996 0.0034999999999999996 … 0.0034999999999999996 0.0034999999999999996; 0.0034999999999999996 0.0034999999999999996 … 0.0034999999999999996 0.0034999999999999996]

julia> @test equi.f0 == newdist.f0 # No force and τ=1 means fnew has to be feq!
Test Passed

julia> @test equi.f2 == newdist.f2
Test Passed
```


# References
In general a very good reference concerning the method at all is the book written by Timm Krüger
-[The Lattice Boltzmann Method: Principles and Practice](https://www.springer.com/gp/book/9783319446479)

## One spatial dimension
Here is the work by Chopard et al. really worth reading
- [Asymmetric lattice Boltzmann model for shallow water flows](https://www.sciencedirect.com/science/article/abs/pii/S0045793013003599)

## Two spatial dimensions
In two dimensions there is the problematic of how to add forces.
One way to justify is to do a Chapman Enskog expansion and match the shallow water equations with a source term.
Apparently there is not a unique solution but a class of possible forcings.
The first approach by Salmon does not include the lattice weights, which results in isotropy problems.
On the other hand Zhous approach is not as acqurate as he says, which was pointed out by Chopard.
Nevertheless here are some references educate yourself which is the fitting one
-[The lattice Boltzmann method as a basis for ocean circulation modeling, by Salmon](https://www.ingentaconnect.com/content/jmr/jmr/1999/00000057/00000003/art00005)
-[Lattice Boltzmann Methods for Shallow Water Flows, by Zhou](https://www.springer.com/gp/book/9783540407461)
-[An evaluation of force terms in the lattice Boltzmann models in simulating shallow water flows over complex topography, by Peng](https://onlinelibrary.wiley.com/doi/pdf/10.1002/fld.4726)

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
  # println(allforces)

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
  # Collision operation for a D2Q9 shallow water lattice Boltzmann.
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

The so called streaming process transports the distribution functions according to thier lattice velocity.
In one dimension this means that the zero velocity distribution stays at its lattice point while the positive velocity advances by one lattice node.
Of course the negative velocity component `.f2` decreases its lattice index by one.

In two dimensions there are nine of these operations according to the `D2Q9` lattice vectors.

# Math
After the collision process the newly calculated distribution functions need to be streamed.
They performe one interation or streaming per time update.
Such the distributions after collision `` f^{\\ast}_i `` are transported to their respective new position

`` f^{\\ast}_i(\\mathbf{x} + \\mathbf{c}_i \\Delta t)``

One has to be careful in the proximity of boundaries, because they could destroy mass which is unphysical.

# Example
```jldoctest
julia> using JuSwalbe

julia> dist = JuSwalbe.DistributionD2Q9(f0=ones(5,5), f1=fill(0.1, (5,5)), f2=fill(0.1, (5,5)),
                                        f3=fill(0.1, (5,5)), f4=fill(0.1, (5,5)), f5=zeros(5,5),
                                        f6=zeros(5,5), f7=zeros(5,5), f8=zeros(5,5))
JuSwalbe.DistributionD2Q9{Array{Float64,2}}
  f0: Array{Float64}((5, 5)) [1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0; … ; 1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0]
  f1: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f2: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f3: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f4: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f5: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f6: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f7: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f8: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
   
julia> dist.f1[2,2] = 1; dist.f2[2,2] = 1; dist.f5[2,2] = 1; dist.f3[2,2] = 1; dist.f8[2,2] = 1 
1

julia> newdist = streamdistperiodic!(dist)
JuSwalbe.DistributionD2Q9{Array{Float64,2}}
  f0: Array{Float64}((5, 5)) [1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0; … ; 1.0 1.0 … 1.0 1.0; 1.0 1.0 … 1.0 1.0]
  f1: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f2: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f3: Array{Float64}((5, 5)) [0.1 1.0 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f4: Array{Float64}((5, 5)) [0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1; … ; 0.1 0.1 … 0.1 0.1; 0.1 0.1 … 0.1 0.1]
  f5: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f6: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f7: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]
  f8: Array{Float64}((5, 5)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0]

julia> dist.f1
5×5 Array{Float64,2}:
 0.1  0.1  0.1  0.1  0.1
 0.1  1.0  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1

julia> newdist.f1
5×5 Array{Float64,2}:
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  1.0  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1
 0.1  0.1  0.1  0.1  0.1

julia> dist.f8
5×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
 0.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> newdist.f8
5×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 
```

# References
Any book on the lattice Boltzmann method should do.
- [The Lattice Boltzmann Method: Principles and Practice](https://www.springer.com/gp/book/9783319446479)

See also [`collisionBGK`](@ref)
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

function streamdistperiodic!(dist::JuSwalbe.DistributionD2Q9{Matrix{T}}) where {T<:Number}
  width, thick = size(dist.f0)
  newdist = zeros(T, (width, thick, 9))
  # Shifting the distributions by the nine lattice velocities
  distarray = dist2array(dist)
  latvel = [(0, 0); (1, 0); (0, 1); (-1, 0); (0, -1); (1, 1); (-1, 1); (-1, -1); (1, -1)]
  for i in 2:9
    newdist[:, :, i] .= circshift(distarray[i, :, :], latvel[i])
  end

  dist = JuSwalbe.DistributionD2Q9(f0=dist.f0, f1=newdist[:, :, 2], f2=newdist[:, :, 3],
                                   f3=newdist[:, :, 4], f4=newdist[:, :, 5], f5=newdist[:, :, 6], 
                                   f6=newdist[:, :, 7], f7=newdist[:, :, 8], f8=newdist[:, :, 9])

  return dist
end