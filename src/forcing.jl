"""
    computeslip(mom, δ, force)

Calculates the force due to slippage on the substrate.

Since the model only alows for forces in `x` and `y` direction but not in `z` it is necessary to regularize the velocity at `` h \\rightarrow 0 ``.
The so called no-slip boundary condition can be in theory implemented using `δ` = 0. 
However as soon as the film dewets and `height` becomes small, the solution diverges.
Choosing `δ` = 1 satisfies the weak slip condition, while choosing it much larger than 1 yields the strong slip solution.

# Math
One and two dimensional approaches are essentially equal.
Behind the force is the idea of a parabolic shape of the velocity along `z`.
With the boundary conditions

`` \\mathbf{u}|_{h=H} = \\mathbf{u}  ``

`` \\mathbf{u}|_{h=0} = 0 ``

The parameter `δ` allows to push the zero velocity into the substrate, such below te computable region.
Solving this equation with some more assumptions yields

`` \\mathbf{F}_{slip} = \\mu\\alpha_{\\delta}(h) \\mathbf{u} = \\frac{6\\mu h \\mathbf{u}}{2 h^2 + 6\\delta h + 3\\delta^2}. ``

# Example 

# References
Here I shamelessly cite the paper from within our group 
- [Lattice Boltzmann method for thin-liquid-film hydrodynamics](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.033313)

See also: [`Inputconstants`](@ref), [`calculatemoments`](@ref)

"""
function computeslip(mom::JuSwalbe.Macroquant{Vector{T}, Vector{T}}, forces::JuSwalbe.Forces{Vector{T}}, input::JuSwalbe.Inputconstants) where {T<:Number}
    # Measure the length and allocate a dummy array
    len = length(mom.height)
    slippage = zeros(T, len)
    num = zeros(T, len)
    denom = zeros(T, len)
    # Slippage can be calculated assuming a parabolic velocity profile in z
    num .= 6 .* mom.height .* mom.velocity
    denom .= 2 .* mom.height .+ 6 * mom.height * input.δ .+ 3 * input.δ^2
    slippage .= input.μ * (num ./ denom)
    # Write the result into the forcing array
    forces.slip = slippage
    # Write it out as well
    return slippage
end

function computeslip(mom::JuSwalbe.Macroquant{Matrix{T}, JuSwalbe.Twovector{Matrix{T}}}, forces::JuSwalbe.Forces{JuSwalbe.Twovector{Matrix{T}}}, input::JuSwalbe.Inputconstants) where {T<:Number}
    # Measure the length and allocate a dummy array
    width, thick = size(mom.height)
    δ, μ = T(input.δ), T(input.μ)
    slippagex = zeros(T, (width, thick))
    slippagey = zeros(T, (width, thick))
    numx = zeros(T, (width, thick))
    numy = zeros(T, (width, thick))
    denom = zeros(T, (width, thick))
    # Slippage can be calculated assuming a parabolic velocity profile in z
    numx .= 6 .* mom.height .* mom.velocity.x
    numy .= 6 .* mom.height .* mom.velocity.y
    denom .= 2 .* mom.height .+ 6 * mom.height * δ .+ 3 * δ^2
    slippagex .= μ * ( numx ./ denom )
    slippagey .= μ * ( numy ./ denom )
    # Write the result into the forcing array
    slippage = JuSwalbe.Twovector(x=slippagex, y=slippagey)
    forces.slip = slippage
    # Write it out as well
    return slippage
end

"""
    computethermalcapillarywaves(mom, force, input)

Computes the forcing arising due to the inclusion of thermal capillary waves that undulate the free surface.

Adds a normal distributed term weighted with the thermal energy kbt to the forcing struct.
This addition makes it possible to simulate not only the thin film equation but also the stochastic thin film equation.

# Math
The stochastic thin film equation (STF) much like the thin film equation is derived from the Landau Lifshitz Navier Stokes equation (LLNS).
Assuming the fluctuations are small and not leading order one can perform an integration of the stochastic stresses, which was done by Grün et al.,

`` \\frac{\\partial h}{\\partial t} = \\nabla\\cdot\\Big[\\frac{h^3}{3\\mu}\\nabla(\\Pi(h) - \\gamma \\Delta h) + \\int_0^h (h-y)\\mathcal{S}_{z||}(y)dy]. ``

Taking into account several assumption it can be shown that the above equation is equal to 

`` \\frac{\\partial h}{\\partial t} = \\nabla\\cdot\\Big[\\frac{h^3}{3\\mu}\\nabla(\\Pi(h) - \\gamma \\Delta h) + \\sqrt{\\frac{2k_BT}{3\\mu}h^3}\\mathcal{N}(t)]. ``

Such the integral of `` \\mathcal{S} `` can be expressed with a single multipicative noise term `` \\mathcal{N} `` which is gaussian distributed with zero mean and a variance of one.

This force is adding exactly the term with `` \\mathcal{N} `` to the forcing struct, however it is a little modified due to slippage.

# Example
```jldoctest
julia> using JuSwalbe

julia> mom = simplemoment1d(1000)
JuSwalbe.Macroquant{Array{Float64,1},Array{Float64,1}}
  height: Array{Float64}((1000,)) [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0  …  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  velocity: Array{Float64}((1000,)) [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1  …  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
  pressure: Array{Float64}((1000,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  energy: Array{Float64}((1000,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

julia> force = JuSwalbe.Forces(slip=fill(0.1. 1000), thermal=zeros(1000), h∇p=zeros(1000), bathymetry=zeros(1000))
JuSwalbe.Forces{Array{Float64,1}}
  slip: Array{Float64}((1000,)) [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1  …  0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
  h∇p: Array{Float64}((1000,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  bathymetry: Array{Float64}((1000,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
  thermal: Array{Float64}((1000,)) [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

julia> input = Inputconstants()
JuSwalbe.Inputconstants
  lx: Int64 512
  ly: Int64 512
  maxruntime: Int64 100000
  dumping: Int64 1000
  τ: Float64 1.0
  gravity: Float64 0.0
  γ: Float64 0.01
  δ: Float64 1.0
  μ: Float64 0.16666666666666666
  kbt: Float64 0.0

julia> thermal = computethermalcapillarywaves(mom, force, input); # No k_BT -> no force!

julia> input2 = JuSwalbe.Inputconstants(kbt = 0.001)
JuSwalbe.Inputconstants
  lx: Int64 512
  ly: Int64 512
  maxruntime: Int64 100000
  dumping: Int64 1000
  τ: Float64 1.0
  gravity: Float64 0.0
  γ: Float64 0.01
  δ: Float64 1.0
  μ: Float64 0.16666666666666666
  kbt: Float64 0.001

julia> thermal = computethermalcapillarywaves(mom, force, input2);
```
"""
function computethermalcapillarywaves(mom::JuSwalbe.Macroquant{Vector{T},Vector{T}}, forces::JuSwalbe.Forces{Vector{T}}, input::Inputconstants) where {T<:Number}
    len = length(mom.height)
    thermocap = zeros(T, len)
    slip = deepcopy(forces.slip)
    μ = T(input.μ)
    kbt = T(input.kbt)
    # Generate a Gaussian distribution with zero mean and variance of one 
    gaussian = Normal()
    # Fill an array with random numbers distributed according to a Normal
    gaussianvec = rand(gaussian, len)
    # Compute the forces due to thermal capillary waves
    thermocap .= sqrt.(2 * kbt * μ * slip) .* gaussianvec

    forces.thermal = thermocap
    return thermocap
end

function computethermalcapillarywaves(mom::JuSwalbe.Macroquant{Matrix{T},JuSwalbe.Twovector{Matrix{T}}}, forces::JuSwalbe.Forces{JuSwalbe.Twovector{Matrix{T}}}, input::Inputconstants) where {T<:Number}
    width, thick = size(mom.height)
    thermocap = zeros(T, (width, thick))
    slip = deepcopy(forces.slip)
    μ = T(input.μ)
    kbt = T(input.kbt)
    # Generate a Gaussian distribution with zero mean and a variance of one 
    gaussian = Normal()
    # Fill an array with random numbers distributed according to a gaussian
    gaussianvec = rand(gaussian, (width, thick))
    # Compute the forces due to thermal capillary waves
    thermocap .= sqrt.(2 * kbt * μ * slip) .* gaussianvec

    forces.thermal = thermocap
    return thermocap
end