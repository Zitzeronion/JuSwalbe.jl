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