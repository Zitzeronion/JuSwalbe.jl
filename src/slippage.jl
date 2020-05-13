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

"""
function computeslip(mom::JuSwalbe.Macroquant{Vector{T}, Vector{T}}, forces::JuSwalbe.Forces{Vector{T}}, input::JuSwalbe.Inputconstants) where {T<:Number}
    # Measure the length and allocate a dummy array
    len = length(mom.height)
    slippage = zeros(T, len)
    # Slippage can be calculated assuming a parabolic velocity profile in z
    slippage = input.μ * (6 * mom.height .* mom.velocity ./ (2 * mom.height.^2 .+ 6 * input.δ * mom.height .+ 6 * input.δ^2))
    # Write the result into the forcing array
    forces.slip = slippage
    # Write it out as well
    return slippage
end

function computeslip(mom::JuSwalbe.Macroquant{Matrix{T}, JuSwalbe.Twovector{Matrix{T}}}, forces::JuSwalbe.Forces{JuSwalbe.Twovector{Matrix{T}}}, δ::T=T(1.0), μ::T=T(1/6)) where {T<:Number}
    # Measure the length and allocate a dummy array
    width, thick = size(mom.height)
    slippagex = zeros(T, (width, thick))
    slippagey = zeros(T, (width, thick))
    # Slippage can be calculated assuming a parabolic velocity profile in z
    slippagex = μ * (6 * mom.height .* mom.velocity.x ./ (2 * mom.height.^2 .+ 6 * δ * mom.height .+ 6 * δ^2))
    slippagey = μ * (6 * mom.height .* mom.velocity.x ./ (2 * mom.height.^2 .+ 6 * δ * mom.height .+ 6 * δ^2))
    # Write the result into the forcing array
    slippage = JuSwalbe.Twovector(x=slippagex, y=slippagey)
    forces.slip = slippage
    # Write it out as well
    return slippage
end