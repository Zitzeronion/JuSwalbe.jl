"""
    pressure(mom::JuSwalbe.macroquant64_2d; γ = 0.01, θ = 1/9)

Film pressure of the thin film equation.

This is just the summation of the laplacian of the height `Δh` and the disjoining potential `Π`.
See also: [`Δh`](@ref), [`Π`](@ref)

# Math
The film pressure is simply

`` p_{film} = -\\gamma\\Delta h + \\Pi(h) ``

The sign changes from paper to paper.
Here I choose to be in agreement with Thiele et al., but the sign here is not straightforward. 
I usually choose the sign based on the solution (yes try and error).

# Example
```jldoctest
julia> using JuSwalbe

julia> height = reshape([i for i in 1.0:1.0:16.0],4,4)
4×4 Array{Float64,2}:
 1.0  5.0   9.0  13.0
 2.0  6.0  10.0  14.0
 3.0  7.0  11.0  15.0
 4.0  8.0  12.0  16.0

julia> p = pressure(height)
4×4 Array{Float64,2}:
 -0.200016  -0.0400001   -0.04        0.12
 -0.160002  -7.44536e-8  -1.6082e-8   0.16
 -0.160001  -4.68862e-8  -1.20826e-8  0.16
 -0.12       0.04         0.04        0.2

```
So what does this mean?
A positive pressure is force that drives the film down in height.
While a negative pressure generates a flux towards that location.
In the example a linear increasing height field was used, an equilibrium though is reached when the film is flat.

# References
TBD

"""
function pressure(mom::JuSwalbe.macroquant64_2d; γ::Float64=0.01, θ::Float64=1.0/9.0)
    # Dummy array to store the result
    p = zeros(size(mom.height))
    # All calculation needed here
    p = -γ * Δh(mom) .+ Π(mom, γ=γ, θ=θ)

    return p
end
# For multiple dispatch
function pressure(height::Array{Float64,2}; γ::Float64=0.01, θ::Float64=1.0/9.0)
    # Dummy array to store the result
    p = zeros(size(height))
    # All calculation needed here
    p = -γ * Δh(height) .+ Π(height, γ=γ, θ=θ)

    return p
end

"""
    Π(mom::JuSwalbe.macroquant64_2d, exponents = [9,3], γ = 0.01, θ = 1.0/9.0)

Calculates the disjoining pressure potential for a given surface tension and contact angle at every lattice point.

The disjoing pressure potential enables the simulation to study phenomena like dewetting of thin liquid films.
Only with this potential it is possible to "dewett" without using highly sophisticated boundary conditions.
Although the functional shape of this potential can vary the exponents (9,3) mimic the Lenard-Jones potential for liquid substrate interaction.

# Arguments
- `mom::JuSwalbe.macroquant64_2d` : Macroscopic moments of the simulation, important here `.height` which contains the height field
- `h_star::Float64` : Minimum of the pontential.
- `exponents::Array{Int64,1}` : Stores the exponent for the powerlaw potential
- `γ::Float64` : Surface tension
- `θ::Float64` : Equilibrium contact angle in multiple of π (for accuarcy reasons)

# Math
The potential can be derived from the assumption that a given surface energy demands an equilibrium contact angle for the fluid.
For the exact derivation take a look at the references.

`` \\Pi(h) = \\kappa f(h) = (1-\\cos(\\theta))\\frac{(n-1)(m-1)}{(n-m)h_{\\ast}}\\left[\\left(\frac{h_{\\ast}}{h}\\right)^n - \\left(\frac{h_{\\ast}}{h}\\right)^m\\right] ``

`` \\Pi(h_{\\ast}) = 0 `` 

# Example
```jldoctest
julia> using JuSwalbe

julia> n,m = (4,4)
(4, 4)

julia> moment = JuSwalbe.macroquant64_2d(ones(n,m), JuSwalbe.velocity64_2d(zeros(n,m),zeros(n,m)), zeros(n,m))
JuSwalbe.macroquant64_2d([1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0], JuSwalbe.velocity64_2d([0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]), [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0])

julia> p = Π(moment.height)
4×4 Array{Float64,2}:
 -1.6082e-5  -1.6082e-5  -1.6082e-5  -1.6082e-5
 -1.6082e-5  -1.6082e-5  -1.6082e-5  -1.6082e-5
 -1.6082e-5  -1.6082e-5  -1.6082e-5  -1.6082e-5
 -1.6082e-5  -1.6082e-5  -1.6082e-5  -1.6082e-5

julia> h = reshape([i for i in 1.0:1.0:(n*m)],n,m)
4×4 Array{Float64,2}:
 1.0  5.0   9.0  13.0
 2.0  6.0  10.0  14.0
 3.0  7.0  11.0  15.0
 4.0  8.0  12.0  16.0

julia> p = Π(h)
4×4 Array{Float64,2}:
 -1.6082e-5   -1.28656e-7  -2.20603e-8  -7.31997e-9
 -2.01025e-6  -7.44536e-8  -1.6082e-8   -5.86078e-9
 -5.95628e-7  -4.68862e-8  -1.20826e-8  -4.76503e-9
 -2.51281e-7  -3.14101e-8  -9.30669e-9  -3.92626e-9

```

# References
To get a good understanding:
- [Long-scale evolution of thin liquid films](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.69.931)
- [Wetting and spreading](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.81.739)
- [Dynamics and stability of thin liquid films](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.81.1131)

A rather recent new setup for the shape of `Π` can be found in 
- [Signatures of slip in dewetting polymer films](https://www.pnas.org/content/116/19/9275) 
""" 
function Π(mom::JuSwalbe.macroquant64_2d; h_star=0.1, exponents=[9,3], γ=0.01, θ=1.0/9.0)
    # Theoretical minimum of the wetting potential
    
    # Actual formular of the disjoining potential, long range attracion short range repulsion.
    Π_h = (γ*(1-cospi(θ)) 
          * (exponents[1] - 1)*(exponents[2] - 1)/((exponents[1] - exponents[2])*h_star) 
          * ((h_star ./ mom.height).^exponents[1] .- (h_star ./ mom.height).^exponents[2]))
    
    return Π_h 
end

function Π(height::Array{Float64,2}; h_star=0.1, exponents=[9,3], γ=0.01, θ=1.0/9.0)
    # Theoretical minimum of the wetting potential
    
    # Actual formular of the disjoining potential, long range attracion short range repulsion.
    Π_h = (γ*(1-cospi(θ)) 
          * (exponents[1] - 1)*(exponents[2] - 1) / ((exponents[1] - exponents[2])*h_star) 
          * ((h_star ./ height).^exponents[1] .- (h_star ./ height).^exponents[2]))
    
    return Π_h 
end

"""
    Δh(mom::JuSwalbe.macroquant64_2d)

Calculates the laplacian of the height field with periodic boundaries.

The calculation of the laplacian is central for the thin film evolution.
Only the pressure gradient will induce flow, at least for the non-fluctuating version.
Therefore it is fairly important to have an accurate computation of the laplacian.

# Math

## Two spatial dimensions

In two dimensions I follow the paper from Santosh and Succi.
A nine-point stencil is used and the equation goes as follow

`` \\Delta h_{i,j} = \\frac{1}{6}[4(h_{i+1,j}+h_{i,j+1}+h_{i-1,j}+h_{i,j-1}) + (h_{i+1,j+1}+h_{i+1,j-1}+h_{i-1,j+1}+h_{i+1,j-1}) - 20 h_{i,j}] ``

Where `i` and `j` are the x and y coordinates.
Periodicity is taken care of by a circular padding of the height array.

# Example
```jldoctest
julia> using JuSwalbe

julia> moment = JuSwalbe.macroquant64_2d(reshape(1:16, 4,4), JuSwalbe.velocity64_2d(zeros(4,4), zeros(4,4)), zeros(4,4))
JuSwalbe.macroquant64_2d([1.0 5.0 9.0 13.0; 2.0 6.0 10.0 14.0; 3.0 7.0 11.0 15.0; 4.0 8.0 12.0 16.0], JuSwalbe.velocity64_2d([0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]), [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0])

julia> moment.height
4×4 Array{Float64,2}:
 1.0  5.0   9.0  13.0
 2.0  6.0  10.0  14.0
 3.0  7.0  11.0  15.0
 4.0  8.0  12.0  16.0

julia> laplace = Δh(moment)
4×4 Array{Float64,2}:
 20.0   4.0   4.0  -12.0
 16.0   0.0   0.0  -16.0
 16.0   0.0   0.0  -16.0
 12.0  -4.0  -4.0  -20.0
```

# References
## Two spatial dimensions
There are plenty of papers concerning discret differentail operators.
Many of them are good, although definitly not the best I go here with:
- [Isotropic discrete Laplacian operators from lattice hydrodynamics](https://www.sciencedirect.com/science/article/pii/S0021999112004226)

"""
function Δh(mom::JuSwalbe.macroquant64_2d)
    # Get the size of the problem
    width, thick = size(mom.height)
    # Build a wrap around to mimic periodic boundaries
    wrapper = 1
    toruswrapped = padarray(mom.height, Pad(:circular, wrapper, wrapper))
    # Generate a dummy array for the solution
    Δ = zeros(Float64, size(toruswrapped))
    # Exclude the wrapped boundaries
    for j in 1:width, i in 1:thick
        @inbounds Δ[i,j] = 1.0/6.0 * (4.0 * (toruswrapped[i+1, j] + toruswrapped[i-1, j] + toruswrapped[i, j+1] + toruswrapped[i, j-1])
                                     + 1.0 * (toruswrapped[i+1, j+1] + toruswrapped[i+1, j-1] + toruswrapped[i-1, j+1] + toruswrapped[i-1, j-1])
                                     - 20.0 * toruswrapped[i,j])
    end
    # Return the result without the wrapped boundaries
    return Δ[1:(end-2), 1:(end-2)]
end

function Δh(height::Array{Float64,2})
    # Get the size of the problem
    width, thick = size(height)
    # Build a wrap around to mimic periodic boundaries
    wrapper = 1
    toruswrapped = padarray(height, Pad(:circular, wrapper, wrapper))
    # Generate a dummy array for the solution
    Δ = zeros(Float64, size(toruswrapped))
    # Exclude the wrapped boundaries
    for j in 1:width, i in 1:thick
        @inbounds Δ[i,j] = 1.0/6.0 * (4.0 * (toruswrapped[i+1, j] + toruswrapped[i-1, j] + toruswrapped[i, j+1] + toruswrapped[i, j-1])
                                     + 1.0 * (toruswrapped[i+1, j+1] + toruswrapped[i+1, j-1] + toruswrapped[i-1, j+1] + toruswrapped[i-1, j-1])
                                     - 20.0 * toruswrapped[i,j])
    end
    # Return the result without the wrapped boundaries
    return Δ[1:(end-2), 1:(end-2)]
end