"""
    pressure(mom::JuSwalbe.Macroquant; γ = 0.01, θ = 1/9)

Film pressure of the thin film equation.

This is just the summation of the laplacian of the height `Δh` and the disjoining potential `Π`.
Since this function uses Π one has to be careful of θ as it should be an array.
However it can be either a single element array or as large as the whole domain.
See also: [`Δh`](@ref), [`Π`](@ref)

# Math
The film pressure is simply

`` p_{film} = -\\gamma\\Delta h + \\Pi(h) ``

The sign is not 100% fixed and can change from paper to paper.
Here I choose to be in agreement with Thiele et al., which is known to be a good theoretician. 


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
## Recent (short)
- [Signatures of slip in dewetting polymer films](https://www.pnas.org/content/116/19/9275)

## Review
- [Dynamics and stability of thin liquid films](https://journals.aps.org/rmp/abstract/10.1103/RevModPhys.81.1131)
"""
function pressure(mom::JuSwalbe.Macroquant{Matrix{T}, JuSwalbe.Twovector{Matrix{T}}}; γ::T=T(0.01), θ::Matrix{T}=ones(T,(1,1))*T(1/9)) where {T<:Number}
    # Dummy array to store the result
    p = zeros(T, size(mom.height))
    # All calculation needed here
    p = -γ * Δh(mom) .+ Π(mom, γ=γ, θ=θ)

    return p
end

function pressure(height::Array{T,2}; γ::T=0.01, θ::Matrix{T}=ones(T,(1,1))*T(1/9)) where {T<:Number}
    # Dummy array to store the result
    p = zeros(size(height))
    # All calculation needed here
    p = -γ * Δh(height) .+ Π(height, γ=γ, θ=θ)

    return p
end
# One dimensional implementation
function pressure(mom::JuSwalbe.Macroquant{Vector{T}, Vector{T}}; γ::T=T(0.01), θ::Vector{T}=ones(T,1)*T(1/9)) where {T<:Number}
    # Dummy array to store the result
    p = zeros(T, length(mom.height))
    # All calculation needed here
    p = -γ * Δh(mom) .+ Π(mom, γ=γ, θ=θ)

    return p
end

function pressure(height::Array{T,1}; γ::T=0.01, θ::Vector{T}=ones(T,1)*T(1/9)) where {T<:Number}
    # Dummy array to store the result
    p = zeros(size(height))
    # All calculation needed here
    p = -γ * Δh(height) .+ Π(height, γ=γ, θ=θ)

    return p
end

"""
    Π(mom::JuSwalbe.Macroquant; h_star = 0.1, exponents = [9,3], γ = 0.01, θ = 1.0/9.0)

Calculates the disjoining pressure for a given surface tension and contact angle at every lattice point.

The disjoing pressure potential enables the simulation to study phenomena like dewetting of thin liquid films.
Only with this potential it is possible to "dewett" without using highly sophisticated boundary conditions at `` h\\rightarrow 0 ``.
The functional shape of this potential can vary, ranging from the powerlaw used here to exponentials to a mixture of both.
Using the powerlaw shape the exponents (9,3) mimic the Lenard-Jones potential for the liquid substrate interaction.

# Arguments
- `mom::JuSwalbe.Macroquant` : Macroscopic moments of the simulation, important here `.height` which contains the height field. Allows for array input as well.
- `h_star::T` : Minimum of the pontential. T needs to be a subtype of `Number`! 
- `exponents::Array{Int64,1}` : Exponents for the powerlaw potential
- `γ::T` : Surface tension
- `θ::T` : Equilibrium contact angle in multiple of π (for accuarcy reasons). Needs to be an array when doing patterning. For simple substrate use a one entry array with correct dimension.

# Math
The potential can be derived from the assumption that a given surface energy demands an equilibrium contact angle for the fluid.
For the exact derivation take a look at the references, in principle it is the derivative of the interfacial potential with respect to `h`.

`` \\Phi'(h) = \\Pi(h) `` 

`` \\Pi(h) = \\kappa f(h) = (1-\\cos(\\theta))\\frac{(n-1)(m-1)}{(n-m)h_{\\ast}}\\Bigg[\\Bigg(\frac{h_{\\ast}}{h}\\Bigg)^n - \\Bigg(\frac{h_{\\ast}}{h}\\Bigg)^m\\Bigg] ``

`` \\Pi(h_{\\ast}) = 0 `` 

# Example
```jldoctest
julia> using JuSwalbe

julia> n,m = (4,4)
(4, 4)

julia> moment = JuSwalbe.Macroquant{Matrix{Float64}, JuSwalbe.Twovector{Matrix{Float64}}}(ones(n,m), JuSwalbe.Twovector(zeros(n,m),zeros(n,m)), zeros(n,m),zeros(n,m))
JuSwalbe.Macroquant{Array{Float64,2},JuSwalbe.Twovector{Array{Float64,2}}}([1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0; 1.0 1.0 1.0 1.0], JuSwalbe.Twovector{Array{Float64,2}}([0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]), [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0])

julia> p = Π(moment)
4×4 Array{Float64,2}:
 -1.6082e-5  -1.6082e-5  -1.6082e-5  -1.6082e-5
 -1.6082e-5  -1.6082e-5  -1.6082e-5  -1.6082e-5
 -1.6082e-5  -1.6082e-5  -1.6082e-5  -1.6082e-5
 -1.6082e-5  -1.6082e-5  -1.6082e-5  -1.6082e-5

julia> p = Π(moment, θ=zeros(1,1)) # Fully wetting substrate
4×4 Array{Float64,2}:
 -0.0  -0.0  -0.0  -0.0
 -0.0  -0.0  -0.0  -0.0
 -0.0  -0.0  -0.0  -0.0
 -0.0  -0.0  -0.0  -0.0

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
function Π(mom::JuSwalbe.Macroquant{Matrix{T}, JuSwalbe.Twovector{Matrix{T}}}; h_star::T=T(0.1), exponents=[9,3], γ::T=T(0.01), θ::Matrix{T}=ones(T,(1,1))*T(1/9)) where {T<:Number}   
    # Theoretical minimum of the wetting pontential
    
    # Actual formular of the disjoining potential, long range attracion short range repulsion.
    Π_h = (γ .* (1 .- cospi.(θ)) 
          .* (exponents[1] - 1)*(exponents[2] - 1)/((exponents[1] - exponents[2])*h_star) 
          .* ((h_star ./ mom.height).^exponents[1] .- (h_star ./ mom.height).^exponents[2]))
    
    return Π_h 
end

function Π(height::Array{T,2}; h_star::T=T(0.1), exponents=[9,3], γ::T=T(0.01), θ::Matrix{T}=ones(T,(1,1))*T(1/9)) where {T<:Number}
    # Theoretical minimum of the wetting potential, two dimensional on array
    
    # Actual formular of the disjoining potential, long range attracion short range repulsion.
    Π_h = (γ .* (1 .- cospi.(θ)) 
          .* (exponents[1] - 1)*(exponents[2] - 1) / ((exponents[1] - exponents[2])*h_star) 
          .* ((h_star ./ height).^exponents[1] .- (h_star ./ height).^exponents[2]))
    
    return Π_h 
end
# One dimensional implementation
function Π(mom::JuSwalbe.Macroquant{Vector{T}, Vector{T}}; h_star::T=T(0.1), exponents=[9,3], γ::T=T(0.01), θ::Vector{T}=ones(T,1)*T(1/9)) where {T<:Number}   
    # Theoretical minimum of the wetting pontential
    
    # Actual formular of the disjoining potential, long range attracion short range repulsion.
    Π_h = (γ .* (1 .- cospi.(θ)) 
          .* (exponents[1] - 1)*(exponents[2] - 1)/((exponents[1] - exponents[2])*h_star) 
          .* ((h_star ./ mom.height).^exponents[1] .- (h_star ./ mom.height).^exponents[2]))
    
    return Π_h 
end

function Π(height::Array{T,1}; h_star::T=T(0.1), exponents=[9,3], γ::T=T(0.01), θ::Vector{T}=ones(T,1)*T(1/9)) where {T<:Number}
    # Theoretical minimum of the wetting potential, one dimensional on array
    
    # Actual formular of the disjoining potential, long range attracion short range repulsion.
    Π_h = (γ .* (1 .- cospi.(θ)) 
          .* (exponents[1] - 1)*(exponents[2] - 1) / ((exponents[1] - exponents[2])*h_star) 
          .* ((h_star ./ height).^exponents[1] .- (h_star ./ height).^exponents[2]))
    
    return Π_h 
end

"""
    Δh(mom::JuSwalbe.Macroquant)

Calculates the laplacian of the height field with periodic boundaries.

The calculation of the laplacian is central for the thin film evolution.
Only the pressure gradient will induce a flow, at least for the non-fluctuating version.
Therefore it is fairly important to have an accurate computation of the laplacian.

# Math
The Laplace equation is given by 

`` \\Delta \\rho = 0 ``

The laplace operator is simply Δ = ∇ . ∇. 
Since the film pressure has a contribution Δh we need a discreticed laplace operator.

## Two spatial dimensions

In two dimensions I follow the paper from Santosh and Succi.
A nine-point stencil is used and the equation goes as follow

`` \\Delta h_{i,j} = \\frac{1}{6}[4\\sum_{nn}h_{i,j} + sum_{diag}h_{i,j} - 20 h_{i,j}] ``

Where `i` and `j` are the x and y coordinates.
The index `nn` relates to the nearest neighbors, all four elements which are exactly Δx away from (i,j).
On the other hand the diagonal elements are those four which are have a distance √2Δx to (i,j). 
Periodicity is taken care of by a circular padding of the height array. 

## One spatial dimension

For a single spatial dimension it is mostly sufficient to use the central difference approach.
This approach is given by

`` \\partial_x^2 h_{i} = h_{i-1} - 2h_{i} + h_{i+1}

with `i` being the index along the spatial dimension.

# Example
```jldoctest
julia> using JuSwalbe

julia> height = reshape(collect(1.0:16.0), (4,4))
4×4 Array{Float64,2}:
 1.0  5.0   9.0  13.0
 2.0  6.0  10.0  14.0
 3.0  7.0  11.0  15.0
 4.0  8.0  12.0  16.0

julia> velx = zeros(Float64, (4,4))
4×4 Array{Float64,2}:
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0

julia> vely = zeros(Float64, (4,4))
4×4 Array{Float64,2}:
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0
 
julia> vel = JuSwalbe.Twovector(velx, vely)
JuSwalbe.Twovector{Array{Float64,2}}([0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0])

julia> moment = JuSwalbe.Macroquant(height, vel, zeros(4,4), zeros(4,4))
JuSwalbe.Macroquant{Array{Float64,2},JuSwalbe.Twovector{Array{Float64,2}}}([1.0 5.0 9.0 13.0; 2.0 6.0 10.0 14.0; 3.0 7.0 11.0 15.0; 4.0 8.0 12.0 16.0], JuSwalbe.Twovector{Array{Float64,2}}([0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0]), [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0], [0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0])

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

julia> height = collect(1.0:16.0)
16-element Array{Float64,1}:
  1.0
  2.0
  3.0
  4.0
  5.0
  6.0
  7.0
  8.0
  9.0
 10.0
 11.0
 12.0
 13.0
 14.0
 15.0
 16.0

julia> Δh(height)
16-element Array{Float64,1}:
  16.0
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
   0.0
   0.0
   0.0
   0.0
 -16.0

```

# References
## Two spatial dimensions
There are plenty of papers concerning discret differentail operators.
Many of them are good, although definitly not the best I go here with:
- [Isotropic discrete Laplacian operators from lattice hydrodynamics](https://www.sciencedirect.com/science/article/pii/S0021999112004226)

## One spatial dimension
Almost any reference on discrete differentiation is good enough.

"""
function Δh(mom::JuSwalbe.Macroquant{Matrix{T}, Twovector{Matrix{T}}}) where {T<:Number}
    # Get the size of the problem
    width, thick = size(mom.height)
    # Build a wrap around to mimic periodic boundaries
    wrapper = 1
    toruswrapped = padarray(mom.height, Pad(:circular, wrapper, wrapper))
    # Generate a dummy array for the solution
    Δ = zeros(T, size(toruswrapped))
    # Exclude the wrapped boundaries
    for j in 1:width, i in 1:thick
        @inbounds Δ[i,j] = 1.0/6.0 * (4.0 * (toruswrapped[i+1, j] + toruswrapped[i-1, j] + toruswrapped[i, j+1] + toruswrapped[i, j-1])
                                     + 1.0 * (toruswrapped[i+1, j+1] + toruswrapped[i+1, j-1] + toruswrapped[i-1, j+1] + toruswrapped[i-1, j-1])
                                     - 20.0 * toruswrapped[i,j])
    end
    # Return the result without the wrapped boundaries
    return Δ[1:(end-2), 1:(end-2)]
end

function Δh(height::Array{T,2}) where {T<:Number}
    # Get the size of the problem
    width, thick = size(height)
    # Build a wrap around to mimic periodic boundaries
    wrapper = 1
    toruswrapped = padarray(height, Pad(:circular, wrapper, wrapper))
    # Generate a dummy array for the solution
    Δ = zeros(T, size(toruswrapped))
    # Exclude the wrapped boundaries
    for j in 1:width, i in 1:thick
        @inbounds Δ[i,j] = 1.0/6.0 * (4.0 * (toruswrapped[i+1, j] + toruswrapped[i-1, j] + toruswrapped[i, j+1] + toruswrapped[i, j-1])
                                     + 1.0 * (toruswrapped[i+1, j+1] + toruswrapped[i+1, j-1] + toruswrapped[i-1, j+1] + toruswrapped[i-1, j-1])
                                     - 20.0 * toruswrapped[i,j])
    end
    # Return the result without the wrapped boundaries
    return Δ[1:(end-2), 1:(end-2)]
end

# The one dimensonal cases
function Δh(mom::JuSwalbe.Macroquant{Vector{T}, Vector{T}}) where {T<:Number}
    # Get the size of the problem
    len = length(mom.height)
    # Build a wrap around to mimic periodic boundaries
    wrapper = 1
    # Ensures stability of height
    copyheight = deepcopy(mom.height)
    # Extend the one dimensional array periodic
    pushfirst!(copyheight, last(copyheight))
    push!(copyheight, copyheight[2])
    # Generate a dummy array for the solution
    Δ = zeros(T, size(copyheight))
    # Exclude the wrapped boundaries
    for i in 2:len+1
        @inbounds Δ[i] = copyheight[i-1] - 2*copyheight[i] + copyheight[i+1] 
    end
    # Return the result without the wrapped boundaries
    return Δ[2:(end-1)]
end

function Δh(height::Array{T,1}) where {T<:Number}
    # Get the size of the problem
    len = length(height)
    # Build a wrap around to mimic periodic boundaries
    wrapper = 1
    # Ensures stability of height
    copyheight = deepcopy(height)
    # Extend the one dimensional array periodic
    pushfirst!(copyheight, last(copyheight))
    push!(copyheight, copyheight[2])
    # Generate a dummy array for the solution
    Δ = zeros(T, size(copyheight))
    # Exclude the wrapped boundaries
    for i in 2:len+1
        @inbounds Δ[i] = copyheight[i-1] - 2*copyheight[i] + copyheight[i+1] 
    end
    # Return the result without the wrapped boundaries
    return Δ[2:(end-1)]
end