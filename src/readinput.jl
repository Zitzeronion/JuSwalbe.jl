"""
    inputfile

Abstract type for all kinds of input files
"""
abstract type inputfile end

"""
    Inputconstants = new(lx, ly, maxruntime, dumping, gravity, γ, δ)

Struct containing input parameters.

Contains `.lx` lattice points in x-direction, `.ly` lattice points in y-direction. 
Other fields are `.maxruntime` for the maximal number of time steps and `.dumping` to limit the number of output files.
On top of these there are physical quantities such as `.gravity`, `.γ` and `.δ` 
for the values of gravitational acceleration, fluids surface tension and the slip length.
The example relates to an quadratic lattice 20 times 20 lattice units in area. 
Run for 100 lattice Boltzmann time steps only printing output every 10 time steps.
Having no gravity and a surface tension of 0.01 and a slip length of 1.  

# Example
```jldoctest firsttest
julia> using JuSwalbe

julia> new_input = JuSwalbe.Inputconstants()
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

julia> new_input.γ
0.01
```

# References
See also: [`readinput`](@ref), [`findargument`](@ref), [`computeslip`](@ref)
"""
@with_kw struct Inputconstants <: inputfile
    lx = 512
    ly = 512
    maxruntime = 100000
    dumping = 1000
    τ = 1.0
    gravity = 0.0
    γ = 0.01
    δ = 1.0
    μ = 1 / 3 * (2 - τ) / 2 * τ
end

"""
    readinput(file)

Reads input parameters from a `file`.

The expected amount of parameters can be addressed with [`Inputconstants`](@ref).
For now it expects seven values for different runtime constants.

# Example
```jldoctest secondtest
julia> using JuSwalbe, DelimitedFiles

julia> args = ["Lattice_points_x" 10; "Lattice_points_y" 5; "Max_run_time" 1000; "Output_dump" 100; "Relaxation_rate" 1.0; "gravity" 0.0; "surface_tension" 0.01; "slippage" 1.0] # Generate a text file with input
8×2 Array{Any,2}:
 "Lattice_points_x"    10
 "Lattice_points_y"     5
 "Max_run_time"      1000
 "Output_dump"        100
 "Relaxation_rate"      1.0
 "gravity"              0.0
 "surface_tension"      0.01
 "slippage"             1.0

julia> writedlm("test.txt", args)

julia> test = readinput("test.txt")
JuSwalbe.Inputconstants
  lx: Int64 10
  ly: Int64 5
  maxruntime: Int64 1000
  dumping: Int64 100
  τ: Float64 1.0
  gravity: Float64 0.0
  γ: Float64 0.01
  δ: Float64 1.0
  μ: Float64 0.16666666666666666

julia> test.lx
10

julia> test.γ
0.01

julia> test.γ + test.δ
1.01

julia> isa(test.lx + test.gravity, Int32)
false

julia> rm("test.txt")
```
"""
function readinput(file)
    # Length of the input file
    num = countlines(file)
    # Actually reading of the file and saving to an array
    input = readdlm(file)
    # Extracting of numerical values
    values = []
    arguments = ["Lattice_points_x", 
                "Lattice_points_y", 
                "Max_run_time", 
                "Output_dump", 
                "Relaxation_rate",
                "gravity", 
                "surface_tension", 
                "slippage",
                "viscosity"]

    for i in arguments
        push!(values, findfirst(x -> x == i, vec(input)))
    end
    
    lx = findargument(input, "Lattice_points_x")
    ly = findargument(input, "Lattice_points_y")
    maxruntime = findargument(input, "Max_run_time")
    dumping = findargument(input, "Output_dump")
    τ = findargument(input, "Relaxation_rate")
    gravity = findargument(input, "gravity")
    γ = findargument(input, "surface_tension")
    δ = findargument(input, "slippage")

    runtimeconstants = Inputconstants(lx=lx, ly=ly, maxruntime=maxruntime, dumping=dumping, τ=τ, gravity=gravity, γ=γ, δ=δ)
    
    return runtimeconstants
end

"""
    findargument(arr, str)

Searches for a numerical value based on a str input and returns the value.

# Example
```jldoctest
julia> using JuSwalbe

julia> arr = ["hmm" 1; "yeah" 0.01; "world" 1090]
3×2 Array{Any,2}:
 "hmm"       1
 "yeah"      0.01
 "world"  1090

julia> world = findargument(arr, "world")
1090

```

# References
See also: [`readinput`](@ref)
"""
function findargument(arr::Array{Any,2}, argument::String)
    for (i,arg) in enumerate(arr[:,1])
        if arr[i,1] == argument
            value = arr[i,2]
            return value
        end
    end
end
