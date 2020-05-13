"""
    inputfile

Abstract type for all kinds of input files
"""
abstract type inputfile end

"""
    inputconstants = new(lx, ly, maxruntime, dumping, gravity, γ, δ)

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

julia> new_input = JuSwalbe.inputconstants(20, 20, 100, 10, 0.0, 0.01, 1.0)
JuSwalbe.inputconstants(20, 20, 100, 10, 0.0, 0.01, 1.0)

julia> new_input.γ
0.01
```
"""
@with_kw struct Inputconstants <: inputfile
    lx::Int32 = 512
    ly::Int32 = 512
    maxruntime::Int64 =100000
    dumping::Int32 = 1000
    tau::Float64 = 1.0
    gravity::Float64 = 0.0
    γ::Float64 = 0.01
    δ::Float64 = 1.0
    μ::Float64 = 1/3*()
end

"""
    readinput(file)

Reads input parameters from a `file`.

The expected amount of parameters can be addressed with [`inputconstants`](@ref).
For now it expects seven values for different runtime constants.

# Example
```jldoctest secondtest
julia> using JuSwalbe, DelimitedFiles

julia> args = ["Lattice_points_x" 10; "Lattice_points_y" 5; "Max_run_time" 1000; "Output_dump" 100; "gravity" 0.0; "surface_tension" 0.01; "slippage" 1.0] # Generate a text file with input
7×2 Array{Any,2}:
 "Lattice_points_x"    10
 "Lattice_points_y"     5
 "Max_run_time"      1000
 "Output_dump"        100
 "gravity"              0.0
 "surface_tension"      0.01
 "slippage"             1.0

julia> writedlm("test.txt", args)

julia> test = readinput("test.txt")
JuSwalbe.inputconstants(10, 5, 1000, 100, 0.0, 0.01, 1.0)

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
                "Relaxation_rate"
                "gravity", 
                "surface_tension", 
                "slippage",
                "viscosity"]

    for arg in arguments
        push!(values, findfirst(x -> x == arg, vec(input))) 
    end
    runtimeconstants = Inputconstants(input[values[1], 2], 
                                      input[values[2], 2], 
                                      input[values[3], 2],
                                      input[values[4], 2],
                                      input[values[5], 2],
                                      input[values[6], 2],
                                      input[values[7], 2],
                                      input[values[8], 2],
                                      )
    
    return runtimeconstants
end
