abstract type inputfile end

struct inputconstants <: inputfile
    lx::Int32
    ly::Int32
    maxruntime::Int64
    dumping::Int32
    gravity::Float64
    γ::Float64
    δ::Float64
end

function readinput(file)
    # Length of the input file
    num = countlines(file)
    # Actually reading of the file and saving to an array
    input = readdlm(file)
    # Extracting of numerical values
    values = []
    arguments = ["Lattice_points_x", "Lattice_points_y", "Max_run_time", "Output_dump", "gravity", "surface_tension", "slippage"]
    for arg in arguments
        push!(values, findfirst(x -> x == arg, vec(input))) 
    end
    runtimeconstants = inputconstants(input[values[1], 2], 
                                      input[values[2], 2], 
                                      input[values[3], 2],
                                      input[values[4], 2],
                                      input[values[5], 2],
                                      input[values[6], 2],
                                      input[values[7], 2])
    
    return runtimeconstants
end
