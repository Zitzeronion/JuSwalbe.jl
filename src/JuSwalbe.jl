module JuSwalbe

using DelimitedFiles, JSON

include("call.jl")
include("macroscopic_typs.jl")
include("distribution_types.jl")
include("readinput.jl")
include("equilibriumcalculation.jl")

export call, readinput, calc_equilibrium_distribution

end # module
