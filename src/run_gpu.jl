module GPU_CFD

using CUDA, Revise, Parameters, Images

include("structs.jl")
include("pressure.jl")

export params, velocity, moments
export Π, Π_cuda, Δ, kernel_laplacianperiodic, kernel_gradientperiodic


#============================================================================#
#                          End of module                                     #
#============================================================================#
end