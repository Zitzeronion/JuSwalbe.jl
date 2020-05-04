using Documenter, JuSwalbe

makedocs(
    sitename="JuSwalbe.jl",
    modules=[JuSwalbe])

deploydocs(
    repo = "github.com/Zitzeronion/JuSwalbe.git",
    traget="build"
)