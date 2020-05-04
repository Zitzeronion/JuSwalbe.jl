using Documenter
using JuSwalbe
doctest(JuSwalbe; manual = false)

makedocs(
    sitename="JuSwalbe.jl",
    modules=[JuSwalbe],
    pages=["index.md", "manual.md", "reference.md", "devnotes.md"])

deploydocs(
    repo = "github.com/Zitzeronion/JuSwalbe.git",
)