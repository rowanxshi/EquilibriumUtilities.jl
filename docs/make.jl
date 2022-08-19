push!(LOAD_PATH,"../src/")
using Documenter, EquilibriumUtilities

makedocs(sitename="Equilibrium Utilities")

deploydocs( repo = "github.com/rowanxshi/EquilibriumUtilities.jl.git")