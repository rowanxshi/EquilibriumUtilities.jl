using Documenter, DocumenterInterLinks, EquilibriumUtilities

links = InterLinks("Base" => "https://docs.julialang.org/en/v1/")
makedocs(sitename="Equilibrium Utilities", plugins=[links; ])

deploydocs( repo = "github.com/rowanxshi/EquilibriumUtilities.jl.git")
