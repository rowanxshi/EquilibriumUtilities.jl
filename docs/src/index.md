# EquilibriumUtilities.jl Documentation
A small package of utilities commonly used to compute economic models.

Exports
```@example
using EquilibriumUtilities # hide
join(names(EquilibriumUtilities), ", ") # hide
```

## Iterative convergence
The generalised iterative routine and its parametres:
```@autodocs
Modules = [EquilibriumUtilities]
Pages = ["converge.jl"]
```

### Reactive dampening
**Optional** helpers, mostly implementing a dampening factor that dynamically adjusts as the equilibrium progresses.

The convergence routine above need not use this routine: the user can always define the `update()` function however they would like. However, this optional functionality allows the solver to choose a dampening factor taking into account the convergence path. It will be provided in the `ds::DampenState` `struct` so the user can access it when defining their `update()` function.
```@autodocs
Modules = [EquilibriumUtilities]
Pages = ["dynamicdampen.jl"]
```

## Newton method
```@autodocs
Modules = [EquilibriumUtilities]
Pages = ["newton.jl"]
```

## Small utilities
```@autodocs
Modules = [EquilibriumUtilities]
Pages = ["EquilibriumUtilities.jl", "arrayviews.jl", "prettyprinting.jl"]
```

## `WrappedDict`
```@autodocs
Modules = [EquilibriumUtilities]
Pages = ["WrappedDict.jl"]
```

