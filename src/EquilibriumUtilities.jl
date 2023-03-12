"""
A package of basic utility functions used commonly when computing economic equilibria.

See [`converge`](@ref), [`normalise!`](@ref), [`zero_safe`](@ref), [`newton`](@ref).
"""
module EquilibriumUtilities


"""
	normalise!(v, factor = zero_safe(first(v)))
	
Normalise `v` by some `factor`. Useful for normalising nominal prices (so the first element is the numeraire) or shares (with `factor = sum(v)`).

See also [`zero_safe`](@ref).
"""
function normalise!(v, factor = zero_safe(first(v)))
	v ./= factor
end
"""
	zero_safe(x)

If `x` is zero, return `one(x)`. Otherwise, return `x`. Useful for safely dividing by `x`.
"""
zero_safe(x) = iszero(x) ? one(x) : x

include("ConvergenceState.jl")
include("solvers.jl")
include("WrappedDict.jl")

export newton, converge, v_diff, normalise!, zero_safe, dampen, update!, WrappedDict

end
