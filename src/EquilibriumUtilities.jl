"""
A package of basic utility functions used commonly when computing economic equilibria.

See [`converge`](@ref), [`newton`](@ref), [`normalise!`](@ref), [`zero_safe`](@ref).
"""
module EquilibriumUtilities
import ProgressMeter
import Logging
import Printf

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
function quietly(f::Function)
	out = Logging.with_logger(Logging.NullLogger()) do
		f()
	end
	out
end

include("arrayviews.jl")
include("WrappedDict.jl")
include("dynamicdampen.jl")
include("newton.jl")
include("converge.jl")
include("prettyprinting.jl")

export newton, converge, infnorm_pctdev, normalise!, zero_safe, chunk, issquare, quietly, diag_view, update!, WrappedDict

end
