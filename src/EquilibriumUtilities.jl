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

# ARRAY VIEWS
"""
	chunk(V::AbstractVector, N::Integer)

Chunk a vector `V` into [`Base.view`](@extref)s each of length `N`. Returns a `Tuple` of [`Base.view`](@extref)s.
"""
function chunk(V::AbstractVector, N::Integer)
	chunked_indices = Iterators.partition(eachindex(V), N)
	(view(V, chunk) for chunk in chunked_indices)
end
function diag_view(mat::AbstractMatrix)
	!issquare(mat) && error("matrix is not square: $mat")
	S = first(size(mat))

	@view mat[1:(S+1):end]
end
function offdiag_view(mat::AbstractMatrix)
	!issquare(mat) && error("matrix is not square: $mat")
	N = first(size(mat))
	isdiagind = in(1:(N+1):(N^2))
	offdiaginds = .!(isdiagind.(eachindex(mat)))
	@views mat[offdiaginds]
end
function issquare(mat::AbstractMatrix)
	dims = size(mat)
	isequal(dims...) 
end

include("WrappedDict.jl")
include("ConvergenceState.jl")
include("solvers.jl")
include("prettyprinting.jl")

export newton, converge, infnorm_pctdev, normalise!, zero_safe, chunk, issquare, quietly, diag_view, dampen, update!, WrappedDict

end
