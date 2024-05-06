"""
	chunk(V::AbstractVector, N::Integer)

Chunk a vector `V` into [`Base.view`](@extref)s each of length `N`. Returns a `Tuple` of [`Base.view`](@extref)s.
"""
function chunk(V::AbstractVector, N::Integer)
	chunked_indices = Iterators.partition(eachindex(V), N)
	(view(V, chunk) for chunk in chunked_indices)
end
"""
	diag_view(mat::AbstractMatrix)

Returns a view of the matrix diagonal, assuming one-based indexing.

See also [`offdiag_view`](@ref), [`issquare`](@ref).
"""
function diag_view(mat::AbstractMatrix)
	!issquare(mat) && error("matrix is not square: $mat")
	S = first(size(mat))
	@view mat[1:(S+1):end]
end
"""
	offdiag_view(mat::AbstractMatrix)

Returns a view of the matrix off-diagonals, assuming one-based indexing.

See also [`diag_view`](@ref), [`issquare`](@ref).
"""
function offdiag_view(mat::AbstractMatrix)
	!issquare(mat) && error("matrix is not square: $mat")
	N = first(size(mat))
	isdiagind = in(1:(N+1):(N^2))
	offdiaginds = .!(isdiagind.(eachindex(mat)))
	@views mat[offdiaginds]
end
"""
	issquare(mat::AbstractMatrix)

Is the matrix square?

See also [`diag_view`](@ref), [`offdiag_view`](@ref).
"""
function issquare(mat::AbstractMatrix)
	dims = size(mat)
	isequal(dims...) 
end

