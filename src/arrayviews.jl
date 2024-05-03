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

