function validate(given, needed)
	msg = ""
	for f in needed
		isthere = in(f, given)
		isthere || (msg *= "$f ")
		isthere
	end
	!isempty(msg) && @warn(msg*"not provided")
	nothing
end
validate(given::Dict, needed) = validate(keys(given), needed)

"""
	abstract type WrappedDict{T} <: AbstractDict{Symbol, T}

Eases creating custom `struct`s that are basically dictionaries.[^1] A concrete subtype of `WrappedDict` _must_ have a field `internal_dict::Dict{Symbol, T}`, which is the dictionary.

For an example, see this package's tests with an implementation of a Pollak demand function `struct` as a `WrappedDict`.

[^1]: Essentially implements `keys`, `values`, `length`, `iterate`, `getindex`, `setindex!`, and `get` by delegating to the necessary field `internal_dict`.
"""
abstract type WrappedDict{T} <: AbstractDict{Symbol, T} end

Base.keys(s::WrappedDict) = keys(getfield(s, :internal_dict))
Base.values(s::WrappedDict) = values(getfield(s, :internal_dict))
Base.length(s::WrappedDict) = length(getfield(s, :internal_dict))
Base.iterate(s::WrappedDict, state...) = iterate(getfield(s, :internal_dict), state...)
Base.getindex(s::WrappedDict, k::Symbol) = getindex(getfield(s, :internal_dict), k)
Base.setindex!(s::WrappedDict, value, k::Symbol) = setindex!(getfield(s, :internal_dict), value, k)
function Base.get(s::WrappedDict, k::Symbol, default...) 
	k ∈ keys(s) ? getindex(s, k) : begin
		@warn "no field $k; returning the first value"
		first(values(s))
	end
end
Base.getproperty(s::WrappedDict, k::Symbol) = getindex(s, k)
