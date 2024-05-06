"""
	abstract type WrappedDict{T} <: AbstractDict{Symbol, T}

Eases creating custom `struct`s that are basically dictionaries. (Essentially implements `keys`, `values`, `length`, `iterate`, `getindex`, `setindex!`, `get`, and `get!` by delegating to the necessary field `internal_dict`.)

A concrete subtype of `WrappedDict` _must_ have a field `internal_dict::Dict{Symbol, T}`, which is the dictionary.

For an example, see this package's tests with an implementation of a Pollak demand function `struct` as a `WrappedDict`.

See also [`validate`](@ref).
"""
abstract type WrappedDict{T} <: AbstractDict{Symbol, T} end

Base.keys(s::WrappedDict) = keys(getfield(s, :internal_dict))
Base.values(s::WrappedDict) = values(getfield(s, :internal_dict))
Base.length(s::WrappedDict) = length(getfield(s, :internal_dict))
Base.iterate(s::WrappedDict, state...) = iterate(getfield(s, :internal_dict), state...)
Base.getindex(s::WrappedDict, k::Symbol) = getindex(getfield(s, :internal_dict), k)
Base.setindex!(s::WrappedDict, value, k::Symbol) = setindex!(getfield(s, :internal_dict), value, k)
Base.get(s::WrappedDict, k::Symbol, default) = get(getfield(s, :internal_dict), k, default)
Base.get!(s::WrappedDict, k::Symbol, default) = get!(getfield(s, :internal_dict), k, default)
Base.getproperty(s::WrappedDict, k::Symbol) = getindex(s, k)
Base.setproperty!(s::WrappedDict, k::Symbol, value) = setindex!(s, value, k)
Base.delete!(s::WrappedDict, k::Symbol) = delete!(getfield(s, :internal_dict), k)
Base.propertynames(s::WrappedDict) = keys(getfield(s, :internal_dict))

"""
	validate(given, needed)	

For each key in `needed`, check that it's present in `given`.
"""
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

