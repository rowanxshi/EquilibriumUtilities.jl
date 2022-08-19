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

See also [`validate`](@ref), [`pretty`](@ref).

[^1]: Essentially implements `keys`, `values`, `length`, `iterate`, `getindex`, `setindex!`, `get`, and `get!` by delegating to the necessary field `internal_dict`.
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

function pretty(s::WrappedDict, var_names = keys(s); pad = 8, digits = 4, spacer = 2)
	str = string_header(string.(var_names); pad, spacer)*"\n"
	for line in string_tablify(s, var_names; pad, digits)
		str *= line*"\n"
	end
	str
end
spaced_join(x; pad = 8) = join(rpad.(x, pad), "")
function string_header(var_names; pad = 8, spacer = 2)
	header = spaced_join(var_names; pad)
	underline = repeat("=", pad - spacer)*repeat(" ", spacer)
	header*"\n"*repeat(underline, length(var_names))
end
function string_tablify(s::WrappedDict, var_names = keys(s); pad = 8, digits = 4)
	columns = (getproperty(s, v) for v in var_names)
	(spaced_join(round.(line; digits); pad) for line in zip(columns...))
end
