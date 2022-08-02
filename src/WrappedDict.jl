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

"""
	abstract type WrappedDict{T} <: AbstractDict{Symbol, T}

Eases creating custom `struct`s that are basically dictionaries.[^1] A concrete subtype of `WrappedDict` _must_ have a field `internal_dict::Dict{Symbol, T}`, which is the dictionary. 

[^1]: Essentially implements `keys`, `values`, `length`, `iterate`, `getindex`, `setindex!`, and `get` by delegating to the necessary field `internal_dict`.

## Example

```julia
# set up a custom struct that's a parameterised Pollak demand function
struct Pollak{T <: Real} <: WrappedDict{T}
	internal_dict::Dict{Symbol, T}
	function Pollak{T}(d) where T
		all(in(keys(d)), (:σ, :γ)) || @error "Missing necessary Pollak parameters"
		d[:γ] ≤ zero(T) || @error "Provided γ must be non-positive" 
		new{T}(d)
	end
end

(p::Pollak)(x) = p[:γ]*x^(-p[:σ]) - p[:γ]

function CES!(p::Pollak)
	p[:γ] = 0
	p
end

# ... functions dispatching on ::Pollak
```

"""
abstract type WrappedDict{T} <: AbstractDict{Symbol, T} end

Base.keys(s::WrappedDict) = keys(s.internal_dict)
Base.values(s::WrappedDict) = values(s.internal_dict)
Base.length(s::WrappedDict) = length(s.internal_dict)
Base.iterate(s::WrappedDict, state...) = iterate(s.internal_dict, state...)
Base.getindex(s::WrappedDict, k::Symbol) = getindex(s.internal_dict, k)
Base.setindex!(s::WrappedDict, value, k::Symbol) = setindex!(s.internal_dict, value, k)
function Base.get(s::WrappedDict, k::Symbol, default...) 
	k ∈ keys(s) ? getindex(s, k) : begin
		@warn "no field $k; returning the first value"
		first(values(s))
	end
end
