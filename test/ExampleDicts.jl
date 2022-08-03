d = Dict(:w => rand(10), :L => rand(10))

module ExampleDicts
using Main.EquilibriumUtilities

# set up a custom struct that's a parameterised Pollak demand function
struct Pollak{T <: Real} <: WrappedDict{T}
	internal_dict::Dict{Symbol, T}
	function Pollak{T}(d) where T
		all(in(keys(d)), (:σ, :γ)) || @error "Missing necessary Pollak parameters"
		d[:γ] ≤ zero(T) || @error "Provided γ must be non-positive" 
		new{T}(d)
	end
end

(p::Pollak)(x) = -p[:γ]*x^(-p[:σ]) + p[:γ]

function CES!(p::Pollak)
	p[:γ] = 0
	p
end

p = Pollak{Float64}(Dict(:σ => 3., :γ => -1. ))

# ... functions dispatching on ::Pollak

# set up a custom struct that's a description of several labour markets (e.g. in a trade model)
struct Labour{T <: Vector{<: Real}} <: WrappedDict{T}
	internal_dict::Dict{Symbol, T}
	function Labour{T}(d::Dict{Symbol, T}) where T
		# choose wage in the first country as numeraire
		normalise!(d[:w])
		# normalise total labour to sum to one
		normalise!(d[:L], sum(d[:L]))
		new{T}(d)
	end
end
Labour(d::Dict{Symbol, T}) where {T <: Vector{<: Real}} = Labour{T}(d)
europe = Labour(Main.d)
Base.show(io::IO, ::MIME"text/plain", l::Labour) = print(io, EquilibriumUtilities.pretty(l))

# ... functions dispatching on ::Labour

end