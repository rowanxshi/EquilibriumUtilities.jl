module PollakStruct
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

end