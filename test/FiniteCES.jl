module FiniteCES
using Main.EquilibriumUtilities

const σ = 5 # CES parameter
const c = collect(1:0.1:20) # the marginal costs of the participants
const s_old = Vector{Float64}(undef, length(c)) # vector to hold shares guesses
const s_new = Vector{Float64}(undef, length(c)) # vector to hold implied shares

function init()
	# initial guess: even shares
	fill!(s_old, length(s_old)^(-1))
end
function step_diff()
	# using s_old, compute implied prices and store in s_new
	for firm in eachindex(s_old)
		ε = s_old[firm] + (1 - s_old[firm])*σ
		μ = ε/(ε - 1)
		s_new[firm] = μ*c[firm]^(1 - σ)
	end
	
	# divide by price index for shares
	normalise!(s_new, sum(s_new))
	
	return infnorm_pctdev(s_old, s_new)
end
function update()
	# use the average between guess and implied shares
	up_diff = infnorm_pctdev(s_new, s_old)/2
	s_old .+= s_new
	s_old ./= 2
	up_diff
end

export init, step_diff, update

end
