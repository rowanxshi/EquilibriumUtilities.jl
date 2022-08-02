using Test
import EquilibriumUtilities

@testset "Newton" begin
	function f_f′(x)
		x^2, 2x
	end
	for x0 in (-7., 7., 8.)
		near_zero = EquilibriumUtilities.newton(f_f′, x0; step_tol = eps(x0))
		@test abs(near_zero) ≤ 1e-6
	end
end

@testset "solving finite CES" begin
	include("FiniteCES.jl")
	using Main.FiniteCES
	
	# test roughly log linear for lots of firms
	EquilibriumUtilities.converge(update, step_diff, init)
	function log_lin(s, c)
		y = log.(s)
		X = [log.(c) ones(length(y))]
		β̂ = inv(X'*X)*X'*y 
		ỹ = X*β̂
		SSR = sum((y .- ỹ).^2)/length(y)
		y̅ = sum(y)/length(y)
		VAR = sum((y .- y̅).^2)/length(y)
		SSR/VAR
	end
	@test log_lin(FiniteCES.s_old, FiniteCES.c) ≤ 1e-3

	# not log with few firms
	N = 10
	resize!.((FiniteCES.c, FiniteCES.s_new, FiniteCES.s_old), N)
	FiniteCES.c .= range(1, step = 5, length = N)
	EquilibriumUtilities.converge(update, step_diff, init)
	@test log_lin(FiniteCES.s_old, FiniteCES.c) ≥ 0.1
	
	# test symmetric case
	FiniteCES.c .= 2
	EquilibriumUtilities.converge(update, step_diff, init)
	@test all(==(inv(length(FiniteCES.c))), FiniteCES.s_old)
end

@testset "WrappedDicts" begin
	
end