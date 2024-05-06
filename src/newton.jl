"""
	NewtonParameters

Contains optional parameters for the [`newton`](@ref) solver.

See also [`newton`](@ref).

Fields
===
* `verbose::Bool = false`: Print some info every iteration?
* `step_tol::Float64 = 1e-8`: If the (absolute value of) step size is smaller than this tolerance, stop.
* `f_tol::Float64 = step_tol`: If the (absolute value of) function value is smaller than this tolerance, stop.
* `max_iter::Integer = 750`: Maximum number of iterations before giving up.
* `msg = "No Newton convergence"`: Warning message to display if there isn't convergence within the maximum number of iterations.
* `l::Real = -Inf`: Left bound.
* `r::Real = Inf`: Right bound.
"""
@kwdef struct NewtonParameters{T1 <: Integer, T2, T3 <: Real, T4 <: Real}
	verbose::Bool = false
	step_tol::Float64 = 1e-8
	f_tol::Float64 = step_tol
	max_iter::T1 = 750
	msg::T2 = "No Newton convergence"
	l::T3 = -Inf
	r::T4 = Inf
end
const default_newton_parameters = NewtonParameters()
"""
Simple Newton method implementation.

See also [`NewtonParameters`](@ref) for more keyword argument options.
"""
function newton end

"""
	newton(f_f′, x; kwargs...)
	
`f_f′` should be a function which returns the tuple `(f, f′)` of function value and derivative. `x` is the initial guess.
"""
function newton(f_f′::Function, x; kw...)
	newton(f_f′, x, NewtonParameters(; kw...))
end
function newton(f_f′::Function, x, p::NewtonParameters)
	(; step_tol, f_tol, max_iter, verbose, r, l, msg) = p
	diff = 2step_tol
	f = 2f_tol
	iter = 1
	while iter ≤ max_iter
		iter += 1
		f, f′ = f_f′(x)
		abs(f) < f_tol && break
		diff = -f/f′
		diff = isinf(diff) ? one(diff) : diff
		abs(diff) < step_tol && break
		# @debug Printf.@sprintf "x: %.3e, f: %.3e, diff: %.3e" x f diff
		verbose && println("x: $x, f: $f, diff: $diff")
		x += if (x+diff ≥ r)
			isequal(x, r) && break # give up if even at upper bound, newton wants to go higher
			@debug "newton guess too high"
			r - x # set to upper if newton wants to exceed
		elseif (x+diff ≤ l)
			isequal(x,l) && break
			@debug "newton guess too low"
			l - x
		else
			diff
		end
	end
	if (abs(f) > f_tol) && (abs(diff) > step_tol)
		@error msg*"\n f = $f"
	end
	x
end;
"""
	newton(f, f′, x; kwargs...)
	
`f` should be a function which returns the function value and `f′` its derivative. `x` is the initial guess.
"""
function newton(f::Function, f′::Function, x; newton_params...)
	newton(f, f′, x, NewtonParameters(; newton_params...))
end
function newton(f::Function, f′::Function, x, p::NewtonParameters)
	f_f′(x) = let f = f, f′ = f′
		(f(x), f′(x))
	end
	newton(f_f′, x, p)
end
