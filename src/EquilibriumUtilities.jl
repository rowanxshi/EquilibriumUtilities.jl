"""
A package of basic utility functions used commonly when computing economic equilibria.

See [`converge`](@ref), [`normalise!`](@ref), [`zero_safe`](@ref), [`newton`](@ref).
"""
module EquilibriumUtilities


"""
	normalise!(v, factor = zero_safe(first(v)))
	
Normalise `v` by some `factor`. Useful for normalising nominal prices (so the first element is the numeraire) or shares (with `factor = sum(v)`).

See also [`zero_safe`](@ref).
"""
function normalise!(v, factor = zero_safe(first(v)))
	v ./= factor
end
"""
	zero_safe(x)

If `x` is zero, return `one(x)`. Otherwise, return `x`. Useful for safely dividing by `x`.
"""
zero_safe(x) = iszero(x) ? one(x) : x

"""
Simple Newton method implementation.

## Keywords
* `step_tol = 1e-8` : If the (absolute value of) step size is smaller than this tolerance, stop.
* `f_tol = step_tol` : If the (absolute value of) function value is smaller than this tolerance, stop.
* `max_iter = 750` : Maximum number of iterations before giving up.
* `l = -Inf` : Left bound.
* `r = Inf` : Right bound.
* `msg = "No Newton convergence"` : Warning message to display if there isn't convergence within the maximum number of iterations.
* `verbose::Bool = false` : Print some info every iteration?
"""
function newton end

"""
	newton(f_f′, x; kwargs...)
	
`f_f′` should be a function which returns the tuple `(f, f′)` of function value and derivative. `x` is the initial guess.
"""
function newton(f_f′::Function, x; verbose::Bool = false, step_tol = 1e-8, f_tol = step_tol, max_iter = 750, msg = "No Newton convergence", l = -Inf, r = Inf)
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
		@debug("x: $x, f: $f, diff: $diff")
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
function newton(f::Function, f′::Function, x::Real; newton_params...)
	f_f′(x) = let f = f, f′ = f′
		(f(x), f′(x))
	end
	newton(f_f′, x; newton_params...)
end

"""
	converge(update::Function, step_diff::Function, init::Function; kwargs...)

Iterate until convergence. In particular, the problem is initiated with `init()`. Then, repeatedly apply `step_diff()` which should do an iteration step, then return the `diff`. `converge` will compare the returned difference with `tol`; if it's smaller, then converge is reached and iteration stops. Otherwise, it will apply `update()`, then continue iterating.

## Keywords
* `tol::Real = 1e-6` : If the `diff` returned by `step_diff` is less than this tolerance, stop.
* `max_iter::Integer = 200` : Maximum number of iterations before giving up.
* `msg = "No convergence"` : Warning message to display if there isn't convergence within the maximum number of iterations.
* `history = nothing` : Container in which to store the iteration history of `diff`s. Useful to check speed of convergence. If `nothing` is provided, saves no history.
* `verbose::Bool = false` : Print the `diff` every iteration?

## Example
Consider the finite-firm CES model. We iteratively find market shares among the participants.
```julia
const σ = 5 # CES parameter
const c = collect(1:10) # the marginal costs of the participants
const s_old = Vector{Float64}(undef, length(z)) # vector to hold shares guesses
const s_new = Vector{Float64}(undef, length(z)) # vector to hold implied shares

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
	
	return v_diff(s_old, s_new)
end
function update()
	# use the average between guess and implied shares
	s_old .+= s_new
	s_old ./= 2
end

converge(update, step_diff, init)
```
"""
function converge(update::Function, step_diff::Function, init::Function; history = nothing, tol::Real = 1e-6, max_iter::Integer = 200, msg = "No convergence", verbose::Bool = false)
	init()
	diff = one(tol) + tol; 
	for iter in 1:max_iter
		diff = step_diff()
		diff < tol && break
		verbose && @info diff
		!isnothing(history) && push!(history, diff)
		update()
	end
	diff < tol || @warn msg
end

"""
	dampen(history; kwargs...)

Given an iteration `history`, return a dampening factor. At the moment, just uses the last value of `history`.

## Keywords
* `slow = 0.95` : the dampening factor for slow updating (when `diff` is above 10).
* `med = 0.75` : the dampening factor for medium updating (when `diff` is above 1).
* `fast = 0.5` : the dampening factor for fast updating.
"""
function dampen(history; slow = 0.95, med = 0.75, fast = 0.5)
	(isempty(history) || last(history) > 10) && return slow
	last(history) > 1 && return med
	return fast
end

"""
	v_diff(v1, v2)

Calculate the distance between two vectors as the sum of element-wise absolute difference.
"""
function v_diff(v1, v2)
	sum(zip(v1, v2)) do (x1, x2)
		abs(x1 - x2)
	end
end

export newton, converge, v_diff, normalise!, zero_safe, dampen

end
