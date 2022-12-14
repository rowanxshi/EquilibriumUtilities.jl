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

Iterate until convergence. In particular, the problem is initiated with `init()`. Then, repeatedly apply `step_diff()` which should do an iteration step, then return the `diff`. `converge` will compare the returned difference with `tol`; if it's smaller, then converge is reached and iteration stops. Otherwise, it will apply `update()`, which should return the size of the update.

The function returns a pair of booleans, the first of which signals whether convergence was reached; the second signals whether iteration was aborted because the update was too small.

For an example, see this package's tests which uses this function to solve the finite-firm CES game.

See also [`update!`](@ref), [`dampen`](@ref), [`v_diff`](@ref).

## Keywords
* `diff_tol::Real = 1e-6` : If the value returned by `step_diff` is less than this tolerance, stop.
* `up_tol::Real = zero(diff_tol)` : If the value returned by `update` is less than this tolerance, stop.
* `max_iter::Integer = 200` : Maximum number of iterations before giving up.
* `msg = "No convergence"` : Warning message to display if there isn't convergence within the maximum number of iterations.
* `history = nothing` : Container in which to store the iteration history of `diff`s. Useful to check speed of convergence. If `nothing` is provided, saves no history.
* `verbose::Bool = false` : Print the `diff` every iteration?
"""
function converge(update::Function, step_diff::Function, init::Function; history = nothing, diff_tol::Real = 1e-6, up_tol::Real = zero(diff_tol), max_iter::Integer = 200, msg = "No convergence", verbose::Bool = false)
	init()
	diff = one(diff_tol) + diff_tol;
	small_up = false
	conv = false
	
	for iter in 1:max_iter
		diff = step_diff()
		!isnothing(history) && push!(history, diff)
		verbose && @info diff
		conv = (diff < diff_tol)
		conv && break
		up = update()
		if (up < up_tol) && small_up
			break # break if the update is too small twice in a row
		end
		(up < up_tol) && (small_up = true)
	end
	!conv || @warn msg
	conv, small_up
end

"""
	update!(main, secondary; dampen = 0.75, v_diff = v_diff)

Update `main` according to `secondary`, with a `dampen`ing factor. Useful for iterative algorithms. Once complete, `main` will hold the updated value and `secondary` will hold main's original value (to keep a record of previous iteration). Returns the `v_diff` of the update step (useful for breaking iteration if the step size is too small).
"""
function update!(main, secondary; dampen = 0.75, v_diff = v_diff)
	for (i_main, i_sec) in zip(eachindex(main), eachindex(secondary))
		main[i_main], secondary[i_sec] = secondary[i_sec], main[i_main]
		main[i_main] += dampen*(secondary[i_sec] - main[i_main])
	end
	v_diff(main, secondary)
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
	(isnothing(history) || isempty(history) || last(history) > 10) && return slow
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