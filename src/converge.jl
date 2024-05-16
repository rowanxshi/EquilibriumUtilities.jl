# MAIN FUNCTIONALITY
"""
	ConvergeParameters{T1 <: Real, T2 <: Integer, T3}

Contains optional parameters for the [`converge`](@ref) solver.

See also [`converge`](@ref).

Fields
===
* `diff_tol::T1 = 1e-6`: If the value returned by `step_diff` is less than this tolerance, stop.
* `up_tol::T1 = zero(diff_tol)`: If the value returned by `update` is less than this tolerance, stop.
* `max_iter::T2 = 200`: Maximum number of iterations before giving up.
* `msg::T3 = "No convergence"`: Warning message to display if there isn't convergence within the maximum number of iterations.
* `verbose::Bool = false`: Print the `diff` every iteration?
"""
@kwdef struct ConvergeParameters{T1 <: Real, T2 <: Integer, T3}
	diff_tol::T1 = 1e-3
	up_tol::T1 = zero(diff_tol)
	max_iter::T2 = 200
	msg::T3 = "No convergence"
	verbose::Bool = false
end
"""
	converge(update::Function, step_diff::Function; kw...) -> (converged, stalled_update)

Iterate until convergence. In particular, repeatedly:
* apply `step_diff()` which should do an iteration step, then return the deviation
* apply `update()`, which should return the size of the update
until convergence. The convergence criterion is controlled by keyword arguments; see [`ConvergeParameters`](@ref).

The function returns a pair of booleans, the first of which indicates whether convergence was reached; the second of which indicates whether iteration was aborted because the update was too small.

For an example, see this package's tests which uses this function to solve the finite-firm CES game.

See also [`update!`](@ref).
"""
function converge(update::Function, step_diff::Function; kw...)
	converge(update, step_diff, ConvergeParameters(; kw...))
end
"""
	converge(update::Function, step_diff::Function, p::ConvergeParameters)
"""
function converge(update::A, step_diff::B, p::ConvergeParameters) where {A, B}
	(; diff_tol, up_tol, max_iter, msg, verbose) = p
	diff = one(diff_tol) + diff_tol;
	small_up = false
	conv = false
	
	if verbose && (max_iter > zero(max_iter))
		meter = ProgressMeter.ProgressThresh(diff_tol; color=:blue, output=stdout, showspeed=true)
	end
	for iter in 1:max_iter
		diff = step_diff()
		conv = (diff < diff_tol)
		conv && break
		up = update()
		if (up < up_tol) && small_up
			break # break if the update is too small twice in a row
		end
		verbose && ProgressMeter.update!(meter, diff; showvalues = [(:iter, iter), (:up, up)])
		(up < up_tol) && (small_up = true)
	end
	!conv && @warn msg
	conv, small_up
end

"""
Optional helper function to use when writing an `update()` function for [`converge`](@ref).

Update `main` according to `secondary`, with a `dampen`ing factor. Useful for iterative algorithms. Once complete, `main` will hold the updated value and `secondary` will hold main's original value (to keep a record of previous iteration).

Returns the `norm` of the update step (useful for breaking iteration if the step size is too small).

See also [`converge`](@ref).
"""
function update! end
"""
	update!(main, secondary; norm = infnorm_pctdev, dampen = 0.5)

See [`infnorm_pctdev`](@ref).
"""
function update!(main, secondary; dampen = 0.5, norm = infnorm_pctdev)
	for (i_main, i_sec) in zip(eachindex(main), eachindex(secondary))
		main[i_main], secondary[i_sec] = secondary[i_sec], main[i_main]
		main[i_main] += dampen*(secondary[i_sec] - main[i_main])
	end
	norm(secondary, main)
end
"""
	infnorm_pctdev(v1, v2)

Calculates the difference between two vectors by finding the maximum absolute percentage deviation.
"""
function infnorm_pctdev(v1, v2)
	maximum(abs_pct_dev, zip(v1, v2))
end
abs_pct_dev((n1, n2)) = abs(n2 - n1)/zero_safe(n1)

function bisection(f::Function, interval; converge_kw...)
	diff!() = let f = f, interval = interval
		f(sum(interval)/2)
	end
	update!() = let f = f, interval = interval
		middle = f(sum(interval)/2)
		left = f(first(interval))
		right = f(last(interval))
		new = sum(interval)/2
		if sign(middle) == sign(left)
			old = interval[first(eachindex(interval))]
			interval[first(eachindex(interval))] = new
			abs(old - new)
		else
			old = interval[last(eachindex(interval))]
			interval[last(eachindex(interval))] = new
			abs(old - new)
		end
	end
	!(length(interval) == 2) && error("provided interval doesn't have two points")
	sort!(interval)
	converge(update!, diff!; converge_kw...)
end
