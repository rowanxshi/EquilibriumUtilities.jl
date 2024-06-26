"""
	NewtonParameters

Contains optional parameters for the [`newton`](@ref) solver.

See also [`newton`](@ref).

## Fields
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

"""
	ConvergeParameters

Contains optional parameters for the [`converge`](@ref) solver.

See also [`converge`](@ref).

## Fields
* `diff_tol::Real = 1e-6`: If the value returned by `step_diff` is less than this tolerance, stop.
* `up_tol::Real = zero(diff_tol)`: If the value returned by `update` is less than this tolerance, stop.
* `max_iter::Integer = 200`: Maximum number of iterations before giving up.
* `msg = "No convergence"`: Warning message to display if there isn't convergence within the maximum number of iterations.
* `history = empty!(cs.history)`: Container in which to store the iteration history of `diff`s. Useful to check speed of convergence. If `nothing` is provided, saves no history.
* `verbose::Bool = false`: Print the `diff` every iteration?
"""
@kwdef struct ConvergeParameters{T1, T2 <: Real, T3 <: Real, T4 <: Integer, T5}
	history::T1 = empty!(cs.history)
	diff_tol::T2 = 1e-6
	up_tol::T3 = zero(diff_tol)
	max_iter::T4 = 200
	msg::T5 = "No convergence"
	verbose::Bool = false
end
const default_convergence_parameters = ConvergeParameters()

"""
	converge(update::Function, step_diff::Function, init::Function; kwargs...) -> (converged, stalled_update)

Iterate until convergence. In particular, the problem is initiated with `init()`. Then, repeatedly apply `step_diff()` which should do an iteration step, then return the `diff`. `converge` will compare the returned difference with `tol`; if it's smaller, then converge is reached and iteration stops. Otherwise, it will apply `update()`, which should return the size of the update.

The function returns a pair of booleans, the first of which signals whether convergence was reached; the second of which signals whether iteration was aborted because the update was too small.

For an example, see this package's tests which uses this function to solve the finite-firm CES game.

See also [`ConvergeParameters`](@ref), [`update!`](@ref), [`dampen`](@ref), [`v_diff`](@ref).

"""
function converge(update::Function, step_diff::Function, init::Function; kw...)
	converge(update, step_diff, init, ConvergeParameters(; kw...))
end
function converge(update::Function, step_diff::Function, init::Function, p::ConvergeParameters)
	(; history, diff_tol, up_tol, max_iter, msg, verbose) = p
	init()
	reset!(cs)
	diff = one(diff_tol) + diff_tol;
	small_up = false
	conv = false
	
	if verbose && (max_iter > 0)
		meter = ProgressMeter.ProgressThresh(diff_tol; color=:blue, output=stdout, showspeed=true)
	end
	for iter in 1:max_iter
		diff = step_diff()
		!isnothing(history) && push!(history, diff)
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
	update!(main, secondary; dampen = 1.0, v_diff = v_diff, cs = cs, dampen_kw...)

Update `main` according to `secondary`, with a `dampen`ing factor. Useful for iterative algorithms. Once complete, `main` will hold the updated value and `secondary` will hold main's original value (to keep a record of previous iteration). Returns the `v_diff` of the update step (useful for breaking iteration if the step size is too small).

If `dampen` is set at `1.0`, [`dynamic_dampen`](@ref)ing is used.
"""
function update!(main, secondary; cs = cs, dampen = 1.0, v_diff = v_diff, dampen_kw...)
	_dampen = isone(dampen) ? dynamic_dampen!(cs; dampen_kw...)[:dampen] : dampen
	for (i_main, i_sec) in zip(eachindex(main), eachindex(secondary))
		main[i_main], secondary[i_sec] = secondary[i_sec], main[i_main]
		main[i_main] += _dampen*(secondary[i_sec] - main[i_main])
	end
	v_diff(secondary, main)
end
function update(main, secondary; cs = cs, dampen = 1.0, dampen_kw...)
	_dampen = isone(dampen) ? dynamic_dampen!(cs; dampen_kw...)[:dampen] : dampen
	main, secondary = secondary, main
	main += _dampen*(secondary - main)
	main, secondary
end

"""
	v_diff(v1, v2)

Calculate the distance between two vectors as the sum of element-wise absolute difference.
"""
function v_diff(v1, v2)
	maximum(zip(v1, v2)) do (x1, x2)
		abs(x1 - x2)
	end
end

# DYNAMIC DAMPENING
"""
	dynamic_dampen(dampen, last_loosened, last_tightened, reference_diff, last_deviations, penultimate_deviations, history; kw...)

Update the dampening factor based on convergence path. Strategy:
* increase by `tighten` if `history` [`isdiverging`](@ref) or [`isovershooting`](@ref)
* make no change if any of:
    - within the first `grace_period` iterations
    - it has been under `tighten_wait` iterations since the last tightening
    - it has been under `loosen_wait` iterations since the last loosening
    - the current difference is above `scale*reference_diff`
* otherwise, decrease by `loosen`
"""
function dynamic_dampen(dampen, last_loosened, last_tightened, reference_diff, last_deviations, penultimate_deviations, history; loosen = -0.01, tighten = 0.01, min_dampen = 0.0, max_dampen = 0.999, grace_period = 50, tighten_wait = 30, loosen_wait = 40, scale = 0.925, tail = 25, convex_tol = 0.005, overshooting_share = 0.5)
	if (isdiverging(history) || isovershooting(last_deviations, penultimate_deviations; share = overshooting_share)) && (dampen + tighten ≤ max_dampen)
		dampen += tighten
		last_tightened = zero(last_tightened)
		last_loosened += one(last_loosened)
		reference_diff = minimum(@view(history[end-min(2, length(history)-1):end]))
		@debug "tightened from $(dampen - tighten) to $dampen: $(isdiverging(history) ? "diverging" : "overshooting"); reference diff $(last(history))"
	elseif (length(history) < 50) || (last_tightened < tighten_wait) || (last_loosened < loosen_wait) || (last(history) ≥ scale*reference_diff) || (dampen + loosen < min_dampen)
	# || !isconvex(history; tail, tol = convex_tol)
		last_loosened += one(last_loosened)
		last_tightened += one(last_tightened)
	else
		dampen += loosen
		last_loosened = zero(last_loosened)
		last_tightened += one(last_tightened)
		@debug "loosened from $(dampen - loosen) to $dampen"
	end

	(; dampen, last_loosened, last_tightened, reference_diff)
end
"""
	isdiverging(history)

Checks if iteration is on a bad path, but seeing if two of the past 3 iterations have worsened the difference. Useful for dynamically updating the dampening factor.[^1]

See also [`dynamic_dampen`](@ref).

[^1]: If updating is too aggressive, usually one of two things happens. Either the difference blows up (if updating is _too_ aggressive) or the difference oscillates between improving and worsening (if updating is a little too aggressive).
"""
function isdiverging(history)
	(length(history) ≤ 3) && return false
	bads = @inbounds count(0:2) do i
		(history[end-i] - history[end-i-1]) ≥ zero(history[end-i])
	end
	bads ≥ 2
end
"""
	isovershooting(last_deviations, penultimate_deviations; share = 0.5)

Check if the convergence is overshooting each guess, as measured by a `share` of the deviations flipping signs between the last two iterations. Useful for dynamically updating the dampening factor.[^1]

See also [`dynamic_dampen`](@ref).

[^1]: If dampening isn't strong enough, the devationsiations alternate between positive and negative. Convergence wastes times bouncing back and forth, usually slowing down. Increasing the dampening can drastically speed up convergence by preventing this oscillation.
"""
function isovershooting(last_deviations, penultimate_deviations; share = 0.5)
	(isempty(last_deviations) || isempty(penultimate_deviations)) && return false
	(all(isinf, penultimate_deviations) || all(isinf, last_deviations)) && return false
	overshot = count(zip(last_deviations, penultimate_deviations)) do (last, penultimate)
		signbit(last) ⊻ signbit(penultimate)
	end
	overshot ≥ share*length(last_deviations)
end
"""
	isconvex(history; tail = 25, tol = 0.005)

Check if the last `tail` entries in `history` is convex. Useful for dynamically updating the dampening factory.

See also [`dynamic_dampen`](@ref).
"""
function isconvex(history; tail = 25, tol = 0.005)
	(length(history) < tail) && return false
	excerpt = @view history[end-tail+1:end]
	@inbounds all(1:(tail-1)) do i
		(log(excerpt[i]) - log(excerpt[i+1])) ≥ tol
	end
end

"""
	dampen(history; kwargs...)

Given an iteration `history`, return a dampening factor. At the moment, just uses the last value of `history`.

## Keywords
* `slow = 0.95`: the dampening factor for slow updating (when `diff` is above 10).
* `med = 0.75`: the dampening factor for medium updating (when `diff` is above 1).
* `fast = 0.5`: the dampening factor for fast updating.
"""
function dampen(history; slow = 0.95, med = 0.75, fast = 0.5)
	(isnothing(history) || isempty(history) || last(history) > 10) && return slow
	last(history) > 1 && return med
	return fast
end

function bisection(f::Function, interval; converge_kw...)
	init!() = let interval = interval
		!(length(interval) == 2) && error("provided interval doesn't have two points")
		sort!(interval)
	end
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

	converge(update!, diff!, init!; converge_kw...)
end
