"""
	struct DampenState{T <: AbstractVector{<: Real}}

Fields
===
* `dampenstats::T = [0.85; Inf; Inf; Inf]`: contains, in order,
    * the dampening factor
    * number of iterations since the dampening factor lowered
    * number of iterations since the dampening factor increased
    * the `diff` used when checking whether iteration is far along enough to attempt loosening again
* `last_deviations::T2 = Float64[]` for storing deviations from the last iteration
* `penultimate_deviations::T2 = Float64[]` for storing deviations from the penultimate iteration
* `history::T2 = Float64[]` for storing the history of `diff`s

See also [`dynamicdampen!`](@ref), [`dampenfactor`](@ref), [`dampenfactor!`](@ref).
"""
@kwdef struct DampenState{T <: AbstractVector{<: Real}}
	dampenstats::T = [0.85; Inf; Inf; Inf]
	last_deviations::T = Float64[]
	penultimate_deviations::T = Float64[]
	history::T = Float64[]
end
function initiate_deviations!(ds::DampenState, C::Int)
	resize!(ds.penultimate_deviations, C)
	fill!(ds.penultimate_deviations, Inf)
	resize!(ds.last_deviations, C)
	fill!(ds.last_deviations, Inf)
end
function reset_dampenstats!(ds::DampenState)
	dampenfactor!(ds, 0.85)
	@. ds.dampenstats[2:4] = Inf;
end
const ds = DampenState()
function reset!(ds::DampenState)
	initiate_deviations!(ds, length(ds.last_deviations))
	reset_dampenstats!(ds)
	empty!(ds.history)
end
"""
	dampenfactor(ds::DampenState)

Retrieve the dampen factor from `ds`.

See also [`dampenfactor!`](@ref), [`DampenState`](@ref).
"""
function dampenfactor(ds::DampenState)
	first(ds.dampenstats)
end
"""
	dampenfactor!(ds::DampenState, factor)

Set the dampen factor in `ds`.

See also [`dampenfactor`](@ref), [`DampenState`](@ref).
"""
function dampenfactor!(ds::DampenState, dampen)
	ds.dampenstats[1] = dampen
end

# DYNAMIC DAMPENING
"""
	struct DampenParameters{T1 <: Real, T2 <: Integer}

Parameters for [`dynamicdampen!`](@ref).

Fields
===
* `loosen::T1 = -0.01`: by how much to loosen the dampening factor
* `tighten::T1 = 0.01`: by how much to tighten the dampening factor
* `min_dampen::T1 = 0.0`: the minimum value for the dampen factor
* `max_dampen::T1 = 0.999`: the maximum value for the dampen factor
* `scale::T1 = 0.925`: the scale mentioned in [`dynamicdampen!`](@ref)
* `overshooting_share::T1 = 0.5`: the share for [`isovershooting`](@ref)
* `grace_period::T2 = 50`: how long to wait between adjustments
* `tighten_wait::T2 = 30`: how long to wait between tightenings
* `loosen_wait::T2 = 40`: how long to wait between loosenings
"""
@kwdef struct DampenParameters{T1 <: Real, T2 <: Integer}
	loosen::T1 = -0.01
	tighten::T1 = 0.01
	min_dampen::T1 = 0.0
	max_dampen::T1 = 0.999
	scale::T1 = 0.925
	overshooting_share::T1 = 0.5
	grace_period::T2 = 50
	tighten_wait::T2 = 30
	loosen_wait::T2 = 40
end
"""
	dynamicdampen!(ds::DampenState; kw...)

Update the dampening factor based on convergence path. Strategy:
* increase by `ds.tighten` if `ds.history` [`isdiverging`](@ref) or [`isovershooting`](@ref)
* make no change if any of:
    - within the first `p.grace_period` iterations
    - it has been under `p.tighten_wait` iterations since the last tightening
    - it has been under `p.loosen_wait` iterations since the last loosening
    - the current difference is above `p.scale*ds.reference_diff`
* otherwise, decrease by `loosen`
Fine tune control through keyword arguments; see [`DampenParameters`](@ref).

See also [`DampenState`](@ref), [`DampenParameters`](@ref).
"""
function dynamicdampen!(ds::DampenState; kw...)
	dynamicdampen!(ds, DampenParameters(; kw...))
end
"""
	dynamicdampen!(ds::DampenState, p::DampenParameters)
"""
function dynamicdampen!(ds::DampenState, p::DampenParameters)
	(; dampenstats, last_deviations, penultimate_deviations, history) = ds
	dampen, last_loosened, last_tightened, reference_diff = dampenstats
	(; loosen, tighten, min_dampen, max_dampen, grace_period, tighten_wait, loosen_wait, scale, overshooting_share) = p
	if (isdiverging(history) || isovershooting(last_deviations, penultimate_deviations; share = overshooting_share)) && (dampen + tighten ≤ max_dampen)
		dampenstats[1] += tighten
		dampenstats[2] += one(last_loosened)
		dampenstats[3] = zero(last_tightened)
		dampenstats[4] = minimum(@view(history[end-min(2, length(history)-1):end]))
		@debug "tightened from $(dampen - tighten) to $dampen: $(isdiverging(history) ? "diverging" : "overshooting"); reference diff $(last(history))"
	elseif (length(history) < 50) || (last_tightened < tighten_wait) || (last_loosened < loosen_wait) || (last(history) ≥ scale*reference_diff) || (dampen + loosen < min_dampen)
		dampenstats[2] += one(last_loosened)
		dampenstats[3] += one(last_tightened)
	else
		dampen += loosen
		dampenstats[2] = zero(last_loosened)
		dampenstats[3] += one(last_tightened)
		@debug "loosened from $(dampen - loosen) to $dampen"
	end
	ds
end
"""
	isdiverging(history)

Checks if iteration is on a bad path, by seeing if two of the past three iterations have worsened the difference. If updating is too aggressive, usually one of two things happens:
* the difference blows up (if updating is _too_ aggressive)
* the difference oscillates between improving and worsening (if updating is a little too aggressive)

See also [`dynamicdampen!`](@ref).
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

Check if the convergence is overshooting each guess, as measured by a `share` of the deviations flipping signs between the last two iterations. If dampening isn't strong enough, the deviations alternate between positive and negative, wasting time by bouncing back and forth. Increasing the dampening can drastically speed up convergence by preventing this oscillation.

See also [`dynamicdampen!`](@ref).
"""
function isovershooting(last_deviations, penultimate_deviations; share = 0.5)
	(isempty(last_deviations) || isempty(penultimate_deviations)) && return false
	(all(isinf, penultimate_deviations) || all(isinf, last_deviations)) && return false
	overshot = count(zip(last_deviations, penultimate_deviations)) do (last, penultimate)
		signbit(last) ⊻ signbit(penultimate)
	end
	overshot ≥ share*length(last_deviations)
end

# INTEGRATED METHODS
"""
	converge(update::Function, step_diff::Function, ds::DampenState, p::ConvergeParameters)

Wrapper that `push!`es the calculated `diff` to `ds.history`, so that `ds` is up-to-date.
"""
function converge(update::Function, step_diff::Function, ds::DampenState, p::ConvergeParameters)
	step_diff!() = let step_diff = step_diff, ds = ds
		push!(ds.history, step_diff()) |> last
	end
	converge(update, step_diff!, p)
end
function converge(update::Function, step_diff::Function, ds::DampenState; kw...)
	converge(update, step_diff, ds, ConvergeParameters(; kw...))
end
"""
	update!(main, secondary, ds::DampenState; kw...)

Update using a dynamically-chosen dampening factor. However, if the user supplies the keyword argument `dampen`, it takes precedence over dynamically-chosen dampening factor.

See also [`dynamicdampen!`](@ref), [`DampenState`](@ref).
"""
function update!(main, secondary, ds::DampenState; kw...)
	dampen = dampenfactor(dynamicdampen!(ds))
	update!(main, secondary; dampen, kw...)
end
