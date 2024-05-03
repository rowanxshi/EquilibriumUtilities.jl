"""
	struct DampenState{T1 <: Real, T2 <: AbstractVector{T1}}

## Fields
* `dampen::T1 = 0.85`: the dampening factor
* `last_loosened::T1 = Inf`: number of iterations since the dampening factor lowered
* `last_tightened::T1 = Inf`: number of iterations since the dampening factor increased
* `reference_diff::T1 = Inf`: the `diff` used when checking whether iteration is far along enough to attempt loosening again
* `last_deviations::T2 = Float64[]` for storing deviations from the last iteration
* `penultimate_deviations::T2 = Float64[]` for storing deviations from the penultimate iteration
* `history::T2 = Float64[]` for storing the history of `diff`s
"""
@kwdef struct DampenState{T1 <: Real, T2 <: AbstractVector{T1}}
	dampen::T1 = 0.85
	last_loosened::T1 = Inf
	last_tightened::T1 = Inf
	reference_diff::T1 = Inf
	last_deviations::T2 = Float64[]
	penultimate_deviations::T2 = Float64[]
	history::T2 = Float64[]
end
function initiate_deviations!(cs::DampenState, C::Int)
	resize!(cs.penultimate_deviations, C)
	fill!(cs.penultimate_deviations, Inf)
	resize!(cs.last_deviations, C)
	fill!(cs.last_deviations, Inf)
end
const ds = DampenState()

# DYNAMIC DAMPENING
"""
	struct DampenParameters{T1 <: Real, T2 <: Integer}

Parameters for [`dynamicdampen!`](@ref).

## Fields
* `loosen::T1 = -0.01`: by how much to loosen the dampening factor
* `tighten::T1 = 0.01`: by how much to tighten the dampening factor
* `min_dampen::T1 = 0.0`: the minimum value for the dampen factor
* `max_dampen::T1 = 0.999`: the maximum value for the dampen factor
* `scale::T1 = 0.925`: by how much to scale the reference deviation
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
function dynamicdampen!(cs::DampenState; kw...)
	dynamicdampen!(cs, DampenParameters(; kw...))
end

"""
	dynamicdampen!(cs::DampenState, p::DampenParameters)

Update the dampening factor based on convergence path. Strategy:
* increase by `cs.tighten` if `cs.history` [`isdiverging`](@ref) or [`isovershooting`](@ref)
* make no change if any of:
    - within the first `p.grace_period` iterations
    - it has been under `p.tighten_wait` iterations since the last tightening
    - it has been under `p.loosen_wait` iterations since the last loosening
    - the current difference is above `p.scale*cs.reference_diff`
* otherwise, decrease by `loosen`

See also [`DampenState`](@ref), [`DampenParameters`](@ref).
"""
function dynamicdampen!(cs::DampenState, p::DampenParameters)
	(; dampen, last_loosened, last_tightened, reference_diff, last_deviations, penultimate_deviations, history) = cs
	(; loosen, tighten, min_dampen, max_dampen, grace_period, tighten_wait, loosen_wait, scale, overshooting_share) = p
	if (isdiverging(history) || isovershooting(last_deviations, penultimate_deviations; share = overshooting_share)) && (dampen + tighten ≤ max_dampen)
		dampen += tighten
		last_tightened = zero(last_tightened)
		last_loosened += one(last_loosened)
		reference_diff = minimum(@view(history[end-min(2, length(history)-1):end]))
		@debug "tightened from $(dampen - tighten) to $dampen: $(isdiverging(history) ? "diverging" : "overshooting"); reference diff $(last(history))"
	elseif (length(history) < 50) || (last_tightened < tighten_wait) || (last_loosened < loosen_wait) || (last(history) ≥ scale*reference_diff) || (dampen + loosen < min_dampen)
		last_loosened += one(last_loosened)
		last_tightened += one(last_tightened)
	else
		dampen += loosen
		last_loosened = zero(last_loosened)
		last_tightened += one(last_tightened)
		@debug "loosened from $(dampen - loosen) to $dampen"
	end
	DampenState(; dampen, last_loosened, last_tightened, reference_diff, last_deviations, penultimate_deviations, history)
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
