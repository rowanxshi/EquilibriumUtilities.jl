const _history = Float64[]
const _dampen_state = Dict(:dampen => 0.85, :last_loosened => Inf, :last_tightened => Inf, :reference_diff => Inf)
function dynamic_dampen!(dampen_state, last_deviations, penultimate_deviations, history; kw...)
	new_state = EquilibriumUtilities.dynamic_dampen(dampen_state[:dampen], dampen_state[:last_loosened], dampen_state[:last_tightened], dampen_state[:reference_diff], last_deviations, penultimate_deviations, history; kw...)
	dampen_state[:dampen] = new_state[1]
	dampen_state[:last_loosened] = new_state[2]
	dampen_state[:last_tightened] = new_state[3]
	dampen_state[:reference_diff] = new_state[4]
	dampen_state
end
dynamic_dampen(last_deviations, penultimate_deviations, history; kw...) = dynamic_dampen!(_dampen_state, last_deviations, penultimate_deviations, history; kw...)
function reset_dampen_state!(dampen_state = _dampen_state)
	dampen_state[:last_loosened] = Inf
	dampen_state[:last_tightened] = Inf
	dampen_state[:reference_diff] = Inf
	dampen_state
end
"""
	const needed_dampen_state = (:dampen, :last_loosened, :last_tightend, :reference_diff)
	
The fields needed for the `dampen_state` `Dict{Symbol}` of a [`ConvergenceState`](@ref).
* `dampen`: the dampening factor
* `last_loosened`: number of iterations since the dampening factor lowered
* `last_tightened`: number of iterations since the dampening factor increased
* `reference_diff`: the `diff` used when checking whether iteration is far along enough to attempt loosening again

See also [`dynamic_dampen`](@ref).
"""
const needed_dampen_state = (:dampen, :last_loosened, :last_tightened, :reference_diff)

"""
	ConvergenceState{D <: Dict{Symbol}, V1 <: AbstractVector, V2 <: AbstractVector}

Keeps track of the state of a [`converge`](@ref)nce routine.

Fields
===
* `dampen_state::Dict{Symbol}` contains the fields [`needed_dampen_state`](@ref)
* `last_deviations::V1` for storing deviations from the last iteration
* `penultimate_deviations::V1` for storing deviations from the penultimate iteration
* `history` for storing the history of `diff`s
"""
struct ConvergenceState{D <: Dict{Symbol}, V1 <: AbstractVector, V2 <: AbstractVector}
	dampen_state::D
	last_deviations::V1
	penultimate_deviations::V1
	history::V2

	function ConvergenceState{D, V1, V2}(dampen_state, last_deviations, penultimate_deviations, history) where {D <: Dict{Symbol}, V1 <: AbstractVector, V2 <: AbstractVector}
		(length(last_deviations) != length(penultimate_deviations)) && error("Deviation vectors need to be of the same length to compare them.")
		validate(dampen_state, needed_dampen_state)
		new{D, V1, V2}(dampen_state, last_deviations, penultimate_deviations, history)
	end
end
function ConvergenceState(dampen_state::D, last_devations::V1, penultimate_deviations::V1, history::V2 = _history) where {D <: Dict, V1 <: AbstractVector, V2 <: AbstractVector}
	ConvergenceState{D, V1, V2}(dampen_state, last_devations, penultimate_deviations, history)
end
function ConvergenceState(last_deviations::V1, penultimate_deviations::V1, history::V2 = _history) where {V1 <: AbstractVector, V2 <: AbstractVector}
	ConvergenceState(_dampen_state, last_deviations, penultimate_deviations, history)
end
function reset!(cs::ConvergenceState)
	reset_dampen_state!(cs.dampen_state)
	empty!(cs.history)
end
const cs = ConvergenceState(_dampen_state, Float64[], Float64[], _history)

function initiate_deviations!(cs::EquilibriumUtilities.ConvergenceState, C::Int)
	resize!(cs.penultimate_deviations, C)
	fill!(cs.penultimate_deviations, Inf)
	resize!(cs.last_deviations, C)
	fill!(cs.last_deviations, Inf)
end

dynamic_dampen!(cs::ConvergenceState; kw...) = dynamic_dampen!(cs.dampen_state, cs.last_deviations, cs.penultimate_deviations, cs.history; kw...)
isdiverging(cs::ConvergenceState) = EquilibriumUtilities.isdiverging(cs.history)
isovershooting(cs::ConvergenceState; share = 0.5) = EquilibriumUtilities.isovershooting(cs.last_deviations, cs.penultimate_deviations; share)
isconvex(cs::ConvergenceState; tail = 25, tol = 0.005) = isconvex(cs.history; tail, tol)

