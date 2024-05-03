var documenterSearchIndex = {"docs":
[{"location":"#EquilibriumUtilities.jl-Documentation","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"","category":"section"},{"location":"","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"Modules = [EquilibriumUtilities]\nPages   = [\"tilies.jl\", \"lvers.jl\", \"ict.jl]","category":"page"},{"location":"#EquilibriumUtilities.EquilibriumUtilities","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.EquilibriumUtilities","text":"A package of basic utility functions used commonly when computing economic equilibria.\n\nSee converge, newton, normalise!, zero_safe.\n\n\n\n\n\n","category":"module"},{"location":"#EquilibriumUtilities.ConvergeParameters","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.ConvergeParameters","text":"ConvergeParameters{T1 <: Real, T2 <: Integer, T3}\n\nContains optional parameters for the converge solver.\n\nSee also converge.\n\nFields\n\ndiff_tol::T1 = 1e-6: If the value returned by step_diff is less than this tolerance, stop.\nup_tol::T1 = zero(diff_tol): If the value returned by update is less than this tolerance, stop.\nmax_iter::T2 = 200: Maximum number of iterations before giving up.\nmsg::T3 = \"No convergence\": Warning message to display if there isn't convergence within the maximum number of iterations.\nverbose::Bool = false: Print the diff every iteration?\n\n\n\n\n\n","category":"type"},{"location":"#EquilibriumUtilities.DampenParameters","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.DampenParameters","text":"struct DampenParameters{T1 <: Real, T2 <: Integer}\n\nParameters for dynamicdampen!.\n\nFields\n\nloosen::T1 = -0.01: by how much to loosen the dampening factor\ntighten::T1 = 0.01: by how much to tighten the dampening factor\nmin_dampen::T1 = 0.0: the minimum value for the dampen factor\nmax_dampen::T1 = 0.999: the maximum value for the dampen factor\nscale::T1 = 0.925: the scale mentioned in dynamicdampen!\novershooting_share::T1 = 0.5: the share for isovershooting\ngrace_period::T2 = 50: how long to wait between adjustments\ntighten_wait::T2 = 30: how long to wait between tightenings\nloosen_wait::T2 = 40: how long to wait between loosenings\n\n\n\n\n\n","category":"type"},{"location":"#EquilibriumUtilities.DampenState","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.DampenState","text":"struct DampenState{T1 <: Real, T2 <: AbstractVector{T1}}\n\nFields\n\ndampen::T1 = 0.85: the dampening factor\nlast_loosened::T1 = Inf: number of iterations since the dampening factor lowered\nlast_tightened::T1 = Inf: number of iterations since the dampening factor increased\nreference_diff::T1 = Inf: the diff used when checking whether iteration is far along enough to attempt loosening again\nlast_deviations::T2 = Float64[] for storing deviations from the last iteration\npenultimate_deviations::T2 = Float64[] for storing deviations from the penultimate iteration\nhistory::T2 = Float64[] for storing the history of diffs\n\n\n\n\n\n","category":"type"},{"location":"#EquilibriumUtilities.NewtonParameters","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.NewtonParameters","text":"NewtonParameters\n\nContains optional parameters for the newton solver.\n\nSee also newton.\n\nFields\n\nverbose::Bool = false: Print some info every iteration?\nstep_tol::Float64 = 1e-8: If the (absolute value of) step size is smaller than this tolerance, stop.\nf_tol::Float64 = step_tol: If the (absolute value of) function value is smaller than this tolerance, stop.\nmax_iter::Integer = 750: Maximum number of iterations before giving up.\nmsg = \"No Newton convergence\": Warning message to display if there isn't convergence within the maximum number of iterations.\nl::Real = -Inf: Left bound.\nr::Real = Inf: Right bound.\n\n\n\n\n\n","category":"type"},{"location":"#EquilibriumUtilities.WrappedDict","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.WrappedDict","text":"abstract type WrappedDict{T} <: AbstractDict{Symbol, T}\n\nEases creating custom structs that are basically dictionaries.[1] A concrete subtype of WrappedDict must have a field internal_dict::Dict{Symbol, T}, which is the dictionary.\n\nFor an example, see this package's tests with an implementation of a Pollak demand function struct as a WrappedDict.\n\nSee also validate, pretty.\n\n[1]: Essentially implements keys, values, length, iterate, getindex, setindex!, get, and get! by delegating to the necessary field internal_dict.\n\n\n\n\n\n","category":"type"},{"location":"#EquilibriumUtilities.chunk-Tuple{AbstractVector, Integer}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.chunk","text":"chunk(V::AbstractVector, N::Integer)\n\nChunk a vector V into Base.views each of length N. Returns a Tuple of Base.views.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.converge-Tuple{Function, Function, Function}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.converge","text":"converge(update::Function, step_diff::Function, init::Function; kw...) -> (converged, stalled_update)\n\nIterate until convergence. In particular, the problem is initiated with init(). Then, repeatedly:\n\napply step_diff() which should do an iteration step, then return the deviation\napply update(), which should return the size of the update\n\nuntil convergence. The convergence criterion is controlled by keyword arguments; see ConvergeParameters.\n\nThe function returns a pair of booleans, the first of which indicates whether convergence was reached; the second of which indicates whether iteration was aborted because the update was too small.\n\nFor an example, see this package's tests which uses this function to solve the finite-firm CES game.\n\nSee also update!.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.dynamicdampen!-Tuple{EquilibriumUtilities.DampenState, EquilibriumUtilities.DampenParameters}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.dynamicdampen!","text":"dynamicdampen!(cs::DampenState, p::DampenParameters)\n\nUpdate the dampening factor based on convergence path. Strategy:\n\nincrease by cs.tighten if cs.history isdiverging or isovershooting\nmake no change if any of:\nwithin the first p.grace_period iterations\nit has been under p.tighten_wait iterations since the last tightening\nit has been under p.loosen_wait iterations since the last loosening\nthe current difference is above p.scale*cs.reference_diff\notherwise, decrease by loosen\n\nSee also DampenState, DampenParameters.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.isdiverging-Tuple{Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.isdiverging","text":"isdiverging(history)\n\nChecks if iteration is on a bad path, by seeing if two of the past three iterations have worsened the difference. If updating is too aggressive, usually one of two things happens:\n\nthe difference blows up (if updating is too aggressive)\nthe difference oscillates between improving and worsening (if updating is a little too aggressive)\n\nSee also dynamicdampen!.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.isovershooting-Tuple{Any, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.isovershooting","text":"isovershooting(last_deviations, penultimate_deviations; share = 0.5)\n\nCheck if the convergence is overshooting each guess, as measured by a share of the deviations flipping signs between the last two iterations. If dampening isn't strong enough, the deviations alternate between positive and negative, wasting time by bouncing back and forth. Increasing the dampening can drastically speed up convergence by preventing this oscillation.\n\nSee also dynamicdampen!.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.newton","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.newton","text":"Simple Newton method implementation.\n\nSee also NewtonParameters for more keyword argument options.\n\n\n\n\n\n","category":"function"},{"location":"#EquilibriumUtilities.newton-Tuple{Function, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.newton","text":"newton(f_f′, x; kwargs...)\n\nf_f′ should be a function which returns the tuple (f, f′) of function value and derivative. x is the initial guess.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.newton-Tuple{Function, Function, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.newton","text":"newton(f, f′, x; kwargs...)\n\nf should be a function which returns the function value and f′ its derivative. x is the initial guess.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.normalise!","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.normalise!","text":"normalise!(v, factor = zero_safe(first(v)))\n\nNormalise v by some factor. Useful for normalising nominal prices (so the first element is the numeraire) or shares (with factor = sum(v)).\n\nSee also zero_safe.\n\n\n\n\n\n","category":"function"},{"location":"#EquilibriumUtilities.pretty","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.pretty","text":"pretty(s, var_names = keys(s); kw...)\n\nReturn a string version of s for pretty-printing.\n\nKeywords\n\npad = 8: column width\ndigits = 4: number of digits to round to\nspacer = 2: number of spaces in between columns\n\n\n\n\n\n","category":"function"},{"location":"#EquilibriumUtilities.update!-Tuple{Any, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.update!","text":"update!(main, secondary; norm = infnorm_pctdev, dampen = 0.5)\n\nUpdate main according to secondary, with a dampening factor. Useful for iterative algorithms. Once complete, main will hold the updated value and secondary will hold main's original value (to keep a record of previous iteration). Returns the norm of the update step (useful for breaking iteration if the step size is too small).\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.validate-Tuple{Any, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.validate","text":"validate(given, needed)\n\nFor each key in needed, check that it's present in given.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.zero_safe-Tuple{Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.zero_safe","text":"zero_safe(x)\n\nIf x is zero, return one(x). Otherwise, return x. Useful for safely dividing by x.\n\n\n\n\n\n","category":"method"}]
}
