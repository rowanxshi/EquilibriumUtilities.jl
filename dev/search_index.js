var documenterSearchIndex = {"docs":
[{"location":"#EquilibriumUtilities.jl-Documentation","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"","category":"section"},{"location":"","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"A small package of utilities commonly used to compute economic models.","category":"page"},{"location":"","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"Exports","category":"page"},{"location":"","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"using EquilibriumUtilities # hide\njoin(names(EquilibriumUtilities), \", \") # hide","category":"page"},{"location":"#Iterative-convergence","page":"EquilibriumUtilities.jl Documentation","title":"Iterative convergence","text":"","category":"section"},{"location":"","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"The generalised iterative routine and its parametres:","category":"page"},{"location":"","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"Modules = [EquilibriumUtilities]\nPages = [\"converge.jl\"]","category":"page"},{"location":"#EquilibriumUtilities.ConvergeParameters","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.ConvergeParameters","text":"ConvergeParameters{T1 <: Real, T2 <: Integer, T3}\n\nContains optional parameters for the converge solver.\n\nSee also converge.\n\nFields\n\ndiff_tol::T1 = 1e-6: If the value returned by step_diff is less than this tolerance, stop.\nup_tol::T1 = zero(diff_tol): If the value returned by update is less than this tolerance, stop.\nmax_iter::T2 = 200: Maximum number of iterations before giving up.\nmsg::T3 = \"No convergence\": Warning message to display if there isn't convergence within the maximum number of iterations.\nverbose::Bool = false: Print the diff every iteration?\n\n\n\n\n\n","category":"type"},{"location":"#EquilibriumUtilities.converge-Tuple{Function, Function}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.converge","text":"converge(update::Function, step_diff::Function; kw...) -> (converged, stalled_update)\n\nIterate until convergence. In particular, repeatedly:\n\napply step_diff() which should do an iteration step, then return the deviation\napply update(), which should return the size of the update\n\nuntil convergence. The convergence criterion is controlled by keyword arguments; see ConvergeParameters.\n\nThe function returns a pair of booleans, the first of which indicates whether convergence was reached; the second of which indicates whether iteration was aborted because the update was too small.\n\nFor an example, see this package's tests which uses this function to solve the finite-firm CES game.\n\nSee also update!.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.converge-Union{Tuple{B}, Tuple{A}, Tuple{A, B, EquilibriumUtilities.ConvergeParameters}} where {A, B}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.converge","text":"converge(update::Function, step_diff::Function, p::ConvergeParameters)\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.infnorm_pctdev-Tuple{Any, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.infnorm_pctdev","text":"infnorm_pctdev(v1, v2)\n\nCalculates the difference between two vectors by finding the maximum absolute percentage deviation.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.update!","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.update!","text":"Optional helper function to use when writing an update() function for converge.\n\nUpdate main according to secondary, with a dampening factor. Useful for iterative algorithms. Once complete, main will hold the updated value and secondary will hold main's original value (to keep a record of previous iteration).\n\nReturns the norm of the update step (useful for breaking iteration if the step size is too small).\n\nSee also converge.\n\n\n\n\n\n","category":"function"},{"location":"#EquilibriumUtilities.update!-Tuple{Any, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.update!","text":"update!(main, secondary; norm = infnorm_pctdev, dampen = 0.5)\n\nSee infnorm_pctdev.\n\n\n\n\n\n","category":"method"},{"location":"#Reactive-dampening","page":"EquilibriumUtilities.jl Documentation","title":"Reactive dampening","text":"","category":"section"},{"location":"","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"Optional helpers, mostly implementing a dampening factor that dynamically adjusts as the equilibrium progresses.","category":"page"},{"location":"","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"The convergence routine above need not use this routine: the user can always define the update() function however they would like. However, this optional functionality allows the solver to choose a dampening factor taking into account the convergence path. It will be provided in the ds::DampenState struct so the user can access it when defining their update() function.","category":"page"},{"location":"","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"Modules = [EquilibriumUtilities]\nPages = [\"dynamicdampen.jl\"]","category":"page"},{"location":"#EquilibriumUtilities.DampenParameters","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.DampenParameters","text":"struct DampenParameters{T1 <: Real, T2 <: Integer}\n\nParameters for dynamicdampen!.\n\nFields\n\nloosen::T1 = -0.01: by how much to loosen the dampening factor\ntighten::T1 = 0.01: by how much to tighten the dampening factor\nmin_dampen::T1 = 0.0: the minimum value for the dampen factor\nmax_dampen::T1 = 0.999: the maximum value for the dampen factor\nscale::T1 = 0.925: the scale mentioned in dynamicdampen!\novershooting_share::T1 = 0.5: the share for isovershooting\ngrace_period::T2 = 50: how long to wait between adjustments\ntighten_wait::T2 = 30: how long to wait between tightenings\nloosen_wait::T2 = 40: how long to wait between loosenings\n\n\n\n\n\n","category":"type"},{"location":"#EquilibriumUtilities.DampenState","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.DampenState","text":"struct DampenState{T <: AbstractVector{<: Real}}\n\nFields\n\ndampenstats::T = [0.85; Inf; Inf; Inf]: contains, in order,\nthe dampening factor\nnumber of iterations since the dampening factor lowered\nnumber of iterations since the dampening factor increased\nthe diff used when checking whether iteration is far along enough to attempt loosening again\nlast_deviations::T2 = Float64[] for storing deviations from the last iteration\npenultimate_deviations::T2 = Float64[] for storing deviations from the penultimate iteration\nhistory::T2 = Float64[] for storing the history of diffs\n\nSee also dynamicdampen!, dampenfactor, dampenfactor!.\n\n\n\n\n\n","category":"type"},{"location":"#EquilibriumUtilities.converge-Tuple{Function, Function, EquilibriumUtilities.DampenState, EquilibriumUtilities.ConvergeParameters}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.converge","text":"converge(update::Function, step_diff::Function, ds::DampenState, p::ConvergeParameters)\n\nWrapper that push!es the calculated diff to ds.history, so that ds is up-to-date.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.dampenfactor!-Tuple{EquilibriumUtilities.DampenState, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.dampenfactor!","text":"dampenfactor!(ds::DampenState, factor)\n\nSet the dampen factor in ds.\n\nSee also dampenfactor, DampenState.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.dampenfactor-Tuple{EquilibriumUtilities.DampenState}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.dampenfactor","text":"dampenfactor(ds::DampenState)\n\nRetrieve the dampen factor from ds.\n\nSee also dampenfactor!, DampenState.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.dynamicdampen!-Tuple{EquilibriumUtilities.DampenState, EquilibriumUtilities.DampenParameters}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.dynamicdampen!","text":"dynamicdampen!(ds::DampenState, p::DampenParameters)\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.dynamicdampen!-Tuple{EquilibriumUtilities.DampenState}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.dynamicdampen!","text":"dynamicdampen!(ds::DampenState; kw...)\n\nUpdate the dampening factor based on convergence path. Strategy:\n\nincrease by ds.tighten if ds.history isdiverging or isovershooting\nmake no change if any of:\nwithin the first p.grace_period iterations\nit has been under p.tighten_wait iterations since the last tightening\nit has been under p.loosen_wait iterations since the last loosening\nthe current difference is above p.scale*ds.reference_diff\notherwise, decrease by loosen\n\nFine tune control through keyword arguments; see DampenParameters.\n\nSee also DampenState, DampenParameters.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.isdiverging-Tuple{Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.isdiverging","text":"isdiverging(history)\n\nChecks if iteration is on a bad path, by seeing if two of the past three iterations have worsened the difference. If updating is too aggressive, usually one of two things happens:\n\nthe difference blows up (if updating is too aggressive)\nthe difference oscillates between improving and worsening (if updating is a little too aggressive)\n\nSee also dynamicdampen!.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.isovershooting-Tuple{Any, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.isovershooting","text":"isovershooting(last_deviations, penultimate_deviations; share = 0.5)\n\nCheck if the convergence is overshooting each guess, as measured by a share of the deviations flipping signs between the last two iterations. If dampening isn't strong enough, the deviations alternate between positive and negative, wasting time by bouncing back and forth. Increasing the dampening can drastically speed up convergence by preventing this oscillation.\n\nSee also dynamicdampen!.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.update!-Tuple{Any, Any, EquilibriumUtilities.DampenState}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.update!","text":"update!(main, secondary, ds::DampenState; kw...)\n\nUpdate using a dynamically-chosen dampening factor. However, if the user supplies the keyword argument dampen, it takes precedence over dynamically-chosen dampening factor.\n\nSee also dynamicdampen!, DampenState.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.update-Tuple{Any, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.update","text":"update(x, dev; dampen = zero(dev), rev = false)\n\nUpdate the value of x according to some deviation dev. If dev is positive, then x will increase; if dev is negative, then x will decrease.\n\nUpdating is concave; that is, x is updated at a slower rate the larger dev gets (to prevent over-updating). The maximum amount that x will get updated is by 2*(1 - dampen) times. \n\nKeywords\n\ndampen = zero(dev): in addition to the concave updating, by how much should the updating dampen?\nrev = false: if true, then a positive dev will decrease x, a negative dev will increase x.\n\n\n\n\n\n","category":"method"},{"location":"#Newton-method","page":"EquilibriumUtilities.jl Documentation","title":"Newton method","text":"","category":"section"},{"location":"","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"Modules = [EquilibriumUtilities]\nPages = [\"newton.jl\"]","category":"page"},{"location":"#EquilibriumUtilities.NewtonParameters","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.NewtonParameters","text":"NewtonParameters\n\nContains optional parameters for the newton solver.\n\nSee also newton.\n\nFields\n\nverbose::Bool = false: Print some info every iteration?\nstep_tol::Float64 = 1e-8: If the (absolute value of) step size is smaller than this tolerance, stop.\nf_tol::Float64 = step_tol: If the (absolute value of) function value is smaller than this tolerance, stop.\nmax_iter::Integer = 750: Maximum number of iterations before giving up.\nmsg = \"No Newton convergence\": Warning message to display if there isn't convergence within the maximum number of iterations.\nl::Real = -Inf: Left bound.\nr::Real = Inf: Right bound.\n\n\n\n\n\n","category":"type"},{"location":"#EquilibriumUtilities.newton","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.newton","text":"Simple Newton method implementation.\n\nSee also NewtonParameters for more keyword argument options.\n\n\n\n\n\n","category":"function"},{"location":"#EquilibriumUtilities.newton-Tuple{Function, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.newton","text":"newton(f_f′, x; kwargs...)\n\nf_f′ should be a function which returns the tuple (f, f′) of function value and derivative. x is the initial guess.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.newton-Tuple{Function, Function, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.newton","text":"newton(f, f′, x; kwargs...)\n\nf should be a function which returns the function value and f′ its derivative. x is the initial guess.\n\n\n\n\n\n","category":"method"},{"location":"#Small-utilities","page":"EquilibriumUtilities.jl Documentation","title":"Small utilities","text":"","category":"section"},{"location":"","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"Modules = [EquilibriumUtilities]\nPages = [\"EquilibriumUtilities.jl\", \"arrayviews.jl\", \"prettyprinting.jl\"]","category":"page"},{"location":"#EquilibriumUtilities.EquilibriumUtilities","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.EquilibriumUtilities","text":"A package of basic utility functions used commonly when computing economic equilibria.\n\nSee converge, newton, normalise!, zero_safe.\n\n\n\n\n\n","category":"module"},{"location":"#EquilibriumUtilities.normalise!","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.normalise!","text":"normalise!(v, factor = zero_safe(first(v)))\n\nNormalise v by some factor. Useful for normalising nominal prices (so the first element is the numeraire) or shares (with factor = sum(v)).\n\nSee also zero_safe.\n\n\n\n\n\n","category":"function"},{"location":"#EquilibriumUtilities.quietly-Tuple{Function}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.quietly","text":"quietly(f::Function)\n\nExecute something without logging it (sends logs to Logging.NullLogger()).\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.zero_safe-Tuple{Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.zero_safe","text":"zero_safe(x)\n\nIf x is zero, return one(x). Otherwise, return x. Useful for safely dividing by x.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.chunk-Tuple{AbstractVector, Integer}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.chunk","text":"chunk(V::AbstractVector, N::Integer)\n\nChunk a vector V into Base.views each of length N. Returns a Tuple of Base.views.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.diag_view-Tuple{AbstractMatrix}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.diag_view","text":"diag_view(mat::AbstractMatrix)\n\nReturns a view of the matrix diagonal, assuming one-based indexing.\n\nSee also offdiag_view, issquare.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.issquare-Tuple{AbstractMatrix}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.issquare","text":"issquare(mat::AbstractMatrix)\n\nIs the matrix square?\n\nSee also diag_view, offdiag_view.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.offdiag_view-Tuple{AbstractMatrix}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.offdiag_view","text":"offdiag_view(mat::AbstractMatrix)\n\nReturns a view of the matrix off-diagonals, assuming one-based indexing.\n\nSee also diag_view, issquare.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.pretty","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.pretty","text":"pretty(s, var_names = keys(s); kw...)\n\nReturn a string version of s for pretty-printing.\n\nKeywords\n\npad = 8: column width\ndigits = 4: number of digits to round to\nspacer = 2: number of spaces in between columns\n\n\n\n\n\n","category":"function"},{"location":"#WrappedDict","page":"EquilibriumUtilities.jl Documentation","title":"WrappedDict","text":"","category":"section"},{"location":"","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"Modules = [EquilibriumUtilities]\nPages = [\"WrappedDict.jl\"]","category":"page"},{"location":"#EquilibriumUtilities.WrappedDict","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.WrappedDict","text":"abstract type WrappedDict{T} <: AbstractDict{Symbol, T}\n\nEases creating custom structs that are basically dictionaries. (Essentially implements keys, values, length, iterate, getindex, setindex!, get, and get! by delegating to the necessary field internal_dict.)\n\nA concrete subtype of WrappedDict must have a field internal_dict::Dict{Symbol, T}, which is the dictionary.\n\nFor an example, see this package's tests with an implementation of a Pollak demand function struct as a WrappedDict.\n\nSee also validate.\n\n\n\n\n\n","category":"type"},{"location":"#EquilibriumUtilities.validate-Tuple{Any, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.validate","text":"validate(given, needed)\n\nFor each key in needed, check that it's present in given.\n\n\n\n\n\n","category":"method"}]
}
