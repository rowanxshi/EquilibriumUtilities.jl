var documenterSearchIndex = {"docs":
[{"location":"#EquilibriumUtilities.jl-Documentation","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"","category":"section"},{"location":"","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.jl Documentation","text":"Modules = [EquilibriumUtilities]\nPages   = [\"tilies.jl\", \"lvers.jl\", \"ict.jl]","category":"page"},{"location":"#EquilibriumUtilities.EquilibriumUtilities","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.EquilibriumUtilities","text":"A package of basic utility functions used commonly when computing economic equilibria.\n\nSee converge, normalise!, zero_safe, newton.\n\n\n\n\n\n","category":"module"},{"location":"#EquilibriumUtilities.WrappedDict","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.WrappedDict","text":"abstract type WrappedDict{T} <: AbstractDict{Symbol, T}\n\nEases creating custom structs that are basically dictionaries.[1] A concrete subtype of WrappedDict must have a field internal_dict::Dict{Symbol, T}, which is the dictionary.\n\nFor an example, see this package's tests with an implementation of a Pollak demand function struct as a WrappedDict.\n\nSee also validate, pretty.\n\n[1]: Essentially implements keys, values, length, iterate, getindex, setindex!, get, and get! by delegating to the necessary field internal_dict.\n\n\n\n\n\n","category":"type"},{"location":"#EquilibriumUtilities.converge-Tuple{Function, Function, Function}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.converge","text":"converge(update::Function, step_diff::Function, init::Function; kwargs...)\n\nIterate until convergence. In particular, the problem is initiated with init(). Then, repeatedly apply step_diff() which should do an iteration step, then return the diff. converge will compare the returned difference with tol; if it's smaller, then converge is reached and iteration stops. Otherwise, it will apply update(), which should return the size of the update.\n\nThe function returns a pair of booleans, the first of which signals whether convergence was reached; the second signals whether iteration was aborted because the update was too small.\n\nFor an example, see this package's tests which uses this function to solve the finite-firm CES game.\n\nSee also update!, dampen, v_diff.\n\nKeywords\n\ndiff_tol::Real = 1e-6 : If the value returned by step_diff is less than this tolerance, stop.\nup_tol::Real = zero(diff_tol) : If the value returned by update is less than this tolerance, stop.\nmax_iter::Integer = 200 : Maximum number of iterations before giving up.\nmsg = \"No convergence\" : Warning message to display if there isn't convergence within the maximum number of iterations.\nhistory = nothing : Container in which to store the iteration history of diffs. Useful to check speed of convergence. If nothing is provided, saves no history.\nverbose::Bool = false : Print the diff every iteration?\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.dampen-Tuple{Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.dampen","text":"dampen(history; kwargs...)\n\nGiven an iteration history, return a dampening factor. At the moment, just uses the last value of history.\n\nKeywords\n\nslow = 0.95 : the dampening factor for slow updating (when diff is above 10).\nmed = 0.75 : the dampening factor for medium updating (when diff is above 1).\nfast = 0.5 : the dampening factor for fast updating.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.newton","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.newton","text":"Simple Newton method implementation.\n\nKeywords\n\nstep_tol = 1e-8 : If the (absolute value of) step size is smaller than this tolerance, stop.\nf_tol = step_tol : If the (absolute value of) function value is smaller than this tolerance, stop.\nmax_iter = 750 : Maximum number of iterations before giving up.\nl = -Inf : Left bound.\nr = Inf : Right bound.\nmsg = \"No Newton convergence\" : Warning message to display if there isn't convergence within the maximum number of iterations.\nverbose::Bool = false : Print some info every iteration?\n\n\n\n\n\n","category":"function"},{"location":"#EquilibriumUtilities.newton-Tuple{Function, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.newton","text":"newton(f_f′, x; kwargs...)\n\nf_f′ should be a function which returns the tuple (f, f′) of function value and derivative. x is the initial guess.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.newton-Tuple{Function, Function, Real}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.newton","text":"newton(f, f′, x; kwargs...)\n\nf should be a function which returns the function value and f′ its derivative. x is the initial guess.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.normalise!","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.normalise!","text":"normalise!(v, factor = zero_safe(first(v)))\n\nNormalise v by some factor. Useful for normalising nominal prices (so the first element is the numeraire) or shares (with factor = sum(v)).\n\nSee also zero_safe.\n\n\n\n\n\n","category":"function"},{"location":"#EquilibriumUtilities.pretty","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.pretty","text":"pretty(s::WrappedDict, var_names = keys(s); kw...)\n\nReturn a string version of s for pretty-printing.\n\nKeywords\n\npad = 8: column width\ndigits = 4: number of digits to round to\nspacer = 2: number of spaces in between columns\n\n\n\n\n\n","category":"function"},{"location":"#EquilibriumUtilities.update!-Tuple{Any, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.update!","text":"update!(main, secondary; dampen = 0.75, v_diff = v_diff)\n\nUpdate main according to secondary, with a dampening factor. Useful for iterative algorithms. Once complete, main will hold the updated value and secondary will hold main's original value (to keep a record of previous iteration). Returns the v_diff of the update step (useful for breaking iteration if the step size is too small).\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.v_diff-Tuple{Any, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.v_diff","text":"v_diff(v1, v2)\n\nCalculate the distance between two vectors as the sum of element-wise absolute difference.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.validate-Tuple{Any, Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.validate","text":"validate(given, needed)\n\nFor each key in needed, check that it's present in given.\n\n\n\n\n\n","category":"method"},{"location":"#EquilibriumUtilities.zero_safe-Tuple{Any}","page":"EquilibriumUtilities.jl Documentation","title":"EquilibriumUtilities.zero_safe","text":"zero_safe(x)\n\nIf x is zero, return one(x). Otherwise, return x. Useful for safely dividing by x.\n\n\n\n\n\n","category":"method"}]
}
