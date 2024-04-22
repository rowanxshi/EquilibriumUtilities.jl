abstract type ListPrint end
function Base.show(io::IO, m::MIME"text/plain", s::ListPrint)
	for field in fieldnames(typeof(s))
		Printf.@printf "\n%6s: " field
		show(getfield(s, field))
	end
end
"""
	pretty(s::WrappedDict, var_names = keys(s); kw...)

Return a string version of `s` for pretty-printing.
	
# Keywords
* `pad = 8`: column width
* `digits = 4`: number of digits to round to
* `spacer = 2`: number of spaces in between columns
"""
function pretty(s, var_names = keys(s); pad = 8, digits = 4, spacer = 2)
	str = string_header(string.(var_names); pad, spacer)*"\n"
	for line in string_tablify(s, var_names; pad, digits)
		str *= line*"\n"
	end
	str
end
spaced_join(x; pad = 8) = join(rpad.(x, pad), "")
function string_header(var_names; pad = 8, spacer = 2)
	header = spaced_join(var_names; pad)
	underline = repeat("=", pad - spacer)*repeat(" ", spacer)
	header*"\n"*repeat(underline, length(var_names))
end
function string_tablify(s, var_names = keys(s); pad = 8, digits = 4)
	columns = (getproperty(s, v) for v in var_names)
	(spaced_join(round.(line; digits); pad) for line in zip(columns...))
end

